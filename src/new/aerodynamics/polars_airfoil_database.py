from pathlib import Path

import numpy as np
import polars as pl

from .airfoil_database import AirfoilDatabase
from .parsers.base import PolarParser
from .parsers.csv_parser import CsvParser
from .parsers.txt_parser import TxtParser

_PARSER_REGISTRY: dict[str, type[PolarParser]] = {
    "txt": TxtParser,
    "csv": CsvParser,
}


class PolarsAirfoilDatabase(AirfoilDatabase):
    """Concrete AirfoilDatabase backed by a Polars DataFrame.

    Internal schema: [airfoil: Utf8, reynolds: Float64, aoa: Float64, cl: Float64, cm0: Float64]
    """

    def __init__(
        self,
        df: pl.DataFrame,
        geometry: dict[str, np.ndarray],
    ):
        self._df = df.sort(["airfoil", "reynolds", "aoa"])
        self._geometry = geometry
        self._airfoil_names: frozenset[str] = frozenset(df["airfoil"].unique().to_list())

        # Hot-path caches populated on first access
        self._reynolds_cache: dict[str, np.ndarray] = {}
        self._cl_data_cache: dict[tuple[str, float], np.ndarray] = {}
        self._linear_cache: dict[tuple, dict] = {}

    @classmethod
    def from_folder(cls, folder_path: str, **kwargs) -> "PolarsAirfoilDatabase":
        """Load all airfoil files from a directory.

        Args:
            folder_path: Path to directory containing airfoil files.
            format: Polar file format ('txt' or 'csv'). Defaults to 'txt'.
        """
        fmt = kwargs.get("format", "txt")
        if fmt not in _PARSER_REGISTRY:
            raise KeyError(f"Unknown format '{fmt}'. Available: {list(_PARSER_REGISTRY)}")

        parser = _PARSER_REGISTRY[fmt]()
        folder = Path(folder_path)

        all_dfs: list[pl.DataFrame] = []
        for polar_file in sorted(folder.glob(f"*.{parser.file_extension()}")):
            result = parser.parse(polar_file)
            df_with_name = result.df.with_columns(pl.lit(result.airfoil_name).alias("airfoil"))
            all_dfs.append(df_with_name)

        if all_dfs:
            combined = pl.concat(all_dfs).select(["airfoil", "reynolds", "aoa", "cl", "cm0"])
        else:
            combined = pl.DataFrame(
                schema={
                    "airfoil": pl.Utf8,
                    "reynolds": pl.Float64,
                    "aoa": pl.Float64,
                    "cl": pl.Float64,
                    "cm0": pl.Float64,
                }
            )

        geometry: dict[str, np.ndarray] = {}
        for dat_file in sorted(folder.glob("*.dat")):
            geometry[dat_file.stem] = np.loadtxt(dat_file, skiprows=1)

        return cls(combined, geometry)

    # ------------------------------------------------------------------
    # AirfoilDatabase interface
    # ------------------------------------------------------------------

    def __contains__(self, name: object) -> bool:
        return name in self._airfoil_names

    def get_polar(self, name: str) -> dict:
        """Return {reynolds_float: {"cm0": ..., "clmax": ...}} for WingPool compatibility."""
        if name not in self:
            raise KeyError(f"Airfoil '{name}' not found in database.")
        subset = self._df.filter(pl.col("airfoil") == name)
        grouped = subset.group_by("reynolds").agg(
            [
                pl.col("cm0").first().alias("cm0"),
                pl.col("cl").max().alias("clmax"),
            ]
        )
        return {row["reynolds"]: {"cm0": row["cm0"], "clmax": row["clmax"]} for row in grouped.iter_rows(named=True)}

    def lookup_cl(self, airfoil: str, reynolds: float, aoa_deg: float) -> float:
        re_arr = self._get_reynolds_list(airfoil)

        if reynolds <= re_arr[0]:
            return float(self._interpolate_aoa(airfoil, re_arr[0], aoa_deg))
        if reynolds >= re_arr[-1]:
            return float(self._interpolate_aoa(airfoil, re_arr[-1], aoa_deg))

        idx = int(np.searchsorted(re_arr, reynolds, side="right")) - 1
        re_lo, re_hi = re_arr[idx], re_arr[idx + 1]
        cl_lo = self._interpolate_aoa(airfoil, re_lo, aoa_deg)
        cl_hi = self._interpolate_aoa(airfoil, re_hi, aoa_deg)
        frac = (reynolds - re_lo) / (re_hi - re_lo)
        return float(cl_lo + frac * (cl_hi - cl_lo))

    def get_linear_data(self, airfoil: str, reynolds: float, aoa_min: float = 0.0, aoa_max: float = 8.0) -> dict:
        cache_key = (airfoil, reynolds, aoa_min, aoa_max)
        if cache_key in self._linear_cache:
            return self._linear_cache[cache_key]

        result = {k: float(v) for k, v in self._compute_linear_data(airfoil, reynolds, aoa_min, aoa_max).items()}
        self._linear_cache[cache_key] = result
        return result

    def get_dat(self, name: str) -> np.ndarray:
        if name not in self._geometry:
            raise KeyError(f"No geometry data for airfoil '{name}'.")
        return self._geometry[name]

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    def _get_reynolds_list(self, airfoil: str) -> np.ndarray:
        if airfoil not in self._reynolds_cache:
            re_values = (
                self._df.filter(pl.col("airfoil") == airfoil)
                .select("reynolds")
                .unique()
                .sort("reynolds")["reynolds"]
                .to_numpy()
            )
            self._reynolds_cache[airfoil] = re_values
        return self._reynolds_cache[airfoil]

    def _get_cl_data(self, airfoil: str, reynolds: float) -> np.ndarray:
        key = (airfoil, reynolds)
        if key not in self._cl_data_cache:
            subset = (
                self._df.filter((pl.col("airfoil") == airfoil) & (pl.col("reynolds") == reynolds))
                .select(["aoa", "cl"])
                .sort("aoa")
            )
            self._cl_data_cache[key] = subset.to_numpy()
        return self._cl_data_cache[key]

    def _interpolate_aoa(self, airfoil: str, reynolds: float, aoa_deg: float) -> np.float64:
        cl_data = self._get_cl_data(airfoil, reynolds)
        return np.interp(aoa_deg, cl_data[:, 0], cl_data[:, 1])

    def _compute_linear_data_at_reynolds(self, airfoil: str, reynolds: float, aoa_min: float, aoa_max: float) -> dict:
        subset = self._df.filter(
            (pl.col("airfoil") == airfoil)
            & (pl.col("reynolds") == reynolds)
            & (pl.col("aoa") >= aoa_min)
            & (pl.col("aoa") <= aoa_max)
        )
        aoa_arr = subset["aoa"].to_numpy()
        cl_arr = subset["cl"].to_numpy()
        coefs = np.polyfit(aoa_arr, cl_arr, 1)

        cm0 = self._df.filter((pl.col("airfoil") == airfoil) & (pl.col("reynolds") == reynolds))["cm0"].first()
        clmax = self._df.filter((pl.col("airfoil") == airfoil) & (pl.col("reynolds") == reynolds))["cl"].max()
        return {"cl_alpha": coefs[0], "cl0": coefs[1], "cm0": cm0, "clmax": clmax}

    def _compute_linear_data(self, airfoil: str, reynolds: float, aoa_min: float, aoa_max: float) -> dict:
        re_arr = self._get_reynolds_list(airfoil)

        if reynolds <= re_arr[0]:
            return self._compute_linear_data_at_reynolds(airfoil, re_arr[0], aoa_min, aoa_max)
        if reynolds >= re_arr[-1]:
            return self._compute_linear_data_at_reynolds(airfoil, re_arr[-1], aoa_min, aoa_max)

        idx = int(np.searchsorted(re_arr, reynolds, side="right")) - 1
        re_lo, re_hi = re_arr[idx], re_arr[idx + 1]
        d_lo = self._compute_linear_data_at_reynolds(airfoil, re_lo, aoa_min, aoa_max)
        d_hi = self._compute_linear_data_at_reynolds(airfoil, re_hi, aoa_min, aoa_max)

        frac = (reynolds - re_lo) / (re_hi - re_lo)
        return {key: d_lo[key] + frac * (d_hi[key] - d_lo[key]) for key in ("cl_alpha", "cl0", "cm0", "clmax")}
