from abc import ABC, abstractmethod

import numpy as np


class AirfoilDatabase(ABC):
    """Interface for airfoil polar data sources."""

    @classmethod
    @abstractmethod
    def from_folder(cls, folder_path: str, **kwargs) -> "AirfoilDatabase":
        """Construct database from a directory of airfoil data files.

        Args:
            folder_path: Path to directory containing airfoil files.
            format: Polar file format to parse ('txt', 'csv'). Defaults to 'txt'.
        """
        ...

    @abstractmethod
    def get_polar(self, name: str) -> dict:
        """Return polar data dict for one airfoil keyed by Reynolds number (float).

        Compatible with WingPool._check_reynolds_bounds which iterates the keys.
        """
        ...

    @abstractmethod
    def __contains__(self, name: str) -> bool: ...

    @abstractmethod
    def lookup_cl(self, airfoil: str, reynolds: float, aoa_deg: float) -> float:
        """Interpolate Cl for a given (airfoil, reynolds, aoa) point.

        Uses bilinear interpolation: bracket Reynolds, then interpolate AOA
        within each bracketing Reynolds, then interpolate across Reynolds.
        Out-of-bounds Reynolds or AOA values are clamped to the data range.
        """
        ...

    @abstractmethod
    def get_linear_data(self, airfoil: str, reynolds: float, aoa_min: float = 0.0, aoa_max: float = 8.0) -> dict:
        """Return linear aerodynamic coefficients interpolated at the given Reynolds.

        Fits a line to Cl vs AOA data within [aoa_min, aoa_max] to derive cl_alpha and cl0.
        Results are cached after first computation.

        Returns:
            {"cl_alpha": float, "cl0": float, "cm0": float, "clmax": float}
        """
        ...

    @abstractmethod
    def get_dat(self, name: str) -> np.ndarray:
        """Return (N, 2) airfoil geometry coordinate array from a .dat file.

        Raises:
            KeyError: if no geometry data is available for the given airfoil name.
        """
        ...
