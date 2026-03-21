from pathlib import Path

import polars as pl

from .base import ParseResult, PolarParser


class TxtParser(PolarParser):
    """Parser for the .txt airfoil polar format.

    File structure:
        RE,3E6
        cm0,-0.3855
        -10.3,-0.684
        -8.4,-0.472
        ...
        [blank line]
        RE,6E6
        ...
    """

    def file_extension(self) -> str:
        return "txt"

    def parse(self, file_path: Path) -> ParseResult:
        airfoil_name = file_path.stem
        rows: list[tuple[float, float, float, float]] = []

        with open(file_path, "r") as f:
            while True:
                line = f.readline()
                if line == "":
                    break
                line = line.strip()
                if not line:
                    continue

                reynolds = float(line.split(",")[1])

                cm0_line = f.readline().strip()
                cm0 = float(cm0_line.split(",")[1])

                while True:
                    cl_line = f.readline()
                    if cl_line == "" or cl_line.strip() == "":
                        break
                    parts = cl_line.strip().split(",")
                    rows.append((reynolds, float(parts[0]), float(parts[1]), cm0))

        df = pl.DataFrame(
            {
                "reynolds": [r[0] for r in rows],
                "aoa": [r[1] for r in rows],
                "cl": [r[2] for r in rows],
                "cm0": [r[3] for r in rows],
            },
            schema={
                "reynolds": pl.Float64,
                "aoa": pl.Float64,
                "cl": pl.Float64,
                "cm0": pl.Float64,
            },
        )
        return ParseResult(airfoil_name=airfoil_name, df=df)
