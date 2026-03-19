from abc import ABC, abstractmethod
from dataclasses import dataclass
from pathlib import Path

import polars as pl


@dataclass
class ParseResult:
    """Result of parsing a single airfoil polar file.

    Schema of df: [reynolds: Float64, aoa: Float64, cl: Float64, cm0: Float64]
    cm0 is repeated for every row in a Reynolds block.
    """

    airfoil_name: str
    df: pl.DataFrame


class PolarParser(ABC):
    """Strategy for parsing one airfoil polar file into a Polars DataFrame."""

    @abstractmethod
    def parse(self, file_path: Path) -> ParseResult:
        """Parse a single file and return structured data."""
        ...

    @abstractmethod
    def file_extension(self) -> str:
        """Return the file extension this parser handles (e.g. 'txt')."""
        ...
