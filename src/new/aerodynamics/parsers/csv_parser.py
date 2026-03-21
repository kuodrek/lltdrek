from pathlib import Path

from .base import ParseResult, PolarParser


class CsvParser(PolarParser):
    """Parser for CSV airfoil polar data.

    Expected columns: reynolds, aoa, cl, cm0
    """

    def file_extension(self) -> str:
        return "csv"

    def parse(self, file_path: Path) -> ParseResult:
        raise NotImplementedError("CsvParser is not yet implemented")
