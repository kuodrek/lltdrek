"""Factory for creating AirfoilDatabase instances."""

from .airfoil_database import AirfoilDatabase

_ALLOWED_BACKENDS = ("polars",)


class AirfoilDatabaseFactory:
    """Factory for instantiating AirfoilDatabase implementations."""

    @staticmethod
    def create(folder_path: str, backend: str = "polars", format: str = "txt") -> AirfoilDatabase:
        """Create an AirfoilDatabase instance.

        Args:
            folder_path: Path to directory containing airfoil files.
            backend: Backend implementation to use ('polars'). Defaults to 'polars'.
            format: Polar file format ('txt', 'csv'). Defaults to 'txt'.

        Returns:
            An AirfoilDatabase instance.

        Raises:
            ValueError: If backend is not supported.
        """
        if backend not in _ALLOWED_BACKENDS:
            raise ValueError(f"Invalid backend '{backend}'. Choose from {_ALLOWED_BACKENDS}.")

        if backend == "polars":
            from .polars_airfoil_database import PolarsAirfoilDatabase

            return PolarsAirfoilDatabase.from_folder(folder_path, format=format)
