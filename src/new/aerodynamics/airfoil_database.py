from abc import ABC, abstractmethod


class AirfoilDatabase(ABC):
    """Interface for airfoil polar data sources."""

    @classmethod
    @abstractmethod
    def from_folder(cls, folder_path: str) -> "AirfoilDatabase":
        """Construct database from a directory of airfoil data files."""
        ...

    @abstractmethod
    def get_polar(self, name: str) -> dict:
        """Return polar data dict for one airfoil (keyed by Reynolds number)."""
        ...

    @abstractmethod
    def __contains__(self, name: str) -> bool: ...
