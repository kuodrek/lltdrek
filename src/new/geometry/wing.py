from abc import ABC, abstractmethod
from typing import Optional

import numpy as np


class Wing(ABC):
    """Interface for wing geometry implementations.

    Lifecycle:
        1. Instantiate with geometry parameters — concrete __init__ must call generate_mesh().
        2. WingPool calls _apply_flight_condition(flight_condition) during assembly.

    All properties listed below are guaranteed to be available after the
    corresponding lifecycle method has been called.
    """

    # -------------------------------------------------------------------------
    # Core lifecycle methods
    # -------------------------------------------------------------------------

    @abstractmethod
    def generate_mesh(self) -> None:
        """Build the panel mesh and compute all geometric quantities.

        After this call the following properties are available:
        total_span, total_area, MAC, AR, collocation_points, vertice_points,
        u_a, u_n, u_s, cp_lengths, cp_dsl, cp_areas, cp_chords, cp_macs.
        """
        ...

    @abstractmethod
    def _apply_flight_condition(self, flight_condition) -> None:
        """Assign per-panel airfoils and compute Reynolds numbers from the flight condition.

        Called by WingPool during assembly. After this call the following
        properties are available: cp_reynolds, cp_airfoils.
        """
        ...

    # -------------------------------------------------------------------------
    # Wing metadata — available after __init__
    # -------------------------------------------------------------------------

    @property
    @abstractmethod
    def surface_name(self) -> str:
        """Unique identifier for this wing surface."""
        ...

    @property
    @abstractmethod
    def N_panels(self) -> int:
        """Total number of panels along the span."""
        ...

    @property
    @abstractmethod
    def x_pos(self) -> float:
        """Chordwise (x-axis) position offset of the wing root."""
        ...

    # -------------------------------------------------------------------------
    # Computed wing-level metrics — available after generate_mesh()
    # -------------------------------------------------------------------------

    @property
    @abstractmethod
    def total_span(self) -> float:
        """Total semi-span of the wing."""
        ...

    @property
    @abstractmethod
    def total_area(self) -> float:
        """Total reference area (sum of all panel areas)."""
        ...

    @property
    @abstractmethod
    def MAC(self) -> float:
        """Mean aerodynamic chord of the wing."""
        ...

    @property
    @abstractmethod
    def AR(self) -> float:
        """Aspect ratio: (2 * total_span)^2 / (2 * total_area)."""
        ...

    # -------------------------------------------------------------------------
    # Panel mesh outputs — available after generate_mesh()
    # Read-write: WingPool mutates these on mirrored wing copies.
    # -------------------------------------------------------------------------

    @property
    @abstractmethod
    def collocation_points(self) -> np.ndarray:
        """(N_panels, 3) array of panel collocation point coordinates."""
        ...

    @collocation_points.setter
    @abstractmethod
    def collocation_points(self, value: np.ndarray) -> None: ...

    @property
    @abstractmethod
    def vertice_points(self) -> np.ndarray:
        """(N_panels + 1, 3) array of panel vertex coordinates."""
        ...

    @vertice_points.setter
    @abstractmethod
    def vertice_points(self, value: np.ndarray) -> None: ...

    @property
    @abstractmethod
    def u_a(self) -> np.ndarray:
        """(N_panels, 3) unit vectors collinear with the chord."""
        ...

    @u_a.setter
    @abstractmethod
    def u_a(self, value: np.ndarray) -> None: ...

    @property
    @abstractmethod
    def u_n(self) -> np.ndarray:
        """(N_panels, 3) unit vectors normal to the chord plane."""
        ...

    @u_n.setter
    @abstractmethod
    def u_n(self, value: np.ndarray) -> None: ...

    @property
    @abstractmethod
    def u_s(self) -> np.ndarray:
        """(N_panels, 3) unit vectors perpendicular to the airfoil plane (u_a x u_n)."""
        ...

    @u_s.setter
    @abstractmethod
    def u_s(self, value: np.ndarray) -> None: ...

    @property
    @abstractmethod
    def cp_lengths(self) -> np.ndarray:
        """(N_panels, 3) spanwise length vectors for each panel."""
        ...

    @cp_lengths.setter
    @abstractmethod
    def cp_lengths(self, value: np.ndarray) -> None: ...

    @property
    @abstractmethod
    def cp_dsl(self) -> np.ndarray:
        """(N_panels, 3) dimensionless spanwise length vectors."""
        ...

    @cp_dsl.setter
    @abstractmethod
    def cp_dsl(self, value: np.ndarray) -> None: ...

    @property
    @abstractmethod
    def cp_areas(self) -> np.ndarray:
        """(N_panels,) area of each panel."""
        ...

    @property
    @abstractmethod
    def cp_chords(self) -> np.ndarray:
        """(N_panels,) local chord at each collocation point."""
        ...

    @property
    @abstractmethod
    def cp_macs(self) -> np.ndarray:
        """(N_panels,) mean aerodynamic chord of each panel."""
        ...

    # -------------------------------------------------------------------------
    # Airfoil data — available after setup_airfoil_data()
    # -------------------------------------------------------------------------

    @property
    @abstractmethod
    def cp_reynolds(self) -> np.ndarray:
        """(N_panels,) local Reynolds number at each collocation point."""
        ...

    @property
    @abstractmethod
    def cp_airfoils(self) -> list:
        """Per-panel airfoil references with spanwise merge parameters."""
        ...

    # -------------------------------------------------------------------------
    # Cross-wing references
    # -------------------------------------------------------------------------

    @property
    @abstractmethod
    def parent_wing(self) -> Optional[str]:
        """Surface name of the original wing this is a mirror of, or None."""
        ...

    @parent_wing.setter
    @abstractmethod
    def parent_wing(self, value: Optional[str]) -> None: ...
