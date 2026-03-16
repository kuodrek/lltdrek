from typing import List

from new.aerodynamics.airfoil_database import AirfoilDatabase
from new.aerodynamics.flight_condition import FlightCondition
from new.aerodynamics.velocity import LoopsVelocityCalculator, NumpyVelocityCalculator, VelocityCalculator
from new.geometry.wing import Wing

_ALLOWED_VELOCITY = ("loops", "numpy")


class WingPool:
    """Assembles a collection of Wing surfaces into a solvable system.

    Lifecycle:
        1. __init__ receives wings + airfoil database, calls _apply_flight_condition() on each wing.

    Parameters
    ----------
    velocity : "loops" | "numpy"
        Backend used to pre-compute induced velocity matrices in
        _build_system_induced_velocities. "loops" is the validated baseline;
        "numpy" is the vectorized implementation (stub until implemented).
    """

    def __init__(
        self,
        wings: List[Wing],
        flight_condition: FlightCondition,
        airfoil_db: AirfoilDatabase,
        velocity: str = "loops",
    ) -> None:
        if velocity not in _ALLOWED_VELOCITY:
            raise ValueError(f"Invalid velocity '{velocity}'. Choose from {_ALLOWED_VELOCITY}.")

        self.wings = wings
        self.flight_condition = flight_condition
        self.airfoil_db = airfoil_db
        self.velocity_calculator: VelocityCalculator = (
            LoopsVelocityCalculator() if velocity == "loops" else NumpyVelocityCalculator()
        )

        for wing in self.wings:
            wing._apply_flight_condition(flight_condition)
            self._check_reynolds_bounds(wing)

    def _check_reynolds_bounds(self, wing: Wing) -> None:
        min_cp_re = float(min(wing.cp_reynolds))
        max_cp_re = float(max(wing.cp_reynolds))

        airfoil_names = set(name for _, names in wing.cp_airfoils for name in names)
        reynolds_list = [float(re) for name in airfoil_names for re in self.airfoil_db.get_polar(name)]

        if min_cp_re < min(reynolds_list) or max_cp_re > max(reynolds_list):
            print(
                f"Warning: [{wing.surface_name}] Detected reynolds out of airfoil data bounds."
                " Results may be inaccurate."
            )
