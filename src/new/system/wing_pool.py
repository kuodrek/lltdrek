from typing import List

from new.aerodynamics.airfoil_database import AirfoilDatabase
from new.geometry.wing import Wing
from new.aerodynamics.flight_condition import FlightCondition


class WingPool:
    """Assembles a collection of Wing surfaces into a solvable system.

    Lifecycle:
        1. __init__ receives wings + airfoil database, calls generate_mesh() on each wing.
    """

    def __init__(
        self,
        wings: List[Wing],
        flight_condition: FlightCondition,
        airfoil_db: AirfoilDatabase,
    ) -> None:
        self.wings = wings
        self.flight_condition = flight_condition
        self.airfoil_db = airfoil_db

        for wing in self.wings:
            wing.generate_mesh()
            wing._setup_airfoil_data(flight_condition, airfoil_db)
