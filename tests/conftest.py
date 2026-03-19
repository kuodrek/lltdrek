import numpy as np
import pytest

from new.aerodynamics.airfoil_database import AirfoilDatabase
from new.aerodynamics.flight_condition import FlightCondition
from new.geometry.wing_llt import WingLLT
from new.system.wing_pool import WingPool


class StubAirfoilDB(AirfoilDatabase):
    """Minimal concrete AirfoilDatabase for unit tests.

    Returns a wide Reynolds range so WingPool._check_reynolds_bounds never warns.
    Polar data is empty — only suitable for tests that don't exercise polar lookup.
    """

    def get_polar(self, name: str) -> dict:
        return {1_000: {}, 10_000_000: {}}

    @classmethod
    def from_folder(cls, folder_path: str) -> "StubAirfoilDB":
        return cls()

    def __contains__(self, name: str) -> bool:
        return True

    def lookup_cl(self, airfoil: str, reynolds: float, aoa_deg: float) -> float:
        return 0.0

    def get_linear_data(self, airfoil: str, reynolds: float, aoa_min: float = 0.0, aoa_max: float = 8.0) -> dict:
        return {"cl_alpha": 0.1, "cl0": 0.0, "cm0": 0.0, "clmax": 1.0}

    def get_dat(self, name: str) -> np.ndarray:
        return np.array([[0.0, 0.0], [1.0, 0.0]])


@pytest.fixture
def stub_airfoil_db() -> StubAirfoilDB:
    return StubAirfoilDB()


@pytest.fixture
def simple_wing() -> WingLLT:
    return WingLLT(
        spans=[1.0],
        chords=[0.5, 0.5],
        offsets=[0.0, 0.0],
        twist_angles=[0.0, 0.0],
        dihedral_angles=[0.0],
        airfoils=["NACA4412", "NACA4412"],
        surface_name="main_wing",
        N_panels=4,
    )


@pytest.fixture
def simple_flight_condition() -> FlightCondition:
    return FlightCondition(V_inf=10.0, nu=1.5e-5, rho=1.225, angles_of_attack=[0.0, 5.0], h=10.0)


@pytest.fixture
def simple_pool(simple_wing, simple_flight_condition, stub_airfoil_db) -> WingPool:
    return WingPool(wings=[simple_wing], flight_condition=simple_flight_condition, airfoil_db=stub_airfoil_db)
