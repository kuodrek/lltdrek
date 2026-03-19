import copy
import logging
from typing import Dict, List, Optional, Sequence

import numpy as np

from new.aerodynamics.airfoil_database import AirfoilDatabase
from new.aerodynamics.flight_condition import FlightCondition
from new.aerodynamics.velocity import LoopsVelocityCalculator, NumpyVelocityCalculator, VelocityCalculator
from new.geometry.wing import Wing

_ALLOWED_VELOCITY = ("loops", "numpy")
logger = logging.getLogger("lltdrek.system.WingPool")


class WingPool:
    """Assembles a collection of Wing surfaces into a solvable system.

    Pre-computes freestream and induced velocity distributions in __init__ so
    that SimulationRunner implementations can retrieve them during iteration
    without recomputing.

    Parameters
    ----------
    velocity : "loops" | "numpy"
        Backend used to pre-compute induced velocity matrices in
        _build_system_induced_velocities. "loops" is the validated baseline;
        "numpy" is the vectorized implementation (stub until implemented).
    S_ref : float, optional
        Reference area. Defaults to sum of all pool wing areas.
    c_ref : float, optional
        Reference chord. Defaults to MAC of the first wing.
    moment_ref : sequence of float
        [x, y, z] reference point for moment calculations.
    """

    def __init__(
        self,
        wings: List[Wing],
        flight_condition: FlightCondition,
        airfoil_db: AirfoilDatabase,
        velocity: str = "loops",
        S_ref: Optional[float] = None,
        c_ref: Optional[float] = None,
        moment_ref: Sequence[float] = (0, 0, 0),
    ) -> None:
        if velocity not in _ALLOWED_VELOCITY:
            raise ValueError(f"Invalid velocity '{velocity}'. Choose from {_ALLOWED_VELOCITY}.")

        self.wings = wings
        self.flight_condition = flight_condition
        self.airfoil_db = airfoil_db
        self._moment_ref = np.array(moment_ref)
        self.velocity_calculator: VelocityCalculator = (
            LoopsVelocityCalculator() if velocity == "loops" else NumpyVelocityCalculator()
        )

        for wing in self.wings:
            wing._apply_flight_condition(flight_condition)
            self._check_reynolds_bounds(wing)

        self.pool = self._build_pool()
        self.total_panels = sum(w.N_panels for w in self.pool)

        self.S_ref = S_ref if S_ref is not None else sum(w.total_area for w in self.pool)
        self.c_ref = c_ref if c_ref is not None else self.wings[0].MAC

        self.system_moment_ref = self._build_system_moment_ref()
        self.system_freestream_velocities = self._build_system_freestream_velocities()
        self.system_induced_velocities = self._build_system_induced_velocities()

    @property
    def moment_ref(self) -> np.ndarray:
        return self._moment_ref

    @moment_ref.setter
    def moment_ref(self, value: Sequence) -> None:
        self._moment_ref = np.array(value)
        self.system_moment_ref = self._build_system_moment_ref()

    def _build_pool(self) -> List[Wing]:
        """Create the full pool: each original wing followed by its mirrored copy."""
        pool = []
        for wing in self.wings:
            mirrored = copy.deepcopy(wing)
            mirrored.surface_name = wing.surface_name + "_mirrored"
            mirrored.parent_wing = wing.surface_name

            mirrored.u_a[:, 1] = -1 * mirrored.u_a[:, 1]
            mirrored.u_n[:, 1] = -1 * mirrored.u_n[:, 1]
            mirrored.u_s[:, 0] = -1 * mirrored.u_s[:, 0]
            mirrored.u_s[:, 2] = -1 * mirrored.u_s[:, 2]
            mirrored.collocation_points[:, 1] = -1 * mirrored.collocation_points[:, 1]
            mirrored.vertice_points[:, 1] = -1 * mirrored.vertice_points[:, 1]
            mirrored.cp_lengths[:, 1] = -1 * mirrored.cp_lengths[:, 1]
            mirrored.cp_dsl[:, 0] = -1 * mirrored.cp_dsl[:, 0]
            mirrored.cp_dsl[:, 2] = -1 * mirrored.cp_dsl[:, 2]

            pool.append(wing)
            pool.append(mirrored)
        return pool

    def _build_system_moment_ref(self) -> Dict[str, np.ndarray]:
        """Position vectors from each collocation point to the moment reference."""
        return {wing.surface_name: wing.collocation_points - self._moment_ref for wing in self.pool}

    def _get_angular_velocities(self) -> Dict[str, np.ndarray]:
        """Linear velocity at each panel due to angular rates: r x omega."""
        omega = self.flight_condition.angular_velocity
        return {name: np.cross(ref, omega) for name, ref in self.system_moment_ref.items()}

    def _build_system_freestream_velocities(self) -> Dict[float, Dict[str, np.ndarray]]:
        """Freestream velocity distribution per angle of attack per wing.

        Result structure: {alpha: {surface_name: (N_panels, 3)}}
        """
        angular_velocities = self._get_angular_velocities()
        system_freestream_velocities = {}
        for i, alpha in np.ndenumerate(self.flight_condition.angles_of_attack):
            v_inf = self.flight_condition.v_inf_list[i] * self.flight_condition.V_inf
            wing_freestream_velocities = {}
            for wing in self.pool:
                wing_freestream_velocities[wing.surface_name] = (
                    np.tile(v_inf, (wing.N_panels, 1)) + angular_velocities[wing.surface_name]
                )
            system_freestream_velocities[alpha] = wing_freestream_velocities
        return system_freestream_velocities

    def _build_system_induced_velocities(self) -> Dict[float, Dict[str, Dict[str, np.ndarray]]]:
        """Induced velocity matrices for all (wing_i, wing_j) pairs per alpha.

        Result structure: {alpha: {wing_i_name: {wing_j_name: (N_cp_i, N_panels_j, 3)}}}
        """
        system_induced_velocities = {}
        for i, alpha in np.ndenumerate(self.flight_condition.angles_of_attack):
            freestream = self.system_freestream_velocities[alpha]
            wing_induced_velocities = {}
            for wing_cp in self.pool:
                wing_induced_velocities[wing_cp.surface_name] = {}
                for wing_vp in self.pool:
                    velocity_distribution = self.velocity_calculator.get_induced_velocity_distribution(
                        wing_cp.collocation_points,
                        wing_cp.cp_macs,
                        wing_vp.vertice_points,
                        freestream[wing_cp.surface_name],
                        is_mirrored=(wing_vp.parent_wing is not None),
                        ground_effect=self.flight_condition.ground_effect_check,
                        h=self.flight_condition.h,
                    )
                    wing_induced_velocities[wing_cp.surface_name][wing_vp.surface_name] = velocity_distribution
            system_induced_velocities[alpha] = wing_induced_velocities
        return system_induced_velocities

    def map_solution(self, G: np.ndarray) -> Dict[str, np.ndarray]:
        """Split flat G array (length = total_panels) into per-wing arrays."""
        G_dict = {}
        offset = 0
        for wing in self.pool:
            G_dict[wing.surface_name] = G[offset : offset + wing.N_panels]
            offset += wing.N_panels
        return G_dict

    def calculate_total_velocity(self, alpha: float, G_dict: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """Total velocity per panel: v_freestream + sum_j(v_induced_ij * G_j).

        Result structure: {surface_name: (N_panels, 3)}
        """
        total_velocities = {}
        for wing_i in self.pool:
            v_total = self.system_freestream_velocities[alpha][wing_i.surface_name].copy()
            for wing_j in self.pool:
                v_ij = self.system_induced_velocities[alpha][wing_i.surface_name][wing_j.surface_name]
                v_total += np.einsum("ijk,j->ik", v_ij, G_dict[wing_j.surface_name])
            total_velocities[wing_i.surface_name] = v_total
        return total_velocities

    def calculate_aoa_eff(self, total_velocity_dict: Dict[str, np.ndarray]) -> Dict[str, np.ndarray]:
        """Effective angle of attack per panel: arctan(v.u_n / v.u_a).

        Computed for original wings only; mirrored wings copy the same result.
        Result structure: {surface_name: (N_panels,)}
        """
        system_aoa_eff = {}
        for wing in self.wings:
            v_distr = total_velocity_dict[wing.surface_name]
            aoa_eff = np.arctan(np.einsum("ij,ij->i", v_distr, wing.u_n) / np.einsum("ij,ij->i", v_distr, wing.u_a))
            system_aoa_eff[wing.surface_name] = aoa_eff
            system_aoa_eff[wing.surface_name + "_mirrored"] = aoa_eff
        return system_aoa_eff

    def _check_reynolds_bounds(self, wing: Wing) -> None:
        min_cp_re = float(min(wing.cp_reynolds))
        max_cp_re = float(max(wing.cp_reynolds))

        airfoil_names = set(name for _, names in wing.cp_airfoils for name in names)
        reynolds_list = [float(re) for name in airfoil_names for re in self.airfoil_db.get_polar(name)]

        if min_cp_re < min(reynolds_list) or max_cp_re > max(reynolds_list):
            logger.warning(
                f"[{wing.surface_name}] Detected reynolds out of airfoil data bounds. Results may be inaccurate."
            )
