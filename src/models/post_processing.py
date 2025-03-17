from copy import deepcopy
from dataclasses import dataclass
from typing import List, Optional, Union

import numpy as np

from src.models.simulation import SimulationResult
from src.models.wingpool import WingPool
from src.utils.lookup import get_linear_data_and_clmax

class AerodynamicCenter:
    def __init__(self, x_ac: float, Cm_ac: float, Cm_alpha: float, n_iter: int):
        self.x_ac = x_ac
        self.Cm_ac = Cm_ac
        self.Cm_alpha = Cm_alpha
        self.n_iter = n_iter

    def __repr__(self):
        return f"AerodynamicCenter(x_ac={self.x_ac}, Cm_ac={self.Cm_ac}, Cm_alpha={self.Cm_alpha}, n_iter={self.n_iter})"

class ForceCoefficients:
    def __init__(self, CD: float, CY: float, CL: float):
        self.CD = CD
        self.CY = CY
        self.CL = CL

    def __repr__(self):
        return f"ForceCoefficients(CD={self.CD}, CY={self.CY}, CL={self.CL})"


class MomentCoefficients:
    def __init__(self, Cl: float, Cm: float, Cn: float):
        self.Cl = Cl
        self.Cm = Cm
        self.Cn = Cn

    def __repr__(self):
        return f"MomentCoefficients(Cl={self.Cl}, Cm={self.Cm}, Cn={self.Cn})"


class Coefficients:
    def __init__(self, force_coefficients, moment_coefficients, cl_distribution):
        self.forces = ForceCoefficients(*force_coefficients)
        self.moments = MomentCoefficients(*moment_coefficients)
        self.cl_distribution = cl_distribution

    def __repr__(self):
        return f"Coefficients(forces={self.forces}, moments={self.moments}, cl_distribution={self.cl_distribution})"


@dataclass
class ProcessedSimulationResults:
    simulation_result: SimulationResult
    global_coefficients: Coefficients
    surface_coefficients: dict[str, Coefficients]


class PostProcessing:
    @classmethod
    def _build_nan_results(cls, simulation_result: SimulationResult, wing_pool: WingPool):
        nan_coefficients = lambda N: Coefficients(
            force_coefficients=np.array((np.nan, np.nan, np.nan)),
            moment_coefficients=np.array((np.nan, np.nan, np.nan)),
            cl_distribution=np.full(N, np.nan),
        )
        global_coefficients = {}
        surface_coefficients = {}
        for wing in wing_pool.pool:
            if "_mirrored" in wing.surface_name:
                continue
            surface_coefficients[wing.surface_name] = nan_coefficients(wing.N_panels)
        global_coefficients = nan_coefficients(sum([wing.N_panels for wing in wing_pool.pool]) // 2)
        return ProcessedSimulationResults(simulation_result, global_coefficients, surface_coefficients)

    @classmethod
    def get_coefficients(
        cls,
        wing_pool: WingPool,
        simulation_results: List[SimulationResult],
        S_ref: Optional[float] = None,
        c_ref: Optional[float] = None,
    ) -> List[ProcessedSimulationResults]:
        if not S_ref:
            S_ref = wing_pool.S_ref
        if not c_ref:
            c_ref = wing_pool.c_ref

        processed_simulation_results = []
        for result in simulation_results:
            if not result.convergence_check:
                processed_simulation_results.append(cls._build_nan_results(result, wing_pool))
                continue

            G_dict = result.G_solution

            surface_coefficients: dict[str, Coefficients] = {}
            for wing_i in wing_pool.pool:
                CF = np.zeros(3)
                CM = np.zeros(3)
                cl_distribution = np.zeros(wing_i.N_panels)

                wing_freestreams_velocities = wing_pool.system_freestream_velocities[result.alpha][wing_i.surface_name]
                ref_points_distr = wing_pool.system_moment_ref[wing_i.surface_name]
                G_i = G_dict[wing_i.surface_name]
                for i, _ in enumerate(wing_i.collocation_points):
                    aux_cf = np.zeros(3)

                    airfoil_coefficients = get_linear_data_and_clmax(
                        wing_i.cp_airfoils[i], wing_i.cp_reynolds[i], wing_i.airfoil_data, show_logs=False
                    )
                    cm_i = airfoil_coefficients["cm0"]
                    cl_distribution[i] = 2 * G_i[i]
                    for wing_j in wing_pool.pool:
                        G_j = G_dict[wing_j.surface_name]
                        v_ij_distr = wing_pool.system_induced_velocities[result.alpha][wing_i.surface_name][
                            wing_j.surface_name
                        ]
                        for j, _ in enumerate(wing_j.collocation_points):
                            v_ij = v_ij_distr[i][j]
                            aux_cf += G_j[j] * v_ij
                    aux_cf = (aux_cf + wing_freestreams_velocities[i]) * G_i[i]
                    CF += 2 * np.cross(aux_cf, wing_i.cp_dsl[i]) * wing_i.cp_areas[i] / S_ref
                    CM += (
                        (
                            2 * np.cross(ref_points_distr[i], np.cross(aux_cf, wing_i.cp_dsl[i]))
                            - cm_i * wing_i.cp_macs[i] * wing_i.u_s[i]
                        )
                        * wing_i.cp_areas[i]
                        / (S_ref * c_ref)
                    )

                # Rotate force coefficients to stability axis
                alpha_rad = result.alpha * np.pi / 180
                rotation_matrix = np.array(
                    [
                        [np.cos(alpha_rad), 0, np.sin(alpha_rad)],
                        [0, 1, 0],
                        [-1 * np.sin(alpha_rad), 0, np.cos(alpha_rad)],
                    ]
                )
                CF = np.matmul(rotation_matrix, CF)
                surface_coefficients[wing_i.surface_name] = Coefficients(CF, CM, cl_distribution)

            CF_global = np.zeros(3)
            CM_global = np.zeros(3)
            for _, coefficients in surface_coefficients.items():
                CF_global += np.array([coefficients.forces.CD, coefficients.forces.CY, coefficients.forces.CL])
                CM_global += np.array([coefficients.moments.Cl, coefficients.moments.Cm, coefficients.moments.Cn])
            global_coefficients = Coefficients(CF_global, CM_global, {})
            processed_result = ProcessedSimulationResults(result, global_coefficients, surface_coefficients)
            processed_simulation_results.append(processed_result)
        return processed_simulation_results
    
    @classmethod
    def get_aerodynamic_center(
        cls,
        original_wing_pool: WingPool,
        simulation_results: list[SimulationResult],
        S_ref: Optional[float] = None,
        c_ref: Optional[float] = None,
        max_iter: int = 100,
        residual: float = 1e-6
    ) -> dict:
        wing_pool = deepcopy(original_wing_pool) # Make a copy to avoid updating original object
        if np.any(wing_pool.flight_condition.angular_velocity != 0):
            print("Warning: trying to find aerodynamic center with non-zero angular rates may result in inconsistent values")
        
        if not S_ref:
            S_ref = wing_pool.S_ref
        if not c_ref:
            c_ref = wing_pool.c_ref

        ac_min = min([wing.x_pos for wing in wing_pool.pool])
        ac_max = max([wing.x_pos for wing in wing_pool.pool]+[c_ref])
        ac_start = (ac_min + ac_max) / 2

        ac = ac_start
        wing_pool.moment_ref = [ac, 0, 0]
        ac_check = False

        i = 1
        alpha_values = [result.alpha for result in simulation_results]
        while not ac_check:
            if i > max_iter:
                print("Reached max iteration number")
                break

            coefficients = cls.get_coefficients(wing_pool, simulation_results, S_ref, c_ref)

            CM_values = [coef.global_coefficients.moments.Cm for coef in coefficients]
            CM_poly = np.polyfit(alpha_values, CM_values, 1)

            Cm_alpha = CM_poly[0] # get angular coefficient, or A in y = Ax + b
            if abs(Cm_alpha) <= residual or i == max_iter:
                Cm_ac = CM_poly[1]
                ac_check = True
                x_ac = ac / c_ref
            else:
                if Cm_alpha < 0:
                    ac_min = ac
                    ac = (ac + ac_max) / 2
                else:
                    ac_max = ac
                    ac = (ac_min + ac_max) / 2
                wing_pool.moment_ref = [ac, 0, 0]
            i += 1

        return AerodynamicCenter(x_ac, Cm_ac, Cm_alpha, i)
