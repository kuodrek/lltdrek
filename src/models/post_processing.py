from dataclasses import dataclass
from typing import List, Optional, Union

import numpy as np

from src.models.simulation import SimulationResult
from src.models.wingpool import WingPool
from src.utils.lookup import get_linear_data_and_clmax


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

    def get_wing_coefficients(
        self, wing_pool: WingPool, G_dict: dict[list], alpha_index: int, S_ref: float, c_ref: float
    ) -> dict:
        wing_coefficients = {}

        v_inf = wing_pool.flight_condition.v_inf_list[alpha_index]
        aoa = wing_pool.flight_condition.angles_of_attack[alpha_index]
        aoa_rad = aoa * np.pi / 180

        CF_distr = np.zeros((wing_pool.total_panels, 3))
        CM_distr = np.zeros((wing_pool.total_panels, 3))
        i_glob = 0
        Cl_distr_dict = {}
        for wing_i in wing_pool.pool:
            CF = np.zeros(3)
            CM = np.zeros(3)

            ref_points_distr = wing_pool.system_moment_ref[wing_i.surface_name]
            Cl_distr = np.zeros(wing_i.N_panels)
            G_i = G_dict[wing_i.surface_name]
            for i, _ in enumerate(wing_i.collocation_points):
                aux_cf = np.zeros(3)

                airfoil_coefficients = get_linear_data_and_clmax(
                    wing_i.cp_airfoils[i], wing_i.cp_reynolds[i], wing_i.airfoil_data, show_logs=False
                )
                cm_i = airfoil_coefficients["cm0"]
                Cl_distr[i] = 2 * G_i[i]
                for wing_j in wing_pool.pool:
                    G_j = G_dict[wing_j.surface_name]
                    v_ij_distr = wing_pool.ind_velocities_list[alpha_index][wing_i.surface_name][wing_j.surface_name]
                    for j, _ in enumerate(wing_j.collocation_points):
                        v_ij = v_ij_distr[i][j]
                        aux_cf += G_j[j] * v_ij
                aux_cf = (aux_cf + v_inf) * G_i[i]
                CF += 2 * np.cross(aux_cf, wing_i.cp_dsl[i]) * wing_i.cp_areas[i] / S_ref
                CF_distr[i_glob, :] = 2 * np.cross(aux_cf, wing_i.cp_dsl[i]) * wing_i.cp_areas[i] / S_ref
                CM += (
                    (
                        2 * np.cross(ref_points_distr[i], np.cross(aux_cf, wing_i.cp_dsl[i]))
                        - cm_i * wing_i.cp_macs[i] * wing_i.u_s[i]
                    )
                    * wing_i.cp_areas[i]
                    / (S_ref * c_ref)
                )
                CM_distr[i_glob, :] = (
                    (
                        2 * np.cross(ref_points_distr[i], np.cross(aux_cf, wing_i.cp_dsl[i]))
                        - cm_i * wing_i.cp_macs[i] * wing_i.u_s[i]
                    )
                    * wing_i.cp_areas[i]
                    / (S_ref * c_ref)
                )
                i_glob += 1

            rotation_matrix = np.array(
                [[np.cos(aoa_rad), 0, np.sin(aoa_rad)], [0, 1, 0], [-1 * np.sin(aoa_rad), 0, np.cos(aoa_rad)]]
            )
            CF = np.matmul(rotation_matrix, CF)
            Cl_distr_dict[wing_i.surface_name] = Cl_distr
            wing_coefficients[wing_i.surface_name] = {"CF": CF, "CM": CM, "Cl_distr": Cl_distr_dict}

        for wing_i in wing_pool.wing_list:
            # Somar as parcelas das asas (wing + wing_mirrored) e remover chaves de asas espelhadas
            wing_coefficients[wing_i.surface_name]["CF"] += wing_coefficients[wing_i.surface_name + "_mirrored"]["CF"]
            wing_coefficients[wing_i.surface_name]["CM"] += wing_coefficients[wing_i.surface_name + "_mirrored"]["CM"]
            wing_coefficients.pop(wing_i.surface_name + "_mirrored")

        return wing_coefficients

    def get_CL_max_linear(
        self,
        wing_pool: WingPool,
        G_list: list[dict],
        S_ref: float,
        aoa_range: Union[tuple[float], None] = None,
    ) -> float:
        """
        Obtenção do CL máximo de uma asa (ou sistema de asas) através do método da seção crítica
        Para cada ângulo é verificado o Cl de cada seção e comparado com o Clmax do perfil no seu respectivo reynolds
        O CLmax é alcançado quando alguma seção atingir o Clmax do seu perfil no seu respectivo reynolds
        O input alpha_index_range pode ser utilizado para procurar numa faixa de angulos específica
        por exemplo, se aoa_list [1, 3, 5, 7, 9, 11, 13, 15, 19]
        passando-se aoa_range = [13, 19] será iterado somente entre os ângulos 13 e 19
        """
        c_ref = 1  # O valor de c_ref pode ser qualquer um aqui porque nao vamos pegar coeficiente de momento

        if aoa_range is None:
            aoa_start = wing_pool.flight_condition.angles_of_attack[0]
            aoa_end = wing_pool.flight_condition.angles_of_attack[-1]
        else:
            aoa_start = aoa_range[0]
            aoa_end = aoa_range[-1]

        CLmax_dict = {}
        for wing_i in wing_pool.wing_list:
            CLmax_check = False
            for aoa_idx, aoa in enumerate(wing_pool.flight_condition.angles_of_attack):
                if aoa_start <= aoa <= aoa_end and not CLmax_check:
                    G_dict = G_list[aoa_idx]
                    convergence_check = self.check_for_nan(G_dict)
                    if not convergence_check:
                        continue
                    G_i = G_dict[wing_i.surface_name]
                    for i, _ in enumerate(wing_i.collocation_points):
                        lookup_data = get_linear_data_and_clmax(
                            wing_i.cp_airfoils[i], wing_i.cp_reynolds[i], wing_i.airfoil_data, show_logs=False
                        )
                        Clmax_i = lookup_data["clmax"]
                        if 2 * G_i[i] > Clmax_i:
                            CLmax_check = True
                            aoa_max_idx = aoa_idx - 1 if aoa_idx > 0 else 0
            if not CLmax_check:
                aoa_max_idx = (
                    -1 if aoa_range is None else wing_pool.flight_condition.angles_of_attack.index(aoa_range[-1])
                )
            G_dict = G_list[aoa_max_idx]
            aoa_max = wing_pool.flight_condition.angles_of_attack[aoa_max_idx]
            coefficients = self.get_wing_coefficients(wing_pool, G_dict, aoa_max_idx, S_ref, c_ref)
            CLmax_dict[wing_i.surface_name] = {
                "CLmax": coefficients[wing_i.surface_name]["CF"][2],
                "aoa_max": aoa_max,
                "aoa_max_idx": aoa_max_idx,
            }

        return CLmax_dict

    def get_aerodynamic_center(
        self, wing_pool: WingPool, G_list: list[dict], aoa_1: float, aoa_2: float, S_ref: float, c_ref: float
    ) -> dict:
        aoa_list = wing_pool.flight_condition.angles_of_attack
        if aoa_1 not in aoa_list or aoa_2 not in aoa_list:
            raise Exception("aoa_1 and aoa_2 must be in aoa_list")

        aoa_1_idx = aoa_list.index(aoa_1)
        aoa_2_idx = aoa_list.index(aoa_2)

        ac = 0.25 * c_ref
        ac_min = 0.1 * c_ref
        ac_max = 1.5 * c_ref
        max_iter = 100

        self.ref_point = [ac, 0, 0]
        ac_check = False

        i = 1
        while not ac_check:
            if i > max_iter:
                print("Reached max iteration number")
                break

            self.build_reference_points_dict(wing_pool)

            coefs_1 = self.get_global_coefficients(
                wing_pool, G_list[0], alpha_index=aoa_1_idx, S_ref=S_ref, c_ref=c_ref
            )
            coefs_1 = self.get_coefficients(
                wing_pool,
            )
            Cm_1 = coefs_1["CM"][1]
            coefs_2 = self.get_global_coefficients(
                wing_pool, G_list[1], alpha_index=aoa_2_idx, S_ref=S_ref, c_ref=c_ref
            )
            Cm_2 = coefs_2["CM"][1]

            Cm_alpha = (Cm_2 - Cm_1) / (aoa_2 - aoa_1)
            if abs(Cm_alpha) <= 1e-6 or i == max_iter:
                Cm_ac = Cm_1
                ac_check = True
                x_ac = ac / c_ref
            else:
                if Cm_alpha < 0:
                    ac_min = ac
                    ac = (ac + ac_max) / 2
                    self.ref_point = [ac, 0, 0]
                else:
                    ac_max = ac
                    ac = (ac_min + ac_max) / 2
                    self.ref_point = [ac, 0, 0]
            i += 1

        return {"x_ac": x_ac, "Cm_ac": Cm_ac}
