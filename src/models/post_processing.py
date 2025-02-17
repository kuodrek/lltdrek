from typing import List, Union
import numpy as np
from dataclasses import dataclass, field
from lltdrek.models.wingpool import WingPool
from lltdrek.models.wing import Wing
from lltdrek.models.flight_condition import FlightCondition
from lltdrek.utils.lookup import get_airfoil_data, get_linear_data_and_clmax


@dataclass(repr=False, eq=False, match_args=False, slots=True)
class PostProcessing:
    ref_point: list
    ref_points_dict: dict = field(init=False)


    def __post_init__(self) -> None:
        self.ref_points_dict = None


    def check_for_nan(self, G_dict) -> bool:
        convergence_check = True
        for surface, G_distribution in G_dict.items():
            if True in np.isnan(G_distribution): convergence_check = False

        return convergence_check


    def get_global_coefficients(self, wing_pool: WingPool, G_dict: dict[list], aoa_index: int, S_ref: float, c_ref: float) -> dict:
        if self.ref_points_dict is None:
            self.build_reference_points_dict(wing_pool)
        
        convergence_check = self.check_for_nan(G_dict)
        if not convergence_check:
            CF = np.array((np.nan, np.nan, np.nan))
            CM = np.array((np.nan, np.nan, np.nan))
            Cl_distr_dict = {}
            for wing in wing_pool.complete_wing_pool:
                Cl_distr_dict[wing.surface_name] = np.nan
            return {
                "CF": CF,
                "CM": CM,
                "Cl_distr": Cl_distr_dict
            }

        v_inf = wing_pool.flight_condition.v_inf_list[aoa_index]
        aoa = wing_pool.flight_condition.aoa[aoa_index]
        aoa_rad = aoa * np.pi / 180

        CF = np.zeros(3)
        CM = np.zeros(3)
        CF_distr = np.zeros((wing_pool.total_panels, 3))
        CM_distr = np.zeros((wing_pool.total_panels, 3))
        i_glob = 0
        Cl_distr_dict = {}
        for wing_i in wing_pool.complete_wing_pool:
            ref_points_distr = self.ref_points_dict[wing_i.surface_name]
            Cl_distr = np.zeros(wing_i.N_panels)
            G_i = G_dict[wing_i.surface_name]
            for i, _ in enumerate(wing_i.collocation_points):
                aux_cf = np.zeros(3)

                airfoil_coefficients = get_linear_data_and_clmax(
                    wing_i.cp_airfoils[i],
                    wing_i.cp_reynolds[i],
                    wing_i.airfoil_data
                )
                cm_i = airfoil_coefficients["cm0"]
                Cl_distr[i] = 2 * G_i[i]
                for wing_j in wing_pool.complete_wing_pool:
                    G_j = G_dict[wing_j.surface_name]
                    v_ij_distr = wing_pool.ind_velocities_list[aoa_index][wing_i.surface_name][wing_j.surface_name]
                    for j, _ in enumerate(wing_j.collocation_points):
                        v_ij = v_ij_distr[i][j]
                        aux_cf += G_j[j] * v_ij
                aux_cf = (aux_cf + v_inf) * G_i[i]
                CF += 2 * np.cross(aux_cf, wing_i.cp_dsl[i]) * wing_i.cp_areas[i] / S_ref
                CF_distr[i_glob,:] = 2 * np.cross(aux_cf, wing_i.cp_dsl[i]) * wing_i.cp_areas[i] / S_ref
                CM += (2 * np.cross(ref_points_distr[i], np.cross(aux_cf, wing_i.cp_dsl[i])) - cm_i * wing_i.cp_macs[i] * wing_i.u_s[i]) * wing_i.cp_areas[i] / ( S_ref * c_ref )
                CM_distr[i_glob,:] = (2 * np.cross(ref_points_distr[i], np.cross(aux_cf, wing_i.cp_dsl[i])) - cm_i * wing_i.cp_macs[i] * wing_i.u_s[i]) * wing_i.cp_areas[i] / ( S_ref * c_ref )

                i_glob += 1
            Cl_distr_dict[wing_i.surface_name] = Cl_distr

        # Rebater os coeficientes de forças para o eixo do escoamento
        rotation_matrix = np.array(
            [[np.cos(aoa_rad), 0, np.sin(aoa_rad)],
            [0, 1, 0],
            [-1*np.sin(aoa_rad), 0, np.cos(aoa_rad)]]
        )
        CF = np.matmul(rotation_matrix, CF)
        return {
            "CF": CF,
            "CM": CM,
            "Cl_distr": Cl_distr_dict
        }

    
    def build_reference_points_dict(self, wing_pool: WingPool) -> dict:
        self.ref_points_dict = {}
        for wing_i in wing_pool.wing_list:
            ref_point_distr = np.zeros((wing_i.N_panels, 3))
            ref_point_distr_mirrored  = np.zeros((wing_i.N_panels, 3))
            for i, cp_i in enumerate(wing_i.collocation_points):
                ref_point_distr[i,:] = cp_i - self.ref_point
                ref_point_distr_mirrored[i,:] = cp_i - self.ref_point
                ref_point_distr_mirrored[i,1] = ref_point_distr_mirrored[i,1] * -1

            self.ref_points_dict[wing_i.surface_name] = ref_point_distr
            self.ref_points_dict[wing_i.surface_name+"_mirrored"] = ref_point_distr_mirrored
        
        return self.ref_points_dict


    def get_wing_coefficients(self, wing_pool: WingPool, G_dict: dict[list], aoa_index: int, S_ref: float, c_ref: float) -> dict:
        if self.ref_points_dict is None:
            self.build_reference_points_dict(wing_pool)
        
        wing_coefficients = {}

        convergence_check = self.check_for_nan(G_dict)
        if not convergence_check:
            CF = np.array((np.nan, np.nan, np.nan))
            CM = np.array((np.nan, np.nan, np.nan))
            Cl_distr_dict = {}
            for wing in wing_pool.complete_wing_pool:
                Cl_distr_dict[wing.surface_name] = np.nan
                wing_coefficients[wing.surface_name] = {
                    "CF": CF,
                    "CM": CM,
                    "Cl_distr": Cl_distr_dict
                }
            return wing_coefficients
        
        v_inf = wing_pool.flight_condition.v_inf_list[aoa_index]
        aoa = wing_pool.flight_condition.aoa[aoa_index]
        aoa_rad = aoa * np.pi / 180

        CF_distr = np.zeros((wing_pool.total_panels, 3))
        CM_distr = np.zeros((wing_pool.total_panels, 3))
        i_glob = 0
        Cl_distr_dict = {}
        for wing_i in wing_pool.complete_wing_pool:
            CF = np.zeros(3)
            CM = np.zeros(3)
            
            ref_points_distr = self.ref_points_dict[wing_i.surface_name]
            Cl_distr = np.zeros(wing_i.N_panels)
            G_i = G_dict[wing_i.surface_name]
            for i, _ in enumerate(wing_i.collocation_points):
                aux_cf = np.zeros(3)

                airfoil_coefficients = get_linear_data_and_clmax(
                    wing_i.cp_airfoils[i],
                    wing_i.cp_reynolds[i],
                    wing_i.airfoil_data
                )
                cm_i = airfoil_coefficients["cm0"]
                Cl_distr[i] = 2 * G_i[i]
                for wing_j in wing_pool.complete_wing_pool:
                    G_j = G_dict[wing_j.surface_name]
                    v_ij_distr = wing_pool.ind_velocities_list[aoa_index][wing_i.surface_name][wing_j.surface_name]
                    for j, _ in enumerate(wing_j.collocation_points):
                        v_ij = v_ij_distr[i][j]
                        aux_cf += G_j[j] * v_ij
                aux_cf = (aux_cf + v_inf) * G_i[i]
                CF += 2 * np.cross(aux_cf, wing_i.cp_dsl[i]) * wing_i.cp_areas[i] / S_ref
                CF_distr[i_glob,:] = 2 * np.cross(aux_cf, wing_i.cp_dsl[i]) * wing_i.cp_areas[i] / S_ref
                CM += (2 * np.cross(ref_points_distr[i], np.cross(aux_cf, wing_i.cp_dsl[i])) - cm_i * wing_i.cp_macs[i] * wing_i.u_s[i]) * wing_i.cp_areas[i] / ( S_ref * c_ref )
                CM_distr[i_glob,:] = (2 * np.cross(ref_points_distr[i], np.cross(aux_cf, wing_i.cp_dsl[i])) - cm_i * wing_i.cp_macs[i] * wing_i.u_s[i]) * wing_i.cp_areas[i] / ( S_ref * c_ref )
                i_glob += 1

            rotation_matrix = np.array(
                [[np.cos(aoa_rad), 0, np.sin(aoa_rad)],
                [0, 1, 0],
                [-1*np.sin(aoa_rad), 0, np.cos(aoa_rad)]]
            )
            CF = np.matmul(rotation_matrix, CF)
            Cl_distr_dict[wing_i.surface_name] = Cl_distr
            wing_coefficients[wing_i.surface_name] = {
                "CF": CF,
                "CM": CM,
                "Cl_distr": Cl_distr_dict
            }

        for wing_i in wing_pool.wing_list:
            # Somar as parcelas das asas (wing + wing_mirrored) e remover chaves de asas espelhadas
            wing_coefficients[wing_i.surface_name]["CF"] += wing_coefficients[wing_i.surface_name+"_mirrored"]["CF"]
            wing_coefficients[wing_i.surface_name]["CM"] += wing_coefficients[wing_i.surface_name+"_mirrored"]["CM"]
            wing_coefficients.pop(wing_i.surface_name+"_mirrored")

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
        O input aoa_index_range pode ser utilizado para procurar numa faixa de angulos específica
        por exemplo, se aoa_list [1, 3, 5, 7, 9, 11, 13, 15, 19]
        passando-se aoa_range = [13, 19] será iterado somente entre os ângulos 13 e 19
        """
        c_ref = 1 # O valor de c_ref pode ser qualquer um aqui porque nao vamos pegar coeficiente de momento

        if aoa_range == None:
            aoa_start = wing_pool.flight_condition.aoa[0]
            aoa_end = wing_pool.flight_condition.aoa[-1]
        else:
            aoa_start = aoa_range[0]
            aoa_end = aoa_range[-1]

        CLmax_dict = {}
        for wing_i in wing_pool.wing_list:
            CLmax_check = False
            for aoa_idx, aoa in enumerate(wing_pool.flight_condition.aoa):
                if aoa_start <= aoa <= aoa_end and not CLmax_check:
                    G_dict = G_list[aoa_idx]
                    convergence_check = self.check_for_nan(G_dict)
                    if not convergence_check: continue
                    G_i = G_dict[wing_i.surface_name]
                    for i, _ in enumerate(wing_i.collocation_points):
                        lookup_data = get_linear_data_and_clmax(
                            wing_i.cp_airfoils[i],
                            wing_i.cp_reynolds[i],
                            wing_i.airfoil_data
                        )
                        Clmax_i = lookup_data["clmax"]
                        if 2 * G_i[i] > Clmax_i:
                            CLmax_check = True
                            aoa_max_idx = aoa_idx - 1 if aoa_idx > 0 else 0
            if not CLmax_check:
                aoa_max_idx = -1 if aoa_range == None else wing_pool.flight_condition.aoa.index(aoa_range[-1])
            G_dict = G_list[aoa_max_idx]
            aoa_max = wing_pool.flight_condition.aoa[aoa_max_idx]
            coefficients = self.get_wing_coefficients(wing_pool, G_dict, aoa_max_idx, S_ref, c_ref)
            CLmax_dict[wing_i.surface_name] = {
                "CLmax": coefficients[wing_i.surface_name]["CF"][2],
                "aoa_max": aoa_max,
                "aoa_max_idx": aoa_max_idx
            }

        return CLmax_dict


    def get_aerodynamic_center(self, wing_pool: WingPool, G_list: list[dict], aoa_1: float, aoa_2: float, S_ref: float, c_ref: float) -> dict:
        aoa_list = wing_pool.flight_condition.aoa
        if aoa_1 not in aoa_list or aoa_2 not in aoa_list:
            raise Exception(f"Ambos aoa_1 e aoa_2 precisam estar dentro da lista de ângulos de ataque")

        aoa_1_idx = aoa_list.index(aoa_1)
        aoa_2_idx = aoa_list.index(aoa_2)

        ac = 0.25*c_ref
        ac_min = 0.1*c_ref
        ac_max = 1.5*c_ref
        max_iter = 100
        
        self.ref_point = [ac, 0, 0]
        ac_check = False
        
        i = 1
        while not ac_check:
            if i > max_iter:
                print("Limite máximo de iterações atingido")
                break
            
            self.build_reference_points_dict(wing_pool)

            coefs_1 = self.get_global_coefficients(wing_pool, G_list[0], aoa_index=aoa_1_idx, S_ref=S_ref, c_ref=c_ref)
            Cm_1 = coefs_1["CM"][1]
            coefs_2 = self.get_global_coefficients(wing_pool, G_list[1], aoa_index=aoa_2_idx, S_ref=S_ref, c_ref=c_ref)
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

        return {
            "x_ac": x_ac,
            "Cm_ac": Cm_ac
        }
