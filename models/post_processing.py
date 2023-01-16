from typing import List, Union
import numpy as np
from dataclasses import dataclass, field
from models.wingpool import WingPool
from models.wing import Wing
from models.flight_condition import FlightCondition
from utils.lookup import get_airfoil_data, get_linear_data


@dataclass(repr=False, eq=False, match_args=False)
class PostProcessing:
    ref_point: list
    ref_points_dict: dict = field(init=False)


    def __post_init__(self) -> None:
        self.ref_points_dict = None

    def get_global_coefficients(self, wing_pool: WingPool, G_dict: dict[list], aoa_index: int, S_ref: float, c_ref: float) -> dict:
        if self.ref_points_dict is None:
            self.build_reference_points_dict(wing_pool)

        CF = np.zeros(3)
        CM = np.zeros(3)
        Cl_distr_dict = {}
        for wing_i in wing_pool.complete_wing_pool:
            ref_points_distr = self.ref_points_dict[wing_i.surface_name]
            Cl_distr = np.zeros(wing_i.N_panels)
            G_i = G_dict[wing_i.surface_name]
            for i, _ in enumerate(wing_i.collocation_points):
                aux_cf = np.zeros(3)

                airfoil_coefficients = get_linear_data(
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
                aux_cf = (aux_cf + wing_pool.flight_condition.v_inf_list[aoa_index]) * G_i[i]
                CF += 2 * np.cross(aux_cf, wing_i.cp_dsl[i]) * wing_i.cp_areas[i] / S_ref
                CM += (2 * np.cross(ref_points_distr[i], np.cross(aux_cf, wing_i.cp_dsl[i])) - cm_i * wing_i.cp_chords[i] * wing_i.u_s[i]) * wing_i.cp_areas[i] / ( S_ref * c_ref )
                a=1
            Cl_distr_dict[wing_i.surface_name] = Cl_distr

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


    def get_wing_coefficients(self):
        pass


    def get_CL_max_linear(
        self,
        wing_pool: WingPool,
        G_list: list[dict],
        aoa_range: Union[tuple[float],
        None],
        S_ref: float,
        c_ref: float
    ) -> float:
        """
        Obtenção do CL máximo de uma asa (ou sistema de asas) através do método da seção crítica
        Para cada ângulo é verificado o Cl de cada seção e comparado com o Clmax do perfil no seu respectivo reynolds
        O CLmax é alcançado quando alguma seção atingir o Clmax do seu perfil no seu respectivo reynolds
        O input aoa_index_range pode ser utilizado para procurar numa faixa de angulos específica
        por exemplo, se aoa_list [1, 3, 5, 7, 9, 11, 13, 15, 19]
        passando-se aoa_range = [13, 19] será iterado somente entre os ângulos 13 e 19
        """
        if aoa_range == None:
            aoa_start = wing_pool.flight_condition.aoa[0]
            aoa_end = wing_pool.flight_condition.aoa[-1]
        else:
            aoa_start = aoa_range[0]
            aoa_end = aoa_range[-1]

        CLmax_check = False
        for aoa_idx, aoa in enumerate(wing_pool.flight_condition.aoa):
            if aoa_start <= aoa <= aoa_end and CLmax_check:
                G_dict = G_list[aoa_idx]
                for wing_i in wing_pool.wing_list:
                    G_i = G_dict[wing_i.surface_name]
                    for i, _ in enumerate(wing_i.collocation_points):
                        Cl = get_airfoil_data( # tá errado, tem que ser o CLmax da seção
                            wing_i.cp_airfoils[i],
                            wing_i.cp_reynolds[i],
                            aoa,
                            wing_i.airfoil_data,
                            cl_alpha_check=False
                        )
                        if 2 * G_i[i] > Cl:
                            CLmax_check = True
                            aoa_max_idx = aoa_idx - 1 if aoa_idx > 0 else 0

        G_dict = G_list[aoa_max_idx]
        coefficients = self.get_global_coefficients(wing_pool, G_dict, aoa_max_idx, S_ref, c_ref)
        return coefficients["CF"][2]


    def get_aerodynamic_center(self, wing_pool: WingPool, G_list: list[dict]) -> float:
        pass
