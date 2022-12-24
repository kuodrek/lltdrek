from typing import List, Union
import numpy as np
from dataclasses import dataclass, field
from models.wingpool import WingPool
from models.wing import Wing
from models.flight_condition import FlightCondition

@dataclass(repr=False, eq=False, match_args=False)
class PostProcessing:
    reference_point: list


    def __post_init__(self) -> None:
        pass


    def get_global_coefficients(self, wing_pool: WingPool, G_dict: dict[list], aoa_index: int) -> dict:
        v_inf_array = np.zeros(3)
        aoa_idx = 0
        Sref = 1
        cref = 1
        CF = np.zeros(3)
        CM = np.zeros(3)
        for wing_i in wing_pool.complete_wing_pool:
            aux_cf = np.zeros(3)
            G_i = G_dict[wing_i.surface_name]
            for i, _ in enumerate(wing_i.collocation_points):
                cm_i = 1
                aux_cf += G_i[i] * v_inf_array
                for wing_j in wing_pool.complete_wing_pool:
                    G_j = G_dict[wing_j.surface_name]
                    v_ij_distr = wing_pool.ind_velocities_list[aoa_idx][wing_i.surface_name][wing_j.surface_name]
                    for j, _ in enumerate(wing_j.collocation_points):
                        v_ij = v_ij_distr[i][j]
                        aux_cf += G_i[i] * G_j[j] * v_ij
                CF += 2 * np.cross(aux_cf, wing_i.cp_dsl[i]) * wing_i.cp_areas[i] / Sref
                CM += (2 * np.cross(r_list[i], np.cross(aux_cf, wing_i.cp_dsl[i])) - cm_i * wing_i.chords[i] * wing_i.u_s[i]) * wing_i.cp_areas[i] / ( Sref * cref )
        # CF = CF * 2

        return {
            "CF": CF,
            "CM": CM
        }


    def get_wing_coefficients(self):
        pass


    def get_CL_max_linear(self, wing_pool: WingPool, G_list: list[dict], aoa_range: Union[tuple[float], None]) -> float:
        """
        Obtenção do CL máximo de uma asa (ou sistema de asas) através do método da seção crítica
        Para cada ângulo é verificado o Cl de cada seção e comparado com o Clmax do perfil no seu respectivo reynolds
        O CLmax é alcançado quando alguma seção atingir o Clmax do seu perfil no seu respectivo reynolds
        O input aoa_index_range pode ser utilizado para procurar numa faixa de angulos específica
        por exemplo, se aoa_list [1, 3, 5, 7, 9, 11, 13, 15, 19]
        passando-se aoa_range = [13, 19] será iterado somente entre os ângulos 13 e 19
        """
        # Iterar sobre a lista de G's (da faixa de angulos de wing_pool.flight_condition)
        #   Obter a distribuição de Cl's da(s) asa(s)
        #   Verificar se max(Cl_distr) > Clmax(perfil,  Re)
        #   Caso sim, o alfa_max = aoa anterior 
        #   Caso contrário continuar iterando
        pass


    def cl_distribution(self, wing_pool: WingPool, G_dict: dict[list], aoa_index: int) -> list[float]:
        pass


    def get_aerodynamic_center(self, wing_pool: WingPool, G_list: list[dict]) -> float:
        pass
