from typing import List
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


    def get_global_coefficients(self, wing_pool: WingPool, G_list: list[dict]) -> dict:
        G_dict = {}
        v_inf_array = np.zeros(3)
        aoa_idx = 0
        Sref = 1
        CF = np.zeros(3)
        for wing_i in wing_pool.complete_wing_pool:
            aux_cf = np.zeros(3)
            G_i = G_dict[wing_i.surface_name]
            for i, _ in enumerate(wing_i.collocation_points):
                aux_cf += G_i[i] * v_inf_array
                for wing_j in wing_pool.complete_wing_pool:
                    G_j = G_dict[wing_j.surface_name]
                    v_ij_distr = wing_pool.ind_velocities_list[aoa_idx][wing_i.surface_name][wing_j.surface_name]
                    for j, _ in enumerate(wing_j.collocation_points):
                        v_ij = v_ij_distr[i][j]
                        aux_cf += G_i[i] * G_j[j] * v_ij
            CF += np.cross(aux_cf, wing_i.cp_dsl[i]) * wing_i.cp_areas[i] / Sref
            # implementar cm, que também usa np.cross(aux_cf, wing_i.cp_dsl[i])
        CF = CF * 2


    def get_wing_coefficients(self):
        pass


    def get_CL_max_linear(self, wing_pool: WingPool, G_list: list) -> float:
        pass


    def get_aerodynamic_center(self, wing_pool: WingPool, G_list: list) -> float:
        pass
