from typing import Dict, List
from models.wing import Wing
from models.flight_condition import FlightCondition
from utils import velocity
from dataclasses import dataclass, field

wingpool_dict_template = {
    "surface_name": {
        "surface_name": "v_ind_distr"
    }
}

G_dict = {
    "surface_name": "G_array",
}

@dataclass(frozen=True)
class WingPool:
    wing_list: List[Wing]
    flight_condition: FlightCondition
    velocities_dict: Dict = field(init=False)
    G_dict: Dict = field(init=False)

    def calculate_induced_velocities(self):
        v_inf = self.flight_condition.v_inf_array
        for wing_cp in self.wing_list:
            velocities_dict = {}
            collocation_points = wing_cp.collocation_points
            cp_macs = wing_cp.cp_macs
            for wing_vp in self.wing_list:
                vertice_points = wing_vp.vertice_points
                v_ij_distr = velocity.get_induced_velocity_distribution(
                    collocation_points, cp_macs, vertice_points, v_inf
                )
                velocities_dict[wing_cp.surface_name][wing_vp.surface_name] = v_ij_distr
        self.velocities_dict = velocities_dict
    

    def update_solution(self, G, surface_name: str):
        self.G_dict[surface_name] = G
