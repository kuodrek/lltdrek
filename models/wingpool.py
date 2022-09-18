from typing import List
from models.wing import Wing
from models.flight_condition import FlightCondition
from utils import velocity


wingpool_dict = {
    "surface_name": {
        "surface_name": "v_ind"
    }
}

class WingPool:
    def __init__(self, wing_list: List[Wing], flight_condition: FlightCondition):
        WingPool.wing_list = wing_list
        WingPool.flight_condition = flight_condition
        WingPool.wingpool_velocities = None
    

    def calculate_induced_velocities(self):
        for wing_cp in self.wing_list:
            wingpool_velocities = {}
            for wing_vp in self.wing_list:
                v_ind_vetor = 1
                wingpool_velocities[wing_cp.surface_name][wing_vp.surface_name] = v_ind_vetor
        self.wingpool_velocities = wingpool_velocities
        