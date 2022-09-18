from typing import List
from models.wing import Wing
from models.flight_condition import FlightCondition
from utils import velocity


wingpool_dict = {
    "flight_condition": {
        "ground_effect_check": False,
        "height": "height_value"
    },
    "surfaces": {
    "surface_name": {
        "id": 1,
        "wing": "wing_object",
    }
    }
}

class WingPool:
    def __init__(self, wing_list: List[Wing], flight_condition: FlightCondition):
        WingPool.wing_list = wing_list
        WingPool.flight_condition = flight_condition
        WingPool.wingpool_dict = None
    
    def generate_wingpool_dict(self):
        # Function that returns wingpool_dict
        wingpool_dict = {}
        wingpool_dict["flight_condition"]["ground_effect_check"] = self.flight_condition.ground_effect_check
        wingpool_dict["flight_condition"]["ground_effect_check"] = self.flight_condition.h
        for wing in self.wing_list:
            wingpool_dict["surfaces"][wing.surface_name] = 
            pass
        pass

    def calculate_induced_velocities(self):
        pass