from typing import List
import numpy as np
from dataclasses import dataclass, field
from models.wingpool import WingPool
from models.wing import Wing
from models.flight_condition import FlightCondition

@dataclass
class PostProcessing:
    reference_point: list


    def __post_init__(self) -> None:
        pass


    def get_coefficients(self, wing_pool: WingPool, G_solution: list) -> dict:
        pass


    def get_CL_max_linear(self, wing_pool: WingPool, G_solution: list) -> float:
        pass


    def get_aerodynamic_center(self, wing_pool: WingPool, G_solution: list) -> float:
        pass
