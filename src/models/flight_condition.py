from typing import Union, Sequence
import numpy as np

class FlightCondition:

    def __init__(
        self,
        V_inf: Union[float, int],
        nu: Union[float, int],
        rho: Union[float, int],
        angles_of_attack: Sequence[Union[float, int]],
        h: Union[float, int],
        ground_effect_check: bool = False,
    ):
        self.V_inf = V_inf
        self.nu = nu
        self.rho = rho
        self.angles_of_attack = angles_of_attack
        self.h = h
        self.ground_effect_check = ground_effect_check
        self.v_inf_list = [
            (
                np.cos(alpha * np.pi / 180),
                0,
                np.sin(alpha * np.pi / 180),
            )
            for alpha in self.angles_of_attack
        ]
