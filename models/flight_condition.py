from typing import List, Union
import numpy as np
from dataclasses import dataclass, field


@dataclass
class FlightCondition:
    V_inf: float
    nu: float
    rho: float
    aoa: Union[int, List[float]]
    h: float = 0
    ground_effect_check: bool = False
    v_inf_array: np.ndarray = field(init=False)

    def __post_init__(self):
        if type(self.aoa) == int:
            self.v_inf_array = np.array([np.cos(self.aoa), 0, np.sin(self.aoa)])
        else:
            self.v_inf_array = np.zeros([len(self.aoa), 3])
            for i, aoa_i in enumerate(self.aoa):
                aoa_i = aoa_i * np.pi / 180
                self.v_inf_array[i,:] = [np.cos(aoa_i), 0, np.sin(aoa_i)]
