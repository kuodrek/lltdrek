from typing import List
import numpy as np


class FlightCondition:
    """
    Velocidade do escoamento - Vinf
    Coeficiente de viscosidade dinamica do ar - nu
    densidade do ar - rho
    altura do chao - h
    angulo de ataque - aoa
    """
    def __init__(self, V_inf, nu, rho, aoa: List, h=0, ground_effect_check=False):
        FlightCondition.V_inf = V_inf
        FlightCondition.nu = nu
        FlightCondition.rho = rho
        FlightCondition.aoa = aoa
        FlightCondition.h = h
        FlightCondition.ground_effect_check = ground_effect_check
        FlightCondition.v_inf_array = np.zeros([len(aoa), 3])
        for i, aoa_i in enumerate(aoa):
            aoa_i = aoa_i * np.pi / 180
            self.v_inf_array[i,:] = [np.cos(aoa_i), 0, np.sin(aoa_i)]
