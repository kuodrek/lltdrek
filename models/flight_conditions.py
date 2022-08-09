import numpy as np


class FlightConditions:
    """
    Velocidade do escoamento - Vinf
    Coeficiente de viscosidade dinamica do ar - nu
    densidade do ar - rho
    altura do chao - h
    angulo de ataque - aoa
    """
    def __init__(self, V_inf, nu, rho, aoa, h=0):
        FlightConditions.V_inf = V_inf
        FlightConditions.nu = nu
        FlightConditions.rho = rho
        FlightConditions.aoa = aoa
        FlightConditions.h = h
        FlightConditions.v_inf_array = np.zeros([len(aoa), 3])
        for i, aoa_i in enumerate(aoa):
            aoa_i = aoa_i * np.pi / 180
            self.v_inf_array[i,:] = [np.cos(aoa_i), 0, np.sin(aoa_i)]
