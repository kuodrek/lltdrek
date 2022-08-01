class FlightConditions:
    """
    Velocidade do escoamento - Vinf
    Coeficiente de viscosidade dinamica do ar - nu
    densidade do ar - rho
    altura do chao - h
    angulo de ataque - aoa
    """
    def __init__(self, V_inf, nu, rho, h=0):
        FlightConditions.V_inf = V_inf
        FlightConditions.nu = nu
        FlightConditions.rho = rho
        FlightConditions.h = h
