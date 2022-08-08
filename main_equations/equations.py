import numpy as np
from models.wing import Wing
from models.simulation import Simulation
from models.flight_conditions import FlightConditions

def main_equation(Wing: Wing, Simulation: Simulation, FlightConditions: FlightConditions):
    N = Wing.N_panels
    for i in range(N):
        # Calcular velocidades induzidas
        vij = 1
        # Calcular angulos de ataque de cada seçao
        aoa_i = 1
        # Calcular os Cl's de cada seção
        cl_i = 1
        # Calculo do valor residual da funçao
        residual = 0

    return residual

def corrector_equation(Wing: Wing, Simulation: Simulation, FlightConditions: FlightConditions):
    delta_G = 0
    return delta_G
    