import numpy as np
from models.wing import Wing
from models.simulation import Simulation
from models.flight_conditions import FlightConditions
from utils import velocity

def main_equation(Wing: Wing, Simulation: Simulation, FlightConditions: FlightConditions):
    N_panels = Wing.N_panels
    collocation_points = Wing.collocation_points
    vertice_points = Wing.vertice_points
    v_inf = 1
    mac = 1
    for i in range(N_panels):
        # Calcular velocidades induzidas
        cp_i = collocation_points[i]
        for j in range(N_panels):
            vp_i = vertice_points[j]
            vp_ii = vertice_points[j+1]
            v_ij = velocity.get_induced_velocity(cp_i, vp_i, vp_ii, v_inf, mac)
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
    