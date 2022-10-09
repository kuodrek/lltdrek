import numpy as np
import numpy.linalg as npla
from models.wing import Wing
from models.simulation import Simulation
from models.flight_condition import FlightCondition
from utils import velocity

def main_equation(Wing: Wing, FlightCondition: FlightCondition, v_ij_distr, G):
    N_panels = Wing.N_panels
    collocation_points = Wing.collocation_points
    vertice_points = Wing.vertice_points
    u_n = Wing.u_n
    u_a = Wing.u_a
    cp_dsl = Wing.cp_dsl
    cp_macs = Wing.cp_macs # Verificar se usa o mac de cada painel ou o MAC global
    v_inf = FlightCondition.v_inf_array

    residual = np.zeros(N_panels)
    velocity_sum = np.zeros(3)
    for i in range(N_panels):
        v_ij_panel = v_ij_distr[i,:,:]
        velocity_sum = v_ij_panel * G
        
        aoa_i = np.arctan(np.dot(velocity_sum+v_inf, u_n[i]) / np.dot(velocity_sum+v_inf, u_a[i]))
        # Calcular o Cl da seção
        cl_i = 1
        residual[i] = 2 * npla.norm( np.cross(velocity_sum+v_inf, cp_dsl[i]) )*G[i] - cl_i

    return residual


def corrector_equation(Wing: Wing, FlightCondition: FlightCondition, G):
    delta_G = None
    N_panels = Wing.N_panels
    collocation_points = Wing.collocation_points
    vertice_points = Wing.vertice_points
    u_n = Wing.u_n
    u_a = Wing.u_a
    cp_dsl = Wing.cp_dsl
    cp_macs = Wing.cp_macs # Verificar se usa o mac de cada painel ou o MAC global
    v_inf = FlightCondition.v_inf_array

    # v_local = v_inf + v_ind_i * G

    w_i = 1
    v_ni = 1
    v_ai = 1
    cl_alfa_i = 1

    return delta_G
    