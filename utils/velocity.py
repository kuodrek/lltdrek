import numpy as np
import numpy.linalg as npla
from models.wing import Wing


# Distribution of induced velocities
def get_induced_velocity_distribution(collocation_points, cp_macs, vertice_points, v_inf_array, surface_name):
    v_ij_distr = np.zeros([len(collocation_points),len(vertice_points)-1,3])
    len_cp = len(collocation_points)
    len_vp = len(vertice_points)
    for i in range(len_cp):
        cp_i = collocation_points[i]
        mac_i = cp_macs[i]
        for j in range(len_vp - 1):
            if "_mirrored" in surface_name:
                # Pelo fato da asa ser espelhada, o sentido do vetor de posição também é invertido
                vp_jj = vertice_points[j]
                vp_j = vertice_points[j+1]
            else:
                vp_j = vertice_points[j]
                vp_jj = vertice_points[j+1]
            v_ij = get_induced_velocity(cp_i, vp_j, vp_jj, mac_i, v_inf_array)
            v_ij_distr[i,j,:] = v_ij

    return v_ij_distr


# Function to calculate the induced velocity caused by a vortex panel in a point
def get_induced_velocity(collocation_point, vertice_point_1, vertice_point_2, mac, v_inf_array): 
    velocity_ij = np.zeros(3)

    ri1j = collocation_point - vertice_point_1
    ri2j = collocation_point - vertice_point_2

    ri1j_abs = npla.norm(ri1j)
    ri2j_abs = npla.norm(ri2j)
    r1_cross_prod = np.cross(v_inf_array, ri1j)
    r2_cross_prod = np.cross(v_inf_array, ri2j)
    r1_dot_prod = np.dot(v_inf_array, ri1j)
    r2_dot_prod = np.dot(v_inf_array, ri2j)
    r12_cross_prod = np.cross(ri1j, ri2j)
    r12_dot_prod = np.dot(ri1j, ri2j)

    velocity_ij = mac / (4 * np.pi) * \
        ( r2_cross_prod / (ri2j_abs*(ri2j_abs-r2_dot_prod)) \
        - r1_cross_prod / (ri1j_abs*(ri1j_abs-r1_dot_prod)) )
    
    bound_vortex_den = (ri1j_abs*ri2j_abs*(ri1j_abs*ri2j_abs+r12_dot_prod))

    if bound_vortex_den != 0:
        velocity_ij +=  mac / (4 * np.pi) * \
            ( (ri1j_abs+ri2j_abs) * r12_cross_prod / bound_vortex_den )
    
    return velocity_ij


# Induced velocity calculation considering ground effect through reflection method
def get_induced_velocity_ground_effect(collocation_point, vertice_point_1, vertice_point_2, v_inf_array, mac, h):
    velocity_ij = np.zeros(3)

    v_inf_ge = v_inf_array
    v_inf_ge[2] = -v_inf_ge[2]

    ri1j = collocation_point - vertice_point_1
    ri1j[2] = collocation_point[2] - (2*h + 2*vertice_point_1[2])
    ri2j = collocation_point - vertice_point_2
    ri2j[2] = collocation_point[2] - (2*h + 2*vertice_point_2[2])

    ri1j_abs = npla.norm(ri1j)
    ri2j_abs = npla.norm(ri2j)
    r1_cross_prod = np.cross(v_inf_ge, ri1j)
    r2_cross_prod = np.cross(v_inf_ge, ri2j)
    r1_dot_prod = np.dot(v_inf_ge, ri1j)
    r2_dot_prod = np.dot(v_inf_ge, ri2j)
    r12_cross_prod = np.cross(ri1j, ri2j)
    r12_dot_prod = np.dot(ri1j, ri2j)

    bound_vortex_den = (ri1j_abs*ri2j_abs*(ri1j_abs*ri2j_abs+r12_dot_prod))

    velocity_ij = mac / (4 * np.pi) * \
        ( r2_cross_prod / (ri2j_abs*(ri2j_abs-r2_dot_prod)) \
        - r1_cross_prod / (ri1j_abs*(ri1j_abs-r1_dot_prod)) )
    
    if bound_vortex_den != 0:
        velocity_ij +=  mac / (4 * np.pi) * \
            ( (ri1j_abs+ri2j_abs) * r12_cross_prod / bound_vortex_den )
    
    return velocity_ij
