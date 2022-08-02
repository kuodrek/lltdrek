import numpy as np

def get_div(collocation_point, vertice_point_1, vertice_point_2, v_inf, mac, h, same_panel_check): # Function to calculate induced velocities between 
    ri1j = np.zeros(3)
    ri2j = np.zeros(3)
    velocity_ij = np.zeros(3)

    for k in range(3):
        ri2j[k] = collocation_point[k] - vertice_point_2[k]
        ri1j[k] = collocation_point[k] - vertice_point_1[k]

    r1_cross_product = np.cross(v_inf, ri1j)
    r2_cross_product = np.cross(v_inf, ri2j)
    r1_dot_product = np.dot(v_inf, ri1j)
    r2_dot_product = np.dot(v_inf, ri2j)

    for k in range(3):
        velocity_ij[k] = mac / (4 * np.pi) * \
        ( r2_cross_product[k] / (np.absolute(ri2j[k])*(np.absolute(ri2j[k])-r2_dot_product)) \
        - r1_cross_product[k] / (np.absolute(ri1j[k])*(np.absolute(ri1j[k])-r1_dot_product)) )