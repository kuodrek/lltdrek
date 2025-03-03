import numpy as np
import numpy.linalg as npla


def get_induced_velocity_distribution(
    collocation_points: np.ndarray,
    cp_macs: np.ndarray,
    vertice_points: np.ndarray,
    v_inf_array: np.ndarray,
    surface_name: str,
    ground_effect_check: bool,
    h: float
    ) -> np.ndarray:
    """
    Distribution of induced velocities
    """
    len_cp = len(collocation_points)
    len_vp = len(vertice_points)
    v_ij_distr = np.zeros([len_cp,len_vp-1,3])

    for i in range(len_cp):
        cp_i = collocation_points[i]
        mac_i = cp_macs[i]
        for j in range(len_vp - 1):
            if "_mirrored" in surface_name:
                # Pelo fato da asa ser espelhada, o sentido do vetor posição também é invertido
                vp_jj = vertice_points[j]
                vp_j = vertice_points[j+1]
            else:
                vp_j = vertice_points[j]
                vp_jj = vertice_points[j+1]
            # print(f'{i}, {j}')
            v_ij, same_panel_check = get_induced_velocity(cp_i, vp_j, vp_jj, mac_i, v_inf_array)
            if ground_effect_check:
                # O método do plano reflexivo inverte a asa, então vp_j e vp_jj precisam ser invertidos para que os sentidos dos vórtices estejam corretos
                v_ij_ge = get_induced_velocity_ground_effect(cp_i, vp_jj, vp_j, mac_i, v_inf_array, h, same_panel_check)
                v_ij += v_ij_ge
            v_ij_distr[i,j,:] = v_ij

    return v_ij_distr


def get_induced_velocity(
    collocation_point: np.ndarray,
    vertice_point_1: np.ndarray,
    vertice_point_2: np.ndarray,
    mac: np.ndarray,
    v_inf_array: np.ndarray
    ) -> np.ndarray: 
    """
    Function to calculate the induced velocity caused by a vortex panel in a point
    """
    velocity_ij = np.zeros(3)

    ri1j = collocation_point - vertice_point_1
    ri2j = collocation_point - vertice_point_2

    ri1j_abs = npla.norm(ri1j)
    # print(f'ri1j_abs: {ri1j_abs}')
    ri2j_abs = npla.norm(ri2j)
    # print(f'ri2j_abs: {ri2j_abs}')
    r1_cross_prod = np.cross(v_inf_array, ri1j)
    # print(f'r1_cross_prod: {r1_cross_prod}')
    r2_cross_prod = np.cross(v_inf_array, ri2j)
    # print(f'r2_cross_prod: {r2_cross_prod}')
    r1_dot_prod = np.dot(v_inf_array, ri1j)
    # print(f'r1_dot_prod: {r1_dot_prod}')
    r2_dot_prod = np.dot(v_inf_array, ri2j)
    # print(f'r2_dot_prod: {r2_dot_prod}')
    r12_cross_prod = np.cross(ri1j, ri2j)
    # print(f'r12_cross_prod: {r12_cross_prod}')
    r12_dot_prod = np.dot(ri1j, ri2j)
    # print(f'r12_dot_prod: {r12_dot_prod}')

    velocity_ij = mac / (4 * np.pi) * \
        ( r2_cross_prod / (ri2j_abs*(ri2j_abs-r2_dot_prod)) \
        - r1_cross_prod / (ri1j_abs*(ri1j_abs-r1_dot_prod)) )
    # print(mac / (4 * np.pi))
    # print(r2_cross_prod / (ri2j_abs*(ri2j_abs-r2_dot_prod)))
    # print(r1_cross_prod / (ri1j_abs * (ri1j_abs - r1_dot_prod)))

    bound_vortex_den = (ri1j_abs*ri2j_abs*(ri1j_abs*ri2j_abs+r12_dot_prod))

    if not np.isclose(
        bound_vortex_den, 0, atol=1e-10
    ):  # Only add this if the panel is not inducing a velocity on itself
        bound_vortex_velocity =  mac / (4 * np.pi) * \
            ( (ri1j_abs+ri2j_abs) * r12_cross_prod / bound_vortex_den )
        velocity_ij += bound_vortex_velocity
        same_panel_check = False
    else:
        same_panel_check = True
    # print(velocity_ij)
    return velocity_ij, same_panel_check


def get_induced_velocity_ground_effect(collocation_point: np.ndarray,
    vertice_point_1: np.ndarray,
    vertice_point_2: np.ndarray,
    mac: float,
    v_inf_array: np.ndarray,
    h: float,
    same_panel_check: bool 
    ) -> np.ndarray:
    """
    Induced velocity calculation considering ground effect through reflection method
    """
    velocity_ij = np.zeros(3)

    v_inf_ge = np.array(v_inf_array)
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
    
    if not same_panel_check: # só adiciona essa parcela se a velocidade induzida não for em relação a si mesmo
        velocity_ij +=  mac / (4 * np.pi) * \
            ( (ri1j_abs+ri2j_abs) * r12_cross_prod / bound_vortex_den )
    
    return velocity_ij
