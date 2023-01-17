import numpy as np
import numpy.linalg as npla
from models.wingpool import WingPool
from utils.lookup import get_airfoil_data, get_linear_data
from utils.timeit import timeit


def calculate_main_equation_simplified(
    v_inf_array: np.ndarray,
    aoa_idx: int,
    wing_pool: WingPool,
    matrix_dim: int
) -> np.ndarray:
    A_matrix = np.zeros([matrix_dim, matrix_dim])
    B_matrix = np.zeros(matrix_dim)

    i_glob = 0
    for wing_i in wing_pool.complete_wing_pool:
            for i, _ in enumerate(wing_i.collocation_points):
                A_matrix[i_glob][i_glob] = 2 * npla.norm(np.cross(v_inf_array, wing_i.cp_dsl[i]))
                linear_data = get_linear_data(wing_i.cp_airfoils[i], wing_i.cp_reynolds[i], wing_i.airfoil_data)
                Cl_0_i = linear_data["cl0"] * np.pi / 180
                Cl_alpha_i = linear_data["cl_alpha"] * 180 / np.pi
                al_0_i = -Cl_0_i/Cl_alpha_i

                B_matrix[i_glob] = Cl_alpha_i*(np.dot(v_inf_array, wing_i.u_n[i])-al_0_i)

                j_glob = 0
                for wing_j in wing_pool.complete_wing_pool:
                    v_ij_distr = wing_pool.ind_velocities_list[aoa_idx][wing_i.surface_name][wing_j.surface_name]
                    for j, _ in enumerate(wing_j.collocation_points):
                        v_ij = v_ij_distr[i][j]
                        A_matrix[i_glob][j_glob] += -1*Cl_alpha_i*np.dot(v_ij, wing_i.u_n[i])
                        j_glob += 1
                i_glob += 1
    G_solution = npla.solve(A_matrix, B_matrix)
    return G_solution


def calculate_main_equation(
    total_velocity_dict: dict,
    aoa_eff_dict: dict,
    G_dict: dict,
    wing_pool: WingPool,
    matrix_dim: int
) -> np.ndarray:
    """
    Main system of equations -
    F(G) = R
    """
    R_array = np.zeros(matrix_dim)
    i_glob = 0
    for wing in wing_pool.complete_wing_pool:
        if "_mirrored" in wing.surface_name: continue
        N_panels = wing.N_panels
        aoa_eff_distr = aoa_eff_dict[wing.surface_name]
        G_list = G_dict[wing.surface_name]
        total_velocity_distr = total_velocity_dict[wing.surface_name]
        for i, _ in enumerate(wing.collocation_points):
            Cl_i = get_airfoil_data(
                wing.cp_airfoils[i],
                wing.cp_reynolds[i],
                aoa_eff_distr[i] * 180 / np.pi,
                wing.airfoil_data,
                cl_alpha_check = False
            )
            R_array[i_glob] = 2 * npla.norm(np.cross(total_velocity_distr[i], wing.cp_dsl[i])) * G_list[i] \
                - Cl_i
            R_array[N_panels+i_glob] = R_array[i_glob]
            i_glob += 1
    return R_array


# @timeit
def calculate_corrector_equation(
    R_array: np.ndarray,
    total_velocity_dict: dict,
    aoa_eff_dict: dict,
    G_dict: dict,
    aoa_idx: int,
    wing_pool: WingPool,
    matrix_dim: int
) -> np.ndarray:
    """
    Newton corrector system of equations
    [J]delta_G = -R
    """
    J_matrix = np.zeros([matrix_dim, matrix_dim])
    i_glob = 0
    for wing_i in wing_pool.complete_wing_pool:
        G_distr = G_dict[wing_i.surface_name]
        total_velocity_distr = total_velocity_dict[wing_i.surface_name]
        aoa_eff_distr = aoa_eff_dict[wing_i.surface_name]
        for i, _ in enumerate(wing_i.collocation_points):
            N_panels = wing_i.N_panels
            j_glob = 0
            w_i = np.cross(total_velocity_distr[i], wing_i.cp_dsl[i])
            w_i_abs = npla.norm(w_i)
            u_n_i = wing_i.u_n[i]
            u_a_i = wing_i.u_a[i]
            v_n_i = np.dot(total_velocity_distr[i], wing_i.u_n[i])
            v_a_i = np.dot(total_velocity_distr[i], wing_i.u_a[i])
            Cl_alpha_i = get_airfoil_data(
                    wing_i.cp_airfoils[i],
                    wing_i.cp_reynolds[i],
                    aoa_eff_distr[i] * 180 / np.pi,
                    wing_i.airfoil_data,
                    cl_alpha_check = True
                ) * 180 / np.pi
            for wing_j in wing_pool.complete_wing_pool:
                if "_mirrored" in wing_j.surface_name: continue
                v_ij_distr = wing_pool.ind_velocities_list[aoa_idx][wing_i.surface_name][wing_j.surface_name]
                v_ij_distr_mirrored = wing_pool.ind_velocities_list[aoa_idx][wing_i.surface_name][wing_j.surface_name+"_mirrored"]
                for j, _ in enumerate(wing_j.collocation_points):
                    v_ij = v_ij_distr[i][j]
                    v_ij_m = v_ij_distr_mirrored[i][j]

                    coef_ij = 2 * np.dot(w_i,np.cross(v_ij, wing_i.cp_dsl[i]))*G_distr[i] / w_i_abs \
                        - Cl_alpha_i * (v_a_i * np.dot(v_ij, u_n_i) - v_n_i * np.dot(v_ij, u_a_i)) /(v_a_i ** 2 + v_n_i ** 2)
                    
                    coef_ij_m = 2 * np.dot(w_i,np.cross(v_ij_m, wing_i.cp_dsl[i]))*G_distr[i] / w_i_abs \
                        - Cl_alpha_i * (v_a_i * np.dot(v_ij_m, u_n_i) - v_n_i * np.dot(v_ij_m, u_a_i)) /(v_a_i ** 2 + v_n_i ** 2)
    
                    if wing_i.parent_wing == wing_j.surface_name and i_glob == (N_panels+j_glob):
                        coef_ij_m += 2 * w_i_abs
                        if "_mirrored" in wing_i.surface_name: # Essa parcela está sendo desconsiderada no llt original. Por quê?
                            coef_ij = 0

                    if wing_i.surface_name == wing_j.surface_name and i == j:
                        coef_ij += 2 * w_i_abs
                        coef_ij_m = 0
                    
                    J_matrix[i_glob][j_glob] = coef_ij
                    J_matrix[i_glob][N_panels+j_glob] = coef_ij_m
                    j_glob += 1
            i_glob += 1

    delta_G = npla.solve(J_matrix, -R_array)
    return delta_G


def calculate_corrector_equation_backup(
    R_array: np.ndarray,
    total_velocity_dict: dict,
    aoa_eff_dict: dict,
    G_dict: dict,
    aoa_idx: int,
    wing_pool: WingPool,
    matrix_dim: int
) -> np.ndarray:
    """
    Newton corrector system of equations
    [J]delta_G = -R
    """
    J_matrix = np.zeros([matrix_dim, matrix_dim])
    i_glob = 0
    for wing_i in wing_pool.complete_wing_pool:
        G_distr = G_dict[wing_i.surface_name]
        total_velocity_distr = total_velocity_dict[wing_i.surface_name]
        aoa_eff_distr = aoa_eff_dict[wing_i.surface_name]
        for i, _ in enumerate(wing_i.collocation_points):
            j_glob = 0
            w_i = np.cross(total_velocity_distr[i], wing_i.cp_dsl[i])
            w_i_abs = npla.norm(w_i)
            u_n_i = wing_i.u_n[i]
            u_a_i = wing_i.u_a[i]
            v_n_i = np.dot(total_velocity_distr[i], wing_i.u_n[i])
            v_a_i = np.dot(total_velocity_distr[i], wing_i.u_a[i])
            Cl_alpha_i = get_airfoil_data(
                    wing_i.cp_airfoils[i],
                    wing_i.cp_reynolds[i],
                    aoa_eff_distr[i] * 180 / np.pi,
                    wing_i.airfoil_data,
                    cl_alpha_check = True
                ) * 180 / np.pi
            for wing_j in wing_pool.complete_wing_pool:
                v_ij_distr = wing_pool.ind_velocities_list[aoa_idx][wing_i.surface_name][wing_j.surface_name]
                for j, _ in enumerate(wing_j.collocation_points):
                    v_ij = v_ij_distr[i][j]
                    coef_ij = 2 * np.dot(w_i,np.cross(v_ij, wing_i.cp_dsl[i]))*G_distr[i] / w_i_abs \
                        - Cl_alpha_i * (v_a_i * np.dot(v_ij, u_n_i) - v_n_i * np.dot(v_ij, u_a_i)) /(v_a_i ** 2 + v_n_i ** 2)

                    if wing_i.surface_name == wing_j.surface_name and i == j:
                        coef_ij += 2 * w_i_abs
                    
                    J_matrix[i_glob][j_glob] = coef_ij
                    j_glob += 1
            i_glob += 1

    delta_G = npla.solve(J_matrix, -R_array)
    return delta_G
