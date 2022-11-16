from dataclasses import dataclass
from models.wingpool import WingPool
from models.flight_condition import FlightCondition
import numpy as np
import numpy.linalg as npla

@dataclass
class Simulation:
    damp: float = 0.7
    max_iter: int = 1000
    max_residual: float = 0.001
    linear_check: bool = False

    def __post_init__(self):
        pass
    
    def calculate_main_equation_simplified(self, wing_pool: WingPool, flight_condition: FlightCondition):
        '''
        TODO: simplified main equation
        '''
        v_inf = flight_condition.v_inf_array

        matrix_dim = sum([wing.N_panels for wing in wing_pool.complete_wing_pool])
        A_matrix = np.zeros([matrix_dim, matrix_dim])
        B_matrix = np.zeros([matrix_dim, 1])
        i_glob = 0
        ind_velocities_dict = wing_pool.ind_velocities_dict
        for wing_i in wing_pool.complete_wing_pool:
                total_velocity_list = wing_pool.total_velocity_dict[wing_i.surface_name]
                G_list_i = wing_pool.G_dict[wing_i.surface_name]
                cp_dsl = wing_i.cp_dsl
                u_n = wing_i.u_n
                for i, _ in enumerate(wing_i.collocation_points):
                    G_i = G_list_i[i]
                    u_n_i = u_n[i]
                    cross_product = 2 * npla.norm(np.cross(v_inf, cp_dsl[i]))
                    Cl_0_i = 1
                    Cl_alpha_i = 1
                    al_0_i = -Cl_0_i/Cl_alpha_i

                    B_matrix[i_glob] = Cl_alpha_i*(np.dot(v_inf,u_n_i)-al_0_i)

                    j_glob = 0
                    i_glob += 1
                    for wing_j in wing_pool.complete_wing_pool:
                        v_ji_distr = ind_velocities_dict[wing_i.surface_name][wing_j.surface_name]
                        for j, _ in enumerate(wing_j.collocation_points):
                            v_ji = v_ji_distr[i][j][:]
                            sum_term = Cl_alpha_i*np.dot(v_ji, u_n_i)
                            A_matrix[i_glob][j_glob] = 1
                            j_glob += 1
        return 0

    def calculate_main_equation(self, wing_pool: WingPool):
        '''
        TODO: solver1 system of equations
        F(G) = R
        '''
        array_dim = sum([wing.N_panels for wing in wing_pool.wing_list])
        R_array = np.zeros([array_dim, 1])
        total = 1
        Cl_i = 1
        for wing in wing_pool.wing_list:
                G_list = wing_pool.G_dict[wing.surface_name]
                total_velocity_list = wing_pool.total_velocity_dict[wing.surface_name]
                for i, _ in enumerate(wing.collocation_points):
                    total_dim_velocity = total
                    norm_value = 1
                    R_array = 2 * norm_value * G_list[i] - Cl_i
        return 0

    def calculate_corrector_equation(self, wing_pool: WingPool, R_array):
        '''
        Newton corrector system of equations
        [J]delta_G = -R
        '''
        matrix_dim = sum([wing.N_panels for wing in wing_pool.complete_wing_pool])
        J_matrix = np.zeros([matrix_dim, matrix_dim])

        i_glob = 0
        for wing_i in wing_pool.complete_wing_pool:
            G_list = wing_pool.G_dict[wing_i.surface_name]
            total_velocity_list = wing_pool.total_velocity_dict[wing_i.surface_name]
            for i, _ in enumerate(wing_i.collocation_points):
                j_glob = 0
                total_dim_velocity = total_velocity_list[i]
                v_ji_distr_i = wing_pool.ind_velocities_dict[wing_i.surface_name]
                G_i = G_list[i]

                w_i = np.cross(total_dim_velocity, wing_i.cp_dsl[i])
                w_i_abs = npla.norm(w_i)
                u_n_i = wing_i.u_n[i]
                u_a_i = wing_i.u_a[i]
                v_n_i = np.dot(total_dim_velocity, wing_i.u_n[i])
                v_a_i = np.dot(total_dim_velocity, wing_i.u_a[i])
                Cl_alpha_i = 1
                w_i_norm = npla.norm(w_i)
                for wing_j in wing_pool.complete_wing_pool:
                    v_ji_distr = v_ji_distr_i[wing_j.surface_name]
                    for j, _ in enumerate(wing_j.collocation_points):
                        v_ji = v_ji_distr[j][:]
                        coef_ij = 2 * np.dot(w_i,np.cross(v_ji, wing_i.cp_dsl[i]))*G_i / w_i_norm \
                            - Cl_alpha_i * (v_a_i * np.dot(v_ji, u_n_i) - v_n_i * np.dot(v_ji, u_a_i)) /(v_a_i ** 2 + v_n_i ** 2)

                        if wing_i.surface_name == wing_j.surface_name and i != j:
                            coef_ij += 2 * w_i_abs
                        
                        J_matrix[i_glob][j_glob] = coef_ij
                        j_glob += 1
                i_glob += 1

        delta_G = npla.solve(J_matrix, -R_array)
        return delta_G
