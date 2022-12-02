from dataclasses import dataclass, field
from models.wingpool import WingPool
import numpy as np
import numpy.linalg as npla


@dataclass
class Simulation:
    wing_pool: WingPool
    damping_factor: float = 0.7
    max_iter: int = 1000
    max_residual: float = 1e-3
    linear_check: bool = False
    matrix_dim: int = field(init=False)

    def __post_init__(self):
        self.matrix_dim = sum([wing.N_panels for wing in self.wing_pool.complete_wing_pool])
    
    def calculate_main_equation_simplified(self, v_inf_array: np.ndarray, ind_velocities_dict: dict, aoa_eff_dict: dict) -> np.ndarray:
        A_matrix = np.zeros([self.matrix_dim, self.matrix_dim])
        B_matrix = np.zeros([self.matrix_dim, 1])

        i_glob = 0
        for wing_i in self.wing_pool.complete_wing_pool:
                G_list_i = self.wing_pool.G_dict[wing_i.surface_name]
                aoa_eff_distr = aoa_eff_dict[wing_i.surface_name]
                cp_dsl = wing_i.cp_dsl
                u_n = wing_i.u_n
                for i, _ in enumerate(wing_i.collocation_points):
                    G_i = G_list_i[i]
                    u_n_i = u_n[i]
                    A_matrix[i_glob][i_glob] = 2 * npla.norm(np.cross(v_inf_array, cp_dsl[i]))
                    Cl_0_i = 1
                    Cl_alpha_i = 1 * aoa_eff_distr[i]
                    al_0_i = -Cl_0_i/Cl_alpha_i

                    B_matrix[i_glob] = Cl_alpha_i*(np.dot(v_inf_array,u_n_i)-al_0_i)

                    j_glob = 0
                    i_glob += 1
                    for wing_j in self.wing_pool.complete_wing_pool:
                        v_ji_distr = ind_velocities_dict[wing_i.surface_name][wing_j.surface_name]
                        for j, _ in enumerate(wing_j.collocation_points):
                            v_ji = v_ji_distr[i][j][:]
                            A_matrix[i_glob][j_glob] += -1*Cl_alpha_i*np.dot(v_ji, u_n_i)
                            j_glob += 1
        G_solution = npla.solve(A_matrix, B_matrix)
        return G_solution

    def calculate_main_equation(self, total_velocity_dict: dict, aoa_eff_dict: dict) -> np.ndarray:
        '''
        Main system of equations -
        F(G) = R
        '''
        array_dim = sum([wing.N_panels for wing in self.wing_pool.wing_list])
        R_array = np.zeros([array_dim, 1])
        for wing in self.wing_pool.wing_list:
                aoa_eff_distr = aoa_eff_dict[wing.surface_name]
                G_list = self.wing_pool.G_dict[wing.surface_name]
                total_velocity_distr = total_velocity_dict[wing.surface_name]
                for i, _ in enumerate(wing.collocation_points):
                    Cl_i = 1 * aoa_eff_distr[i] # WIP
                    norm_value = np.absolute(np.cross(total_velocity_distr[i], wing.cp_dsl[i]))
                    R_array[i] = 2 * norm_value * G_list[i] - Cl_i
        return R_array

    def calculate_corrector_equation(self, R_array: np.ndarray) -> np.ndarray:
        '''
        Newton corrector system of equations
        [J]delta_G = -R
        '''
        J_matrix = np.zeros([self.matrix_dim, self.matrix_dim])

        i_glob = 0
        for wing_i in self.wing_pool.complete_wing_pool:
            G_list = self.wing_pool.G_dict[wing_i.surface_name]
            total_velocity_list = self.wing_pool.total_velocity_dict[wing_i.surface_name]
            for i, _ in enumerate(wing_i.collocation_points):
                j_glob = 0
                total_dim_velocity = total_velocity_list[i]
                v_ji_distr_i = self.wing_pool.ind_velocities_dict[wing_i.surface_name]
                G_i = G_list[i]

                w_i = np.cross(total_dim_velocity, wing_i.cp_dsl[i])
                w_i_abs = npla.norm(w_i)
                u_n_i = wing_i.u_n[i]
                u_a_i = wing_i.u_a[i]
                v_n_i = np.dot(total_dim_velocity, wing_i.u_n[i])
                v_a_i = np.dot(total_dim_velocity, wing_i.u_a[i])
                Cl_alpha_i = 1
                w_i_norm = npla.norm(w_i)
                for wing_j in self.wing_pool.complete_wing_pool:
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

    def run_simulation(self):
        '''
        Método para rodar a simulação principal do lltdrek
        O chute inicial da lista G já está presente na wing_pool
        A simulação usa o método da convergencia acelerada naturalmente,
        onde o resultado da simulação anterior é o primeiro chute da simulação atual,
        acelerando a convergência
        '''
        iteration = 0
        G_solution_list = []
        for aoa in self.wing_pool.flight_condition.aoa:
            G_list_history = []
            while True:
                R_array = self.calculate_main_equation()
                delta_G = self.calculate_corrector_equation(R_array)

                G = G + delta_G * self.damping_factor
                G_list_history.append(G)

                self.wing_pool.update_solution(G)
                # É preciso chamar os métodos que atualizam as velocidades totais e aoa_eff
                aoa_eff_dict = self.wing_pool.calculate_aoa_eff()
                total_velocity_dict = self.wing_pool.calculate_total_velocity()

                iteration += 1
                if iteration > self.max_iter:
                    G_solution_list.append(np.nan)
                    break
                if R_array.max < self.max_residual:
                    G_solution_list.append(self.wing_pool.G_dict)
                    break
