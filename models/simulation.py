from dataclasses import dataclass, field
from models.wingpool import WingPool
import numpy as np
import numpy.linalg as npla
from utils.lookup import get_airfoil_data, get_linear_data

@dataclass
class Simulation:
    wing_pool: WingPool
    damping_factor: float = 0.7
    max_iter: int = 300
    max_residual: float = 1e-3
    linear_check: bool = False
    matrix_dim: int = field(init=False)


    def __post_init__(self):
        self.matrix_dim = sum([wing.N_panels for wing in self.wing_pool.complete_wing_pool])
    

    def calculate_main_equation_simplified(self, v_inf_array: np.ndarray, ind_velocities_dict: dict) -> np.ndarray:
        A_matrix = np.zeros([self.matrix_dim, self.matrix_dim])
        B_matrix = np.zeros(self.matrix_dim)

        i_glob = 0
        for wing_i in self.wing_pool.complete_wing_pool:
                for i, _ in enumerate(wing_i.collocation_points):
                    A_matrix[i_glob][i_glob] = 2 * npla.norm(np.cross(v_inf_array, wing_i.cp_dsl[i]))
                    linear_data = get_linear_data(wing_i.cp_airfoils[i], wing_i.cp_reynolds[i], wing_i.airfoil_data)
                    Cl_0_i = linear_data["cl0"] * np.pi / 180
                    Cl_alpha_i = linear_data["cl_alpha"] * 180 / np.pi
                    al_0_i = -Cl_0_i/Cl_alpha_i

                    B_matrix[i_glob] = Cl_alpha_i*(np.dot(v_inf_array, wing_i.u_n[i])-al_0_i)

                    j_glob = 0
                    for wing_j in self.wing_pool.complete_wing_pool:
                        v_ij_distr = ind_velocities_dict[wing_i.surface_name][wing_j.surface_name]
                        for j, _ in enumerate(wing_j.collocation_points):
                            v_ij = v_ij_distr[i][j]
                            A_matrix[i_glob][j_glob] += -1*Cl_alpha_i*np.dot(v_ij, wing_i.u_n[i])
                            j_glob += 1
                    i_glob += 1
        G_solution = npla.solve(A_matrix, B_matrix)
        return G_solution


    def calculate_main_equation(self, total_velocity_dict: dict, aoa_eff_dict: dict, G_dict: dict) -> np.ndarray:
        """
        Main system of equations -
        F(G) = R
        """
        R_array = np.zeros(self.matrix_dim)
        for wing in self.wing_pool.complete_wing_pool:
                aoa_eff_distr = aoa_eff_dict[wing.surface_name]
                G_list = G_dict[wing.surface_name]
                total_velocity_distr = total_velocity_dict[wing.surface_name]
                for i, _ in enumerate(wing.collocation_points):
                    Cl_i = 1 * aoa_eff_distr[i] # TODO: Chamar a função de lookup
                    Cl_i = get_airfoil_data(
                        wing.cp_airfoils[i],
                        wing.cp_reynolds[i],
                        aoa_eff_distr[i],
                        wing.airfoil_data,
                        cl_alpha_check = False
                    )
                    norm_value = npla.norm(np.cross(total_velocity_distr[i], wing.cp_dsl[i]))
                    R_array[i] = 2 * norm_value * G_list[i] - Cl_i
        return R_array


    def calculate_corrector_equation(
        self,
        R_array: np.ndarray,
        total_velocity_dict: dict,
        aoa_eff_dict: dict,
        G_dict: dict,
        ind_velocity_dict: dict
        ) -> np.ndarray:
        """
        Newton corrector system of equations
        [J]delta_G = -R
        """
        J_matrix = np.zeros([self.matrix_dim, self.matrix_dim])

        i_glob = 0
        for wing_i in self.wing_pool.complete_wing_pool:
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
                        aoa_eff_distr[i],
                        wing_i.airfoil_data,
                        cl_alpha_check = True
                    ) * 180 / np.pi
                w_i_norm = npla.norm(w_i)
                for wing_j in self.wing_pool.complete_wing_pool:
                    v_ij_distr = ind_velocity_dict[wing_i.surface_name][wing_j.surface_name]
                    for j, _ in enumerate(wing_j.collocation_points):
                        v_ij = v_ij_distr[i][j] # TODO: verificar se índices estão corretos
                        coef_ij = 2 * np.dot(w_i,np.cross(v_ij, wing_i.cp_dsl[i]))*G_distr[i] / w_i_norm \
                            - Cl_alpha_i * (v_a_i * np.dot(v_ij, u_n_i) - v_n_i * np.dot(v_ij, u_a_i)) /(v_a_i ** 2 + v_n_i ** 2)

                        if wing_i.surface_name == wing_j.surface_name and i != j:
                            coef_ij += 2 * w_i_abs
                        
                        J_matrix[i_glob][j_glob] = coef_ij
                        j_glob += 1
                i_glob += 1

        delta_G = npla.solve(J_matrix, -R_array)
        return delta_G


    def run_simulation(self) -> list:
        """
        Método para rodar a simulação principal do lltdrek
        O chute inicial da lista G já está presente na wing_pool
        A simulação usa o método da convergencia acelerada naturalmente,
        onde o resultado da simulação anterior é o primeiro chute da simulação atual,
        acelerando a convergência

        tem bastante coisa pra aprimorar aqui nesse método
        """
        print(f"Running simulation for angles between {self.wing_pool.flight_condition.aoa[0]} and {self.wing_pool.flight_condition.aoa[-1]}")
        print(f"Linear simulation check: {self.linear_check}")
        G_solution_list = []
        aoa_history_list = []
        for idx, aoa in enumerate(self.wing_pool.flight_condition.aoa):
            iteration = 0
            G_history_list = []

            if idx == 0:
                # Solve linear system to get a better approximation for G
                G = self.calculate_main_equation_simplified(
                    v_inf_array=self.wing_pool.flight_condition.v_inf_list[idx],
                    ind_velocities_dict=self.wing_pool.ind_velocities_list[idx],
                    )
                G_dict = self.wing_pool.update_solution(G_solution=G)


            if self.linear_check:
                print(f"Found solution for angle {aoa}")
                G_solution_list.append(G)
            else:
                while True:
                    # Pré-calcular distribuição de alfas e velocidade total por painel
                    total_velocity_dict = self.wing_pool.calculate_total_velocity(
                        v_inf_array=self.wing_pool.flight_condition.v_inf_list[idx],
                        G_dict=G_dict
                        )
                    aoa_eff_dict = self.wing_pool.calculate_aoa_eff(total_velocity_dict)

                    R_array = self.calculate_main_equation(total_velocity_dict, aoa_eff_dict, G_dict)
                    delta_G = self.calculate_corrector_equation(
                        R_array,
                        total_velocity_dict,
                        aoa_eff_dict,
                        G_dict,
                        self.wing_pool.ind_velocities_list[idx]
                    )

                    G = G + delta_G * self.damping_factor                  
                    G_dict = self.wing_pool.update_solution(G)

                    G_history_list.append(G)
                    iteration += 1
                    if iteration > self.max_iter:
                        G_solution_list.append(np.nan)
                        print(f"Reached max iterations for angle {aoa}")
                        break
                    if abs(R_array.max()) < self.max_residual:
                        G_solution_list.append(self.wing_pool.G_dict)
                        print(f"Found solution for angle {aoa}")
                        break
            aoa_history_list.append(G_history_list)
        return G_solution_list
