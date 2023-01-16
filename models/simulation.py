from dataclasses import dataclass, field
from models.wingpool import WingPool
import numpy as np
from utils.timeit import timeit
from simulation.main_equations import (
    calculate_corrector_equation,
    calculate_main_equation,
    calculate_main_equation_simplified
)


@dataclass(repr=False, eq=False, match_args=False)
class Simulation:
    wing_pool: WingPool
    damping_factor: float = 0.7
    max_iter: int = 300
    max_residual: float = 1e-3
    linear_check: bool = False
    matrix_dim: int = field(init=False)


    def __post_init__(self):
        self.matrix_dim = sum([wing.N_panels for wing in self.wing_pool.complete_wing_pool])


    def run_simulation(self) -> list:
        """
        Método para rodar a simulação principal do lltdrek
        A simulação usa o método da convergencia acelerada como default,
        onde o resultado da simulação anterior é o primeiro chute da simulação atual,
        acelerando a convergência
        Outra opção é rodar a equação simplificada como primeiro chute

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
            #     # Solve linear system to get a better approximation for G
            #     G = calculate_main_equation_simplified(
            #         v_inf_array=self.wing_pool.flight_condition.v_inf_list[idx],
            #         aoa_idx=idx,
            #         wing_pool=self.wing_pool,
            #         matrix_dim=self.matrix_dim
            #         )
                G = [0.1 for _ in range(self.matrix_dim)]
                G_dict = self.wing_pool.G_dict
                G_dict = self.wing_pool.update_solution(G_solution=G)
                # Pré-calcular distribuição de alfas e velocidade total por painel
                total_velocity_dict = self.wing_pool.calculate_total_velocity(
                    v_inf_array=self.wing_pool.flight_condition.v_inf_list[idx],
                    G_dict=G_dict
                    )
                aoa_eff_dict = self.wing_pool.calculate_aoa_eff(total_velocity_dict)
            

            if self.linear_check and idx > 0:
                G = calculate_main_equation_simplified(
                    v_inf_array=self.wing_pool.flight_condition.v_inf_list[idx],
                    aoa_idx=idx,
                    wing_pool=self.wing_pool,
                    matrix_dim=self.matrix_dim
                    )
                print(f"Found solution for angle {aoa}")
                G_solution_list.append(G)
            else:
                while True:
                    R_array = calculate_main_equation(
                        total_velocity_dict,
                        aoa_eff_dict,
                        G_dict,
                        self.wing_pool,
                        self.matrix_dim
                    )
                    delta_G = calculate_corrector_equation(
                        R_array,
                        total_velocity_dict,
                        aoa_eff_dict,
                        G_dict,
                        idx,
                        self.wing_pool,
                        self.matrix_dim
                    )

                    G = G + delta_G * self.damping_factor                  
                    G_dict = self.wing_pool.update_solution(G)

                    # Pré-calcular distribuição de alfas e velocidade total por painel
                    total_velocity_dict = self.wing_pool.calculate_total_velocity(
                        v_inf_array=self.wing_pool.flight_condition.v_inf_list[idx],
                        G_dict=G_dict
                        )
                    aoa_eff_dict = self.wing_pool.calculate_aoa_eff(total_velocity_dict)

                    # G_history_list.append(G)
                    iteration += 1
                    if iteration > self.max_iter:
                        G_solution_list.append(np.nan)
                        print(f"Reached max iterations for angle {aoa}")
                        break
                    if abs(R_array.max()) < self.max_residual:
                        G_solution_list.append(self.wing_pool.G_dict)
                        print(f"Found solution for angle {aoa}")
                        print(f"number of iterations: {iteration}")
                        break
            aoa_history_list.append(G_history_list)
        return G_solution_list
