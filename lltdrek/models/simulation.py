from dataclasses import dataclass, field
import numpy as np
from lltdrek.models.wingpool import WingPool
from lltdrek.utils.timeit import timeit
from lltdrek.simulation.main_equations import (
    calculate_corrector_equation,
    calculate_main_equation,
    calculate_main_equation_simplified
)


@dataclass(repr=False, eq=False, match_args=False, slots=True)
class Simulation:

    simulation_modes = [
        "linear_first",
        "latest_solution"
    ]

    wing_pool: WingPool
    damping_factor: float = 0.7
    max_iter: int = 150
    max_residual: float = 1e-3
    linear_check: bool = False
    show_logs: bool = True
    simulation_mode: str = "latest_solution"
    matrix_dim: int = field(init=False)


    def __post_init__(self):
        self.matrix_dim = sum([wing.N_panels for wing in self.wing_pool.complete_wing_pool])
        if self.simulation_mode not in self.simulation_modes:
            raise ValueError(f"Valor de simulation_mode ({self.simulation_mode}) inválido. Valores aceitos: {self.simulation_modes}")


    @timeit
    def run_simulation(self) -> list:
        """
        Método para rodar a simulação principal do lltdrek
        A simulação usa o método da convergencia acelerada como default,
        onde o resultado da simulação anterior é o primeiro chute da simulação atual,
        acelerando a convergência
        """
        print(f"Running simulation for angles between {self.wing_pool.flight_condition.aoa[0]} and {self.wing_pool.flight_condition.aoa[-1]}") if self.show_logs is True else None
        print(f"Linear simulation check: {self.linear_check}") if self.show_logs is True else None
        G_solution_list = []
        for idx, aoa in enumerate(self.wing_pool.flight_condition.aoa):
            iteration = 1

            if idx == 0: 
                G = np.ones(self.matrix_dim) * 0.1
                G_dict = self.wing_pool.G_dict
            if self.simulation_mode == "linear_first":
                # Solve linear system to get a better approximation for G
                G = calculate_main_equation_simplified(
                    v_inf_array=self.wing_pool.flight_condition.v_inf_list[idx],
                    aoa_idx=idx,
                    wing_pool=self.wing_pool,
                    matrix_dim=self.matrix_dim
                    )
                G_linear_solution = G # debugging purposes
                G_dict = self.wing_pool.update_solution(G_solution=G)
                a=1

            total_velocity_dict = self.wing_pool.calculate_total_velocity(
                aoa_idx=idx,
                G_dict=G_dict
                )
            aoa_eff_dict = self.wing_pool.calculate_aoa_eff(total_velocity_dict)
            

            if self.linear_check:
                G = calculate_main_equation_simplified(
                    v_inf_array=self.wing_pool.flight_condition.v_inf_list[idx],
                    aoa_idx=idx,
                    wing_pool=self.wing_pool,
                    matrix_dim=self.matrix_dim
                    )
                G_dict = self.wing_pool.update_solution(G)
                print(f"Found solution for angle {aoa}") if self.show_logs is True else None
                G_solution_list.append(G_dict)
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
                        self.matrix_dim,
                    )

                    if iteration > self.max_iter:
                        G_solution_list.append(np.nan)
                        if "last_successful_solution" in locals():
                            G_dict = last_successful_solution_dict
                        else:
                            G_dict = self.wing_pool.G_dict
                        print(f"Reached max iterations for angle {aoa}") if self.show_logs is True else None
                        break
                    if abs(R_array.max()) < self.max_residual:
                        G_solution_list.append(G_dict)
                        print(f"Found solution for angle {aoa}") if self.show_logs is True else None
                        print(f"number of iterations: {iteration}") if self.show_logs is True else None
                        last_successful_solution_dict = G_dict
                        break
                    else:
                        G = G + delta_G * self.damping_factor                  
                        G_dict = self.wing_pool.update_solution(G)

                        # Pré-calcular distribuição de alfas e velocidade total por painel
                        total_velocity_dict = self.wing_pool.calculate_total_velocity(
                            aoa_idx=idx,
                            G_dict=G_dict
                            )
                        aoa_eff_dict = self.wing_pool.calculate_aoa_eff(total_velocity_dict)
                        iteration += 1

        return G_solution_list
