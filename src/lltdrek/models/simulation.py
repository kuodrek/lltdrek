from enum import Enum
from dataclasses import dataclass
import numpy as np
from lltdrek.models.wingpool import WingPool
from lltdrek.simulation.main_equations import (
    calculate_corrector_equation,
    calculate_main_equation,
    calculate_main_equation_simplified
)
from lltdrek.models.types import AngleOfAttack, DVSMap

class SimulationModes(Enum):
    """Ways of running a simulation
    
    - `LINEAR_FIRST`: For each alpha in a non-linear simulation, solve the linear set of equations first
    as a first guess for the non-linear problem
    - `LATEST_SOLUTION`: Use the solution of previous alpha as a first guess for the next alpha

    If simulation is linear, the mode doesn't matter
    """
    LINEAR_FIRST = "linear_first"
    LATEST_SOLUTION = "latest_solution"


@dataclass(repr=False, eq=False, match_args=False, slots=True)
class SimulationResult:
    """
    :param alpha: Angle of attack of simulation
    :type alpha: AngleOfAttack
    :param G_solution: Dictionary of panel's dimensionless vortex strength.
    Each key is a surface from related `WingPool`
    :type G_solution: DSVMap
    :param residual: Residual array of simulation
    :param convergence_check: True if simulation converged
    """
    alpha: AngleOfAttack
    G_solution: DVSMap
    residual: np.ndarray
    convergence_check: bool

    def __repr__(self):
        return f"""SimulationResult(alpha={self.alpha}, convergence_check={self.convergence_check})"""


@dataclass(repr=False, eq=False, match_args=False, slots=True)
class Simulation:
    simulation_modes = [
        e.value for e in SimulationModes
    ]  # Allowed values for simulation_modes

    damping_factor: float = 0.7
    max_iter: int = 150
    max_residual: float = 1e-3
    linear_check: bool = False
    show_logs: bool = True
    simulation_mode: str = SimulationModes.LATEST_SOLUTION.value

    def __post_init__(self):
        if self.simulation_mode not in self.simulation_modes:
            raise ValueError(f"Valor de simulation_mode ({self.simulation_mode}) inválido. Valores aceitos: {self.simulation_modes}")

    @classmethod
    def _get_matrix_dimension(cls, wing_pool: WingPool):
        return sum([wing.N_panels for wing in wing_pool.pool])

    # @timeit
    def run(self, wing_pool: WingPool) -> list[SimulationResult]:
        """
        Método para rodar a simulação principal do lltdrek
        A simulação usa o método da convergencia acelerada como default,
        onde o resultado da simulação anterior é o primeiro chute da simulação atual,
        acelerando a convergência
        """
        matrix_dim = self._get_matrix_dimension(wing_pool)

        print(f"Running simulation for angles between {wing_pool.flight_condition.angles_of_attack[0]} and {wing_pool.flight_condition.angles_of_attack[-1]}") if self.show_logs is True else None
        print(f"Linear simulation check: {self.linear_check}") if self.show_logs is True else None
        G_solution_list = []
        for idx, alpha in enumerate(wing_pool.flight_condition.angles_of_attack):
            iteration = 1
            freestream_velocities = wing_pool.system_freestream_velocities[alpha]
            if idx == 0: 
                G = np.ones(matrix_dim) * 0.1
                G_dict = wing_pool.G_dict
            if self.simulation_mode == "linear_first":
                # Solve linear system to get a better approximation for G
                G = calculate_main_equation_simplified(
                    freestream_velocities=freestream_velocities,
                    alpha=alpha,
                    wing_pool=wing_pool,
                    matrix_dim=matrix_dim,
                    show_logs=self.show_logs
                )
                G_dict = wing_pool.map_solution(G)

            total_velocity_dict = wing_pool.calculate_total_velocity(
                    alpha=alpha,
                    G_dict=G_dict
                )
            aoa_eff_dict = wing_pool.calculate_aoa_eff(total_velocity_dict)

            while True:
                R_array = calculate_main_equation(
                    total_velocity_dict,
                    aoa_eff_dict,
                    G_dict,
                    wing_pool,
                    matrix_dim,
                    self.linear_check,
                    self.show_logs
                )
                delta_G = calculate_corrector_equation(
                    R_array,
                    total_velocity_dict,
                    aoa_eff_dict,
                    G_dict,
                    alpha,
                    wing_pool,
                    matrix_dim,
                    self.linear_check,
                    self.show_logs
                )
                if iteration > self.max_iter:
                    G_solution = np.ones(matrix_dim) * np.nan
                    G_dict = wing_pool.map_solution(G=G_solution)
                    G_solution_list.append(SimulationResult(
                        alpha,
                        G_dict,
                        R_array,
                        convergence_check=False
                    ))
                    if "last_successful_solution" in locals():
                        G_dict = last_successful_solution_dict
                    else:
                        G_dict = wing_pool.G_dict
                    print(f"Reached max iterations for angle {alpha}") if self.show_logs is True else None
                    break
                if abs(R_array.max()) < self.max_residual:
                    G_solution_list.append(SimulationResult(
                        alpha,
                        G_dict,
                        R_array,
                        convergence_check=True
                    ))
                    print(f"Found solution for angle {alpha}") if self.show_logs is True else None
                    print(f"number of iterations: {iteration}") if self.show_logs is True else None
                    last_successful_solution_dict = G_dict
                    break
                else:
                    G = G + delta_G * self.damping_factor                  
                    G_dict = wing_pool.map_solution(G)

                    # Pre calculate alpha distribution and total velocity for each panel
                    total_velocity_dict = wing_pool.calculate_total_velocity(
                        alpha=alpha,
                        G_dict=G_dict
                        )
                    aoa_eff_dict = wing_pool.calculate_aoa_eff(total_velocity_dict)
                    # print(f"aoa_eff_dict: {aoa_eff_dict['asa'] * 180 / np.pi}")
                    iteration += 1

        return G_solution_list
