from dataclasses import dataclass

import numpy as np

from new.solver.newton_raphson import NewtonRaphsonMixin
from new.solver.simulation_runner import SimulationResult, SimulationRunner
from new.system.wing_pool import WingPool


@dataclass
class IterationState:
    """Per-iteration state for the Newton-Raphson loop.

    Grouping these here makes it easy to distinguish:
        self.state.X      — iteration state (changes every iteration)
        self._wing_pool   — domain context (set once per run() call)
        self.damping_factor / max_iter / ... — runner config (set at construction)
    """

    G: np.ndarray
    G_dict: dict[str, np.ndarray]
    total_velocity_dict: dict[str, np.ndarray]
    aoa_eff_dict: dict[str, np.ndarray]


class NonlinearLoopsRunner(SimulationRunner, NewtonRaphsonMixin):
    """Newton-Raphson solver using Python for-loops over panels.

    This is the validated baseline implementation, ported directly from the
    legacy Simulation.run(). Used as the reference when benchmarking the
    numpy-vectorized runner.

    warm_start=True: solves the linear system first to get an initial G per alpha.
    warm_start=False: starts from G = 0.1 * ones, reuses previous alpha solution after first.
    """

    def run(self, wing_pool: WingPool) -> list[SimulationResult]:
        self._wing_pool = wing_pool
        results = []
        for alpha in wing_pool.flight_condition.alphas:
            self._alpha = alpha
            R, converged = self._run_newton_raphson()
            results.append(SimulationResult(alpha=alpha, G=self.state.G_dict, residual=R, converged=converged))
        return results

    # --- NewtonRaphsonMixin hooks ---

    def _init_state(self) -> None:
        G = np.ones(self._wing_pool.matrix_dim) * 0.1  # TODO: warm_start logic
        G_dict = self._wing_pool.map_solution(G)
        total_velocity_dict = self._wing_pool.calculate_total_velocity(self._alpha, G_dict)
        aoa_eff_dict = self._wing_pool.calculate_aoa_eff(total_velocity_dict)
        self.state = IterationState(G, G_dict, total_velocity_dict, aoa_eff_dict)

    def _update_state(self, delta: np.ndarray) -> None:
        G = self.state.G + self.damping_factor * delta
        G_dict = self._wing_pool.map_solution(G)
        total_velocity_dict = self._wing_pool.calculate_total_velocity(self._alpha, G_dict)
        aoa_eff_dict = self._wing_pool.calculate_aoa_eff(total_velocity_dict)
        self.state = IterationState(G, G_dict, total_velocity_dict, aoa_eff_dict)

    def calculate_residual(self) -> np.ndarray:
        return self._compute_residual(
            self.state.G_dict,
            self.state.total_velocity_dict,
            self.state.aoa_eff_dict,
            self._wing_pool,
        )

    def _compute_residual(
        self,
        _G_dict: dict[str, np.ndarray],
        _total_velocity_dict: dict[str, np.ndarray],
        _aoa_eff_dict: dict[str, np.ndarray],
        _wing_pool: WingPool,
    ) -> np.ndarray:
        # TODO: port of legacy calculate_main_equation (for-loops over panels)
        # Pure function of its inputs — unit-testable without runner setup
        raise NotImplementedError

    def calculate_correction(self, residual: np.ndarray) -> np.ndarray:
        return self._compute_correction(
            residual,
            self.state.G_dict,
            self.state.total_velocity_dict,
            self.state.aoa_eff_dict,
            self._wing_pool,
            self._alpha,
        )

    def _compute_correction(
        self,
        residual: np.ndarray,
        G_dict: dict[str, np.ndarray],
        total_velocity_dict: dict[str, np.ndarray],
        aoa_eff_dict: dict[str, np.ndarray],
        wing_pool: WingPool,
        alpha: float,
    ) -> np.ndarray:
        # TODO: port of legacy calculate_corrector_equation (builds Jacobian, np.linalg.solve)
        # Pure function of its inputs — unit-testable without runner setup
        raise NotImplementedError
