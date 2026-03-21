from new.solver.simulation_runner import SimulationResult, SimulationRunner
from new.system.wing_pool import WingPool


class LinearRunner(SimulationRunner):
    """Solves the linearized LLT system once per angle of attack.

    Uses calculate_main_equation_simplified — no iteration, no convergence loop.
    Ignores: damping_factor, max_iter, max_residual, warm_start.
    """

    def run(self, wing_pool: WingPool) -> list[SimulationResult]:
        raise NotImplementedError
