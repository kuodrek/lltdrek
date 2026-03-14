from new.solver.simulation_runner import SimulationResult, SimulationRunner
from new.system.wing_pool import WingPool


class NonlinearLoopsRunner(SimulationRunner):
    """Newton-Raphson solver using Python for-loops over panels.

    This is the validated baseline implementation, ported directly from the
    legacy Simulation.run(). Used as the reference when benchmarking the
    numpy-vectorized runner.

    warm_start=True: solves the linear system first to get an initial G per alpha.
    warm_start=False: starts from G = 0.1 * ones, reuses previous alpha solution after first.
    """

    def run(self, wing_pool: WingPool) -> list[SimulationResult]:
        raise NotImplementedError
