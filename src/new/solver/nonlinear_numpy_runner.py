from new.solver.simulation_runner import SimulationResult, SimulationRunner
from new.system.wing_pool import WingPool


class NonlinearNumpyRunner(SimulationRunner):
    """Newton-Raphson solver using fully vectorized numpy operations.

    Implements the same equations as NonlinearLoopsRunner with no Python-level
    loops over panels. Validate output against NonlinearLoopsRunner before use.

    TODO: implement vectorized matrix assembly and Newton-Raphson loop.
    """

    def run(self, wing_pool: WingPool) -> list[SimulationResult]:
        raise NotImplementedError
