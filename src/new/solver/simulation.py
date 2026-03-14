from new.solver.linear_runner import LinearRunner
from new.solver.nonlinear_loops_runner import NonlinearLoopsRunner
from new.solver.nonlinear_numpy_runner import NonlinearNumpyRunner
from new.solver.simulation_runner import SimulationResult, SimulationRunner
from new.system.wing_pool import WingPool

_ALLOWED_EQUATIONS = ("linear", "nonlinear")
_ALLOWED_IMPLEMENTATIONS = ("loops", "numpy")


class Simulation:
    """User-facing entry point for running LLT simulations.

    Selects the appropriate SimulationRunner based on three orthogonal params:
    - equations: which mathematical model to use ("linear" | "nonlinear")
    - implementation: computational backend for nonlinear solves ("loops" | "numpy")
    - warm_start: for nonlinear only — use linear solution as initial G each alpha

    Example usage:
        sim = Simulation(equations="nonlinear", implementation="loops", warm_start=True)
        results = sim.run(wing_pool)
    """

    def __init__(
        self,
        equations: str = "nonlinear",
        implementation: str = "loops",
        warm_start: bool = False,
        damping_factor: float = 0.7,
        max_iter: int = 150,
        max_residual: float = 1e-3,
        show_logs: bool = True,
    ) -> None:
        if equations not in _ALLOWED_EQUATIONS:
            raise ValueError(f"Invalid equations '{equations}'. Choose from {_ALLOWED_EQUATIONS}.")
        if implementation not in _ALLOWED_IMPLEMENTATIONS:
            raise ValueError(f"Invalid implementation '{implementation}'. Choose from {_ALLOWED_IMPLEMENTATIONS}.")

        self.equations = equations
        self.implementation = implementation
        self.warm_start = warm_start
        self.damping_factor = damping_factor
        self.max_iter = max_iter
        self.max_residual = max_residual
        self.show_logs = show_logs

        self._runner: SimulationRunner = self._build_runner()

    def _build_runner(self) -> SimulationRunner:
        kwargs = dict(
            damping_factor=self.damping_factor,
            max_iter=self.max_iter,
            max_residual=self.max_residual,
            show_logs=self.show_logs,
            warm_start=self.warm_start,
        )
        if self.equations == "linear":
            return LinearRunner(**kwargs)
        if self.implementation == "numpy":
            return NonlinearNumpyRunner(**kwargs)
        return NonlinearLoopsRunner(**kwargs)

    def run(self, wing_pool: WingPool) -> list[SimulationResult]:
        """Run the simulation over all angles of attack defined in wing_pool.flight_condition."""
        return self._runner.run(wing_pool)
