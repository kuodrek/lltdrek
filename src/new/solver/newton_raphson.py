from abc import ABC, abstractmethod

import numpy as np


class NewtonRaphsonMixin(ABC):
    """Domain-agnostic Newton-Raphson iteration engine (mixin).

    Provides a concrete _run_newton_raphson() template method and four abstract
    hooks that the consuming class must implement. The mixin knows nothing about
    WingPool, alpha, or aerodynamics — all domain wiring lives in the runner.

    Relies on the consuming class having:
        - self.max_iter (int)
        - self.max_residual (float)
        - self.damping_factor (float)

    These are provided by SimulationRunner.__init__.
    """

    @abstractmethod
    def _init_state(self) -> None:
        """Set up iteration state (e.g. self.state) before the loop starts."""
        ...

    @abstractmethod
    def _update_state(self, delta: np.ndarray) -> None:
        """Apply the damped Newton correction and rebuild any derived state."""
        ...

    @abstractmethod
    def calculate_residual(self) -> np.ndarray:
        """Compute R(G). Reads current state from self."""
        ...

    @abstractmethod
    def calculate_correction(self, residual: np.ndarray) -> np.ndarray:
        """Compute delta_G given the current residual. Reads current state from self."""
        ...

    def _run_newton_raphson(self) -> tuple[np.ndarray, bool]:
        """Run the Newton-Raphson loop for a single solve.

        Returns:
            (residual, converged) — the final residual array and whether the
            solver converged within max_iter iterations.
        """
        self._init_state()
        R = self.calculate_residual()
        for _ in range(self.max_iter):
            if np.max(np.abs(R)) < self.max_residual:
                return R, True
            delta = self.calculate_correction(R)
            self._update_state(delta)
            R = self.calculate_residual()
        return R, False
