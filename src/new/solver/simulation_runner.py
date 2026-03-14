from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from new.system.wing_pool import WingPool


@dataclass
class SimulationResult:
    """Result of a simulation for a single angle of attack.

    :param alpha: Angle of attack.
    :param G_solution: Dimensionless vortex strength per panel, keyed by surface name.
    :param residual: Residual array at convergence (or at max_iter).
    :param convergence_check: True if the solver converged within max_iter.
    """

    alpha: float
    G_solution: dict
    residual: np.ndarray
    convergence_check: bool


class SimulationRunner(ABC):
    """Abstract base class for simulation strategies.

    All concrete runners share the same __init__ signature so that
    Simulation can construct any runner uniformly. Unused params are ignored.
    """

    def __init__(
        self,
        damping_factor: float,
        max_iter: int,
        max_residual: float,
        show_logs: bool,
        warm_start: bool,
    ) -> None:
        self.damping_factor = damping_factor
        self.max_iter = max_iter
        self.max_residual = max_residual
        self.show_logs = show_logs
        self.warm_start = warm_start

    @abstractmethod
    def run(self, wing_pool: WingPool) -> list[SimulationResult]:
        """Execute the simulation over all angles of attack in wing_pool.flight_condition."""
        ...
