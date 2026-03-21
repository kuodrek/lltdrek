import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass

import numpy as np

from new.system.wing_pool import WingPool


@dataclass
class SimulationResult:
    """Result of a simulation for a single angle of attack.

    :param alpha: Angle of attack.
    :param G: Dimensionless vortex strength per panel, keyed by surface name.
    :param residual: Residual array at convergence (or at max_iter).
    :param converged: True if the solver converged within max_iter.
    """

    alpha: float
    G: dict[str, np.ndarray]
    residual: np.ndarray
    converged: bool


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
        warm_start: bool,
    ) -> None:
        self.damping_factor = damping_factor
        self.max_iter = max_iter
        self.max_residual = max_residual
        self.warm_start = warm_start
        self.logger = logging.getLogger(f"lltdrek.solver.{type(self).__name__}")

    @abstractmethod
    def run(self, wing_pool: WingPool) -> list[SimulationResult]:
        """Execute the simulation over all angles of attack in wing_pool.flight_condition."""
        ...
