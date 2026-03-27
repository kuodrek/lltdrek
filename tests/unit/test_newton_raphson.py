"""Contract tests for NewtonRaphsonMixin.

Tests the iteration loop behavior using a FakeRunner with canned residuals and
corrections. No aerodynamics — these tests verify the mixin wiring only.
"""

import numpy as np
import pytest

from new.solver.newton_raphson import NewtonRaphsonMixin
from new.solver.simulation_runner import SimulationResult, SimulationRunner
from new.system.wing_pool import WingPool


class FakeRunner(SimulationRunner, NewtonRaphsonMixin):
    """Minimal concrete runner for testing NewtonRaphsonMixin contracts.

    Consumes pre-set _residuals and _corrections lists in order.
    Records call order in _call_log for hook-ordering assertions.
    """

    def run(self, wing_pool: WingPool) -> list[SimulationResult]:
        raise NotImplementedError  # not exercised in contract tests

    def _init_state(self) -> None:
        self._call_log.append("init_state")

    def _update_state(self, delta: np.ndarray) -> None:
        self._call_log.append("update_state")
        self._last_delta = delta

    def calculate_residual(self) -> np.ndarray:
        self._call_log.append("residual")
        return self._residuals.pop(0)

    def calculate_correction(self, residual: np.ndarray) -> np.ndarray:
        self._call_log.append("correction")
        return self._corrections.pop(0)

    # --- helpers to set up canned data before each test ---

    def _setup(self, residuals: list, corrections: list | None = None) -> None:
        self._residuals = [np.atleast_1d(np.array(r, dtype=float)) for r in residuals]
        self._corrections = [np.atleast_1d(np.array(c, dtype=float)) for c in (corrections or [])]
        self._call_log: list[str] = []
        self._last_delta = None


@pytest.fixture
def runner() -> FakeRunner:
    return FakeRunner(damping_factor=0.7, max_iter=10, max_residual=1e-3, warm_start=False)


# --- ABC enforcement ---

def test_missing_init_state_raises_on_instantiation():
    class Incomplete(SimulationRunner, NewtonRaphsonMixin):
        def run(self, wing_pool): ...
        def _update_state(self, delta): ...
        def calculate_residual(self): ...
        def calculate_correction(self, R): ...
        # _init_state not implemented

    with pytest.raises(TypeError):
        Incomplete(damping_factor=0.7, max_iter=10, max_residual=1e-3, warm_start=False)


def test_missing_calculate_residual_raises_on_instantiation():
    class Incomplete(SimulationRunner, NewtonRaphsonMixin):
        def run(self, wing_pool): ...
        def _init_state(self): ...
        def _update_state(self, delta): ...
        def calculate_correction(self, R): ...
        # calculate_residual not implemented

    with pytest.raises(TypeError):
        Incomplete(damping_factor=0.7, max_iter=10, max_residual=1e-3, warm_start=False)


def test_missing_calculate_correction_raises_on_instantiation():
    class Incomplete(SimulationRunner, NewtonRaphsonMixin):
        def run(self, wing_pool): ...
        def _init_state(self): ...
        def _update_state(self, delta): ...
        def calculate_residual(self): ...
        # calculate_correction not implemented

    with pytest.raises(TypeError):
        Incomplete(damping_factor=0.7, max_iter=10, max_residual=1e-3, warm_start=False)


def test_missing_update_state_raises_on_instantiation():
    class Incomplete(SimulationRunner, NewtonRaphsonMixin):
        def run(self, wing_pool): ...
        def _init_state(self): ...
        def calculate_residual(self): ...
        def calculate_correction(self, R): ...
        # _update_state not implemented

    with pytest.raises(TypeError):
        Incomplete(damping_factor=0.7, max_iter=10, max_residual=1e-3, warm_start=False)


# --- Convergence ---

def test_converges_when_residual_below_threshold(runner):
    runner._setup(residuals=[[0.0005]])
    R, converged = runner._run_newton_raphson()

    assert converged is True
    np.testing.assert_array_equal(R, [0.0005])


def test_converges_after_several_iterations(runner):
    runner._setup(
        residuals=[[1.0], [0.5], [0.0005]],
        corrections=[[0.1], [0.1]],
    )
    R, converged = runner._run_newton_raphson()

    assert converged is True
    np.testing.assert_array_equal(R, [0.0005])


# --- Divergence ---

def test_diverges_after_max_iter(runner):
    n = runner.max_iter
    runner._setup(
        residuals=[[1.0]] * (n + 1),
        corrections=[[0.1]] * n,
    )
    R, converged = runner._run_newton_raphson()

    assert converged is False
    np.testing.assert_array_equal(R, [1.0])


# --- Hook call order ---

def test_init_state_called_first(runner):
    runner._setup(residuals=[[0.0001]])
    runner._run_newton_raphson()

    assert runner._call_log[0] == "init_state"


def test_hook_call_order_one_iteration(runner):
    runner._setup(residuals=[[1.0], [0.0001]], corrections=[[0.1]])
    runner._run_newton_raphson()

    assert runner._call_log == ["init_state", "residual", "correction", "update_state", "residual"]


# --- Damping responsibility ---

def test_mixin_passes_raw_delta_to_update_state(runner):
    """The mixin must NOT apply damping — that is the runner's responsibility in _update_state."""
    correction = np.array([0.5])
    runner._setup(residuals=[[1.0], [0.0001]], corrections=[correction])
    runner._run_newton_raphson()

    np.testing.assert_array_equal(runner._last_delta, correction)
