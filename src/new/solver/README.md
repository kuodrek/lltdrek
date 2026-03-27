# Solver Architecture

## 1. Overview

`Simulation` is the user-facing entry point. It acts as a factory that selects the right `SimulationRunner` based on three orthogonal parameters and delegates all computation to it.

```
┌─────────────────────────────────────────────────────────────┐
│  User code                                                  │
│                                                             │
│  sim = Simulation(equations="nonlinear", implementation=...)│
│  results = sim.run(wing_pool)                               │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│  Simulation (factory)                                       │
│                                                             │
│  _build_runner() selects the right SimulationRunner:        │
│                                                             │
│  equations="linear"       → LinearRunner                    │
│  equations="nonlinear"                                      │
│    implementation="loops" → NonlinearLoopsRunner            │
│    implementation="numpy" → NonlinearNumpyRunner            │
└──────────────────────┬──────────────────────────────────────┘
                       │
                       ▼
┌─────────────────────────────────────────────────────────────┐
│  SimulationRunner.run(wing_pool) → list[SimulationResult]   │
└─────────────────────────────────────────────────────────────┘
```

## 2. Class Hierarchy

Runner identity and solver behavior are composed separately:

```
SimulationRunner (ABC)            NewtonRaphsonMixin (ABC)
│ - damping_factor                │ - _run_newton_raphson()  ← concrete
│ - max_iter                      │
│ - max_residual                  │  @abstractmethod:
│ - warm_start                    │  - _init_state()
│                                 │  - _update_state(delta)
│ @abstractmethod:                │  - calculate_residual()
│ - run(wing_pool)                │  - calculate_correction(R)
│                                 │
├── LinearRunner                  │
│     one-shot linear solve       │
│     (no Newton iteration)       │
│                                 │
├── NonlinearLoopsRunner ─────────┘
│     for-loop math (validated baseline)
│
└── NonlinearNumpyRunner ─────────┘
      vectorized math
```

- **`SimulationRunner`** provides identity (what the runner *is*) and shared config.
- **`NewtonRaphsonMixin`** provides behavior (what the runner *does*) — a domain-agnostic iteration loop.
- Concrete runners inherit both and implement the four abstract hooks.

## 3. Newton-Raphson Iteration Flow

```
_run_newton_raphson()
│
├── _init_state()                    ← build self.state (IterationState)
│
├── R = calculate_residual()
│
└── loop (max_iter times):
    │
    ├── if max(|R|) < tolerance:
    │     return (R, converged=True)
    │
    ├── δ = calculate_correction(R)  ← solve Jacobian system
    │
    ├── _update_state(δ)             ← apply damped update, rebuild self.state
    │
    └── R = calculate_residual()
```

The mixin knows only `self.max_iter`, `self.max_residual`, and the four hooks.
All domain concepts (`WingPool`, `alpha`, `G_dict`, velocities) live in the concrete runner.

## 4. IterationState

Each nonlinear runner defines an `IterationState` dataclass to bundle per-iteration variables.
This makes three categories of state visually distinct when reading the code:

| Access pattern       | Category        | Changes when?              |
|----------------------|-----------------|----------------------------|
| `self.state.X`       | Iteration state | Every Newton iteration     |
| `self._wing_pool`    | Domain context  | Once per `run()` call      |
| `self.damping_factor`| Runner config   | Never (set at construction)|

## 5. How to Implement a New Runner

### Case A — New runner using Newton-Raphson (e.g. GPU-accelerated)

```python
from dataclasses import dataclass
import numpy as np
from new.solver.newton_raphson import NewtonRaphsonMixin
from new.solver.simulation_runner import SimulationResult, SimulationRunner
from new.system.wing_pool import WingPool


@dataclass
class IterationState:
    G: np.ndarray
    G_dict: dict[str, np.ndarray]
    total_velocity_dict: dict[str, np.ndarray]
    aoa_eff_dict: dict[str, np.ndarray]


class NonlinearGPURunner(SimulationRunner, NewtonRaphsonMixin):

    def run(self, wing_pool: WingPool) -> list[SimulationResult]:
        self._wing_pool = wing_pool
        results = []
        for alpha in wing_pool.flight_condition.alphas:
            self._alpha = alpha
            R, converged = self._run_newton_raphson()
            results.append(SimulationResult(alpha=alpha, G=self.state.G_dict, residual=R, converged=converged))
        return results

    def _init_state(self) -> None:
        # build self.state
        ...

    def _update_state(self, delta: np.ndarray) -> None:
        # apply damped correction, rebuild self.state
        ...

    def calculate_residual(self) -> np.ndarray:
        # your math here, reading self.state
        ...

    def calculate_correction(self, residual: np.ndarray) -> np.ndarray:
        # your math here, reading self.state
        ...
```

Then register the new runner in `Simulation._build_runner()` in [simulation.py](simulation.py).

### Case B — New solver method (e.g. Broyden's method)

1. Create a new mixin in `src/new/solver/`:

```python
from abc import ABC, abstractmethod
import numpy as np


class BroydenMixin(ABC):

    @abstractmethod
    def _init_state(self) -> None: ...

    @abstractmethod
    def _update_state(self, delta: np.ndarray) -> None: ...

    @abstractmethod
    def calculate_residual(self) -> np.ndarray: ...

    @abstractmethod
    def _init_jacobian_approx(self) -> None: ...

    @abstractmethod
    def _update_jacobian_approx(self, delta: np.ndarray, delta_R: np.ndarray) -> None: ...

    def _run_broyden(self) -> tuple[np.ndarray, bool]:
        self._init_state()
        self._init_jacobian_approx()
        R = self.calculate_residual()
        for _ in range(self.max_iter):
            if np.max(np.abs(R)) < self.max_residual:
                return R, True
            delta = ...  # solve using approximate inverse Jacobian
            self._update_state(delta)
            R_new = self.calculate_residual()
            self._update_jacobian_approx(delta, R_new - R)
            R = R_new
        return R, False
```

2. Create a runner that uses it:

```python
class NonlinearBroydenRunner(SimulationRunner, BroydenMixin):
    def run(self, wing_pool):
        ...
```

3. Register in `Simulation._build_runner()`.
