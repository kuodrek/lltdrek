# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

`lltdrek` is a Python library implementing a **nonlinear Lifting Line Theory (LLT)** solver for aerodynamic analysis of wing systems. It simulates lift, drag, and moment coefficients across angles of attack, supporting multi-surface configurations (e.g., wing + tail), ground effect, and angular velocity (damping derivatives).

## Current branch: `feature/system-redesign`

The codebase is being redesigned. **All new work happens in `src/new/`.** The legacy code in
`src/lltdrek/` is the reference for porting — do not extend it. After migration is complete,
`src/lltdrek/` will be removed and this file updated.

## Commands

**Install for development:**
```bash
pip install -e ".[dev]"
# or just install dependencies:
pip install -r requirements.txt
```

**Run tests:**
```bash
pytest
# single test file:
pytest tests/integration/validation/test_validation_1.py
```

**Lint / format (via pre-commit):**
```bash
pre-commit run --all-files
```
- Formatter: `black` with `--line-length=120`
- Import sorter: `isort`
- Linter: `flake8` with `--max-line-length=120 --ignore=E731`

---

## Redesign: `src/new/`

### Goals (from `anotacoes.md`)

- Wing and the object that assembles wings are separate concerns
- `AirfoilDatabase` is a standalone object — `Wing` does not own airfoil data
- `WingPool` owns airfoil data attachment and must validate database coverage before running
- `Simulation` responsibilities are split across focused `SimulationRunner` subclasses (Strategy Pattern)
- Velocity functions and data-loading functions still need to be refactored (pending)

### Design principles

- **ABCs for all interfaces.** Concrete implementations are named with a suffix: `WingLLT`, `NonlinearLoopsRunner`, etc.
- **Orthogonal parameters.** Each parameter on a class controls one independent axis of behavior. Do not combine axes into a single flat enum.
- **Strategy Pattern for swappable implementations.** The user-facing class (e.g., `Simulation`) is a factory that selects and holds the right strategy.
- **`WingPool` owns the assembly lifecycle.** Wings are passive geometry objects; they do not know about flight conditions or airfoil databases until `WingPool` attaches that data.

### Current structure

```
src/new/
├── aerodynamics/
│   └── airfoil_database.py     # AirfoilDatabase ABC
├── geometry/
│   ├── wing.py                 # Wing ABC (lifecycle + properties)
│   └── wing_llt.py             # WingLLT — concrete LLT implementation
├── system/
│   └── wing_pool.py            # WingPool — assembles wings + flight condition + airfoil DB
└── solver/
    ├── simulation.py           # Simulation — user-facing factory (equations, implementation, warm_start)
    ├── simulation_runner.py    # SimulationRunner ABC + SimulationResult dataclass
    ├── linear_runner.py        # LinearRunner — one-shot linear solve (stub)
    ├── nonlinear_loops_runner.py  # NonlinearLoopsRunner — Newton-Raphson, for-loops (stub)
    └── nonlinear_numpy_runner.py  # NonlinearNumpyRunner — Newton-Raphson, vectorized (stub)
```

### `WingPool` — `src/new/system/wing_pool.py`

Takes `wings: List[Wing]`, `flight_condition: FlightCondition`, `airfoil_db: AirfoilDatabase`.
In `__init__`, calls `wing.generate_mesh()` and `wing._setup_airfoil_data(flight_condition, airfoil_db)` for each wing.

### `Simulation` — `src/new/solver/simulation.py`

Three orthogonal params select the runner internally:

| `equations` | `implementation` | Runner selected |
|---|---|---|
| `"linear"` | (ignored) | `LinearRunner` |
| `"nonlinear"` | `"loops"` | `NonlinearLoopsRunner` |
| `"nonlinear"` | `"numpy"` | `NonlinearNumpyRunner` |

`warm_start: bool` — nonlinear only; uses linear solution as initial G per alpha.

### Still pending

- `AirfoilDatabase` concrete implementation (port from `load_folder` + polar lookup in legacy)
- `FlightCondition` in `src/new/aerodynamics/`
- Velocity utility functions (port from `src/lltdrek/utils/`)
- `NonlinearLoopsRunner.run()` implementation (port from `src/lltdrek/models/simulation.py`)
- `NonlinearNumpyRunner.run()` implementation (new, vectorized)
- `PostProcessing` equivalent
- Wing mirroring logic in `WingPool`
- Airfoil database validation in `WingPool`

---

## Legacy reference: `src/lltdrek/`

Read-only. Used as the source of truth for porting logic.

| Class | Location | Role |
|---|---|---|
| `Wing` | `models/wing.py` | Monolithic geometry + airfoil data |
| `FlightCondition` | `models/flight_condition.py` | V_inf, nu, rho, AOAs, angular rates, ground effect |
| `WingPool` | `models/wingpool.py` | Mirrors wings, pre-computes induced/freestream velocity matrices |
| `Simulation` | `models/simulation.py` | Newton-Raphson loop (to be split into runners) |
| `PostProcessing` | `models/post_processing.py` | G → force/moment coefficients |

Key equations in `src/lltdrek/simulation/main_equations.py`:
- `calculate_main_equation` (line 47) → residual R(G)
- `calculate_corrector_equation` (line 95) → Jacobian solve ΔG
- `calculate_main_equation_simplified` (line 9) → linearized solve (warm-start or LinearRunner)
