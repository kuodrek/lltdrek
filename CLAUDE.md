# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this project is

`lltdrek` is a Python library implementing a **nonlinear Lifting Line Theory (LLT)** solver for aerodynamic analysis of wing systems. It simulates lift, drag, and moment coefficients across angles of attack, supporting multi-surface configurations (e.g., wing + tail), ground effect, and angular velocity (damping derivatives).

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

## Architecture

### Core workflow

The typical usage pattern (see `example.py`):

1. **Load airfoil data** with `load_folder(folder)` — reads `.txt` (Cl polar tables) and `.dat` (geometry) files from a directory.
2. **Define wing geometry** with `Wing(...)` — then call `wing.generate_mesh()` and `wing.setup_airfoil_data(flight_condition, airfoils_data)`.
3. **Define flight conditions** with `FlightCondition(V_inf, nu, rho, angles_of_attack, h, ...)`.
4. **Create a `WingPool`** from a list of wings + flight condition. This automatically mirrors each wing (y-symmetry) and pre-computes all induced velocity matrices.
5. **Run simulation** with `Simulation(...).run(wing_pool)` → returns `list[SimulationResult]`.
6. **Post-process** with `PostProcessing.get_coefficients(wing_pool, results)` → returns `list[ProcessedSimulationResults]` containing global and per-surface `Coefficients`.

### Key classes

| Class | Location | Role |
|---|---|---|
| `Wing` | `models/wing.py` | Geometry + mesh. Stores collocation points, vertice points, panel normals (`u_n`), chord vectors (`u_a`), and span vectors (`u_s`). |
| `FlightCondition` | `models/flight_condition.py` | Atmospheric + flight state (V_inf, nu, rho, AOAs, angular rates, ground effect). |
| `WingPool` | `models/wingpool.py` | Assembles the full system. Builds mirrored wings, pre-computes `system_induced_velocities` and `system_freestream_velocities` for every AOA. |
| `Simulation` | `models/simulation.py` | Newton-Raphson iterator. Solves for dimensionless vortex strength `G` per panel. |
| `PostProcessing` | `models/post_processing.py` | Converts `G` solutions to aerodynamic force/moment coefficients. Also has `get_aerodynamic_center()` via binary search. |

### Simulation solver details (`simulation/main_equations.py`)

The nonlinear solve uses Newton-Raphson iteration:
- `calculate_main_equation` → residual vector `R(G)`
- `calculate_corrector_equation` → Jacobian solve `[J]ΔG = -R`
- `calculate_main_equation_simplified` → linearized version (used as an initial guess when `simulation_mode="linear_first"`)

`SimulationModes.LINEAR_FIRST` solves the linear equations first per alpha to warm-start the nonlinear solver. `SimulationModes.LATEST_SOLUTION` reuses the previous alpha's solution.

### WingPool internals

- Each `Wing` in `wing_list` is automatically paired with a `_mirrored` copy (negated y-coordinates). Mirrored wings are skipped in the main residual loop but contribute induced velocities.
- `system_induced_velocities[alpha][wing_i][wing_j]` is a 3D array indexed `[cp_i, cp_j, 3]` — the induced velocity at collocation point `i` due to panel `j` of `wing_j`.
- `G_dict` maps surface name → 1D numpy array of vortex strengths (one per panel).

### Airfoil data format (`.txt` files)

Each `.txt` file contains polar data for one airfoil across multiple Reynolds numbers. Structure per Reynolds block:
```
reynolds,<Re_value>
cm0,<value>
<aoa>,<cl>
<aoa>,<cl>
...
<blank line>
```
Linear coefficients (`cl_alpha`, `cl0`) are fitted via `np.polyfit` over a configurable AOA range (default 0–8°). Panels spanning two airfoils use linear interpolation weighted by spanwise position (`merge_parameter`).

### Geometry utilities (`utils/geometry.py`)

Panel distribution supports `"linear"` or `"cosine"` spanwise spacing. Euler angles (dihedral, twist, sweep) define the panel orientation matrices (`u_a`, `u_n`, `u_s`).
