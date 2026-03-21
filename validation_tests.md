# Test Architecture for `src/new/`

## Why not test every combination?

The combinatorial explosion problem: with 2 runners × 2 airfoil DB implementations × 2 wing implementations × 2 parameter sets = 16 validation tests. Every new implementation multiplies the matrix. This is unsustainable.

The key insight: **where does the physics actually live?** In this codebase, it's concentrated in three independent layers:

- **Geometry**: induced velocity matrix (wing mesh, Biot-Savart)
- **Aerodynamics**: airfoil polar lookup (Cl/Cd vs α, Re)
- **Solver**: residual `R(G)` and its Jacobian solve

These layers are **orthogonal**. The runner doesn't care how the airfoil DB is implemented — it just calls `get_cl(alpha, Re)`. The airfoil DB doesn't care about wing geometry. This orthogonality is the key to avoiding combinatorial tests.

---

## Strategy: Building-Block Validation

This is called **building-block validation** in NASA/aerospace V&V standards (NASA-STD-7009). Validate at progressively higher levels of integration, not at the full combination matrix. Each block is validated once; higher-level tests assume lower blocks are correct.

Three principles:

1. **Validate ONE reference pipeline against experimental data** — proves the physics model is correct end-to-end
2. **For alternative implementations: prove equivalence to reference, not to experiment** — transitivity gives the guarantee
3. **Intermediate checkpoints for physically critical quantities** — catch bugs in components without running the full pipeline

---

## Reference Pipeline (canonical combination)

| Parameter | Value |
|---|---|
| Runner | `NonlinearLoopsRunner` |
| Wing | `WingLLT` |
| Airfoil DB | reference `AirfoilDatabase` (ported from legacy) |
| `warm_start` | `False` |

All other combinations are tested for consistency against this reference, not against experimental data.

---

## Directory Structure

```
tests/
├── unit/
│   ├── geometry/
│   │   └── test_wing_llt.py         # mesh, cosine distribution, panel geometry
│   ├── aerodynamics/
│   │   └── test_airfoil_database.py # polar lookup, Re interpolation
│   └── solver/
│       └── test_runners.py          # contract tests: all runners satisfy same interface
├── integration/
│   └── test_simulation.py           # WingPool + Simulation wired together, flat plate check
└── validation/
    ├── data/                         # airfoil .dat and experimental .exp.dat files (moved from validation/)
    ├── conftest.py                   # shared Wing + FlightCondition fixtures
    ├── test_artigo1.py               # NACA paper 1 — reference pipeline vs experimental CL/CD
    └── test_artigo2.py               # NACA paper 2 — reference pipeline vs experimental CL/CD
```

Note: `unit/` mirrors `src/new/` structure. `integration/` and `validation/` are organized by scenario, not by source file.

---

## Test Tiers

### 1. Unit Checkpoints — physics in isolation

**`tests/unit/geometry/test_wing_llt.py`**
- Cosine panel distribution matches analytical formula
- Panel chord, span, area are geometrically correct
- Induced velocity from a single vortex filament matches Biot-Savart analytical result

**`tests/unit/aerodynamics/test_airfoil_database.py`**
- `get_cl(alpha=0, Re)` ≈ 0 for symmetric airfoils
- Cl increases linearly with alpha in attached flow region (slope ≈ 2π/rad)
- Interpolation between Re values is continuous

**`tests/unit/solver/test_runners.py`**
- Contract test parametrized over all runners: given the same `WingPool`, all runners return a `SimulationResult` with the required fields
- `LinearRunner` converges in one step (no iteration)
- `NonlinearLoopsRunner` with `warm_start=True` converges in fewer iterations than without

### 2. Integration — components wired together

**`tests/integration/test_simulation.py`**
- `Simulation(equations="nonlinear", implementation="loops")` runs end-to-end on a flat plate (analytical CL ≈ π·α)
- `warm_start=True` result agrees with `warm_start=False` result within tolerance

### 3. Validation — V&V against experimental data

**`tests/validation/test_artigo1.py`** and **`tests/validation/test_artigo2.py`**
- Use the reference pipeline only
- Compare CL(α) and CD(α) curves against experimental values from each paper
- Tolerance: within experimental scatter (e.g., ΔCL < 0.05)
- These are the **only two tests** that use experimental data

### 4. Consistency — for new implementations

When a new implementation is added (e.g., `NonlinearNumpyRunner`, `WingLLTv2`, `AirfoilDatabaseB`):
- Add **one consistency test**: assert new implementation agrees with reference on a canonical case (simple rectangular wing, 5 AOAs)
- No re-running experimental validation — transitivity provides the guarantee
- Lives alongside the relevant unit tests (e.g., new runner → `tests/unit/solver/test_runners.py`)

---

## Summary Table

| Test type | What it proves | Count |
|---|---|---|
| Unit checkpoints | Each component is physically correct in isolation | ~5–10 |
| Integration | Components wire together correctly | ~2–3 |
| Validation | Reference pipeline matches experimental data | 2 (one per paper) |
| Consistency | New implementations agree with reference | 1 per new implementation |

---

## Running Tests

```bash
pytest tests/unit/        # fast, no experimental data needed
pytest tests/integration/ # fast, analytical check only
pytest tests/validation/  # slow — run manually or on release CI
```
