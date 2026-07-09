# Design Spec: `config.yaml` Schema
*Date: 2026-06-28*

## Purpose

`config.yaml` is the single source of truth for a simulation run in the 1D QM Playground. It is read by the Python controller (`controller.py`), which validates it and passes it to each solver subprocess. The schema must cover the full pipeline (TISE, TDSE, Analysis) and remain agnostic of any particular solver implementation.

The canonical example file lives at `BSpline-CPP/1D-QM-Playground/config.yaml`.

---

## Top-Level Structure

```
run           — orchestration flags and output directory
physics       — physical constants (mass, hbar)
bspline       — B-spline basis parameters (shared by all solvers)
potential     — piecewise potential definition
tise          — TISE solver settings (includes continuum sub-block)
tdse          — TDSE solver settings
analysis      — quantities computed from TISE + TDSE output
visualization — what to plot (raw solver output and analysis output)
```

---

## Section Reference

### `run`

| Field | Type | Description |
|---|---|---|
| `run_tise` | bool | Whether to invoke the TISE solver |
| `run_tdse` | bool | Whether to invoke the TDSE solver |
| `run_analysis` | bool | Whether to invoke the Analysis script |
| `output_dir` | string | Root directory for all solver output; subdirs `tise/` and `tdse/` are created here |

---

### `physics`

| Field | Type | Description |
|---|---|---|
| `mass` | float | Particle mass. Set to `1.0` for atomic units. |
| `hbar` | float | Reduced Planck constant. Set to `1.0` for atomic units. |

---

### `bspline`

Defines the B-spline basis used to represent wavefunctions. Shared by all solvers.

| Field | Type | Description |
|---|---|---|
| `n_nodes` | int | Number of knot points on the spatial grid |
| `order` | int | B-spline order k (polynomial degree = k − 1; continuity order = k − 2) |
| `domain` | [float, float] | Spatial domain `[x_min, x_max]` |

---

### `potential`

A piecewise potential V(x). Defined as a YAML list; each element is a **string** encoding a dict with two keys.

#### Format

```
"{'domain': '<interval>', 'function': '<expression>'}"
```

Each string must be parseable by Python's `ast.literal_eval` after extracting the dict, or by an equivalent parser in C++ (yaml-cpp + custom interval/expression parser).

#### `domain` — interval notation

Specifies the subdomain on which this piece of the potential applies. Uses standard mathematical interval notation:

| Syntax | Meaning |
|---|---|
| `[a, b]` | a ≤ x ≤ b (closed on both ends) |
| `(a, b)` | a < x < b (open on both ends) |
| `[a, b)` | a ≤ x < b (closed left, open right) |
| `(a, b]` | a < x ≤ b (open left, closed right) |

`inf` is a valid bound (e.g., `[0, inf)` for the positive half-line).

The pieces in the list should together tile the spatial domain without gaps or overlaps. Behavior on gaps or overlaps is undefined.

#### `function` — mathematical expression

A mathematical expression in the spatial variable `x`. All constants must be **literal numeric values** — no references to other config fields (e.g., `L`, `hbar`, `mass`) are supported. The centrifugal barrier and all other terms must be written explicitly.

Supported operators and functions are determined by the expression parser used in each solver. At minimum, `+`, `-`, `*`, `/`, `^` (power), and parentheses must be supported.

#### Examples

```yaml
# Hydrogen-like, L=1 effective potential: V(x) = -1/x + L(L+1)/(2x^2) with L=1
potential:
  - "{'domain': '(0, 100]', 'function': '-1/x + 1/x^2'}"

# Particle in a box with a rectangular barrier
potential:
  - "{'domain': '[0, 5)',  'function': '0'}"
  - "{'domain': '[5, 6]',  'function': '10'}"
  - "{'domain': '(6, 10]', 'function': '0'}"

# Harmonic oscillator
potential:
  - "{'domain': '(-inf, inf)', 'function': '0.5 * 0.25 * x^2'}"
```

---

### `tise`

Settings for the TISE solver. The solver always computes all bound states below the ionization threshold — this count is output, not input.

| Field | Type | Description |
|---|---|---|
| `n_pts_eigenstate` | int | Number of spatial grid points for eigenstate wavefunction output |
| `error_threshold` | float | Eigenvalue accuracy cutoff for reporting (not a solver tolerance) |

#### `tise.continuum`

| Field | Type | Description |
|---|---|---|
| `enabled` | bool | Whether to compute continuum pseudostates and phase shifts |
| `E_threshold` | float | Lower bound of the continuum spectrum range |
| `E_max` | float | Upper bound of the continuum spectrum range |
| `n_energies` | int | Number of energy grid points in `[E_threshold, E_max]` |
| `n_pts` | int | Spatial grid points per continuum state output |

---

### `tdse`

Settings for the TDSE solver (eigenstate expansion).

#### `tdse.initial_state`

| Field | Type | Description |
|---|---|---|
| `type` | string | `eigenstate` or `gaussian` |
| `index` | int | *(eigenstate only)* 0-indexed bound state from TISE output |
| `position` | float | *(gaussian only)* Wavepacket center x₀ |
| `momentum` | float | *(gaussian only)* Mean momentum k₀ |
| `width` | float | *(gaussian only)* Gaussian width σ |

#### Other `tdse` fields

| Field | Type | Description |
|---|---|---|
| `dt` | float | Output time resolution — a snapshot is written at t = 0, dt, 2dt, … |
| `t_final` | float | Total propagation time |
| `snapshot_interval` | int | Number of dt steps between wavefunction snapshot writes |

---

### `analysis`

Quantities computed from TISE + TDSE output, per `docs/planning/architecture-07-02.md`. All fields require `run_tdse: true` — every quantity is a function of the propagated state $|\Psi(t)\rangle$.

| Field | Type | Description |
|---|---|---|
| `bound_state_populations` | bool | $p_n(t) = |\langle\phi_n|\Psi(t)\rangle|^2$ — population of each bound eigenstate over time |
| `asymptotic_populations` | bool | $p_n(\infty) = |\langle\phi_n|\Psi(t_\text{final})\rangle|^2$ — final-time population of each bound eigenstate |
| `asymptotic_distribution` | bool | $dP_\alpha/dE = |\langle\psi^-_{\alpha E}|\Psi(t_\text{final})\rangle|^2$ — energy-resolved population density in the continuum. Requires `tise.continuum.enabled: true`. |
| `expectation_values.x` | bool | $\langle\hat{x}\rangle(t)$ |
| `expectation_values.p` | bool | $\langle\hat{p}\rangle(t)$ |
| `expectation_values.T` | bool | $\langle\hat{T}\rangle(t)$ — kinetic energy |
| `expectation_values.V` | bool | $\langle\hat{V}\rangle(t)$ — potential energy |
| `expectation_values.H` | bool | $\langle\hat{H}\rangle(t)$ — total energy |
| `interval_probability.enabled` | bool | Whether to compute interval probability (optional quantity) |
| `interval_probability.intervals` | list of [float, float] | $P_{[x_a,x_b]}(t) = \int_{x_a}^{x_b} |\psi(x,t)|^2\, dx$ for each `[x_a, x_b]` pair listed |

---

### `visualization`

What to plot. Spans raw solver output (TISE, TDSE) and `analysis`-computed quantities alike — this is deliberately a separate block from `analysis` because `eigenstates`, `phase_shifts`, and `time_evolution` are not computed by Analysis; they are raw TISE/TDSE output being displayed. Boolean toggles only; plot parameters (ranges, state subsets) are not yet part of the schema.

| Field | Type | Description |
|---|---|---|
| `eigenstates` | bool | TISE: plot $|\phi_n(x)|^2$ for bound states |
| `phase_shifts` | bool | TISE: plot $\delta(E)$ and $d\delta/dE$. Requires `tise.continuum.enabled: true`. |
| `time_evolution` | bool | TDSE: plot/animate $|\Psi(x,t)|^2$. Requires `run_tdse: true`. |
| `bound_state_populations` | bool | Analysis: plot $p_n(t)$ |
| `asymptotic_populations` | bool | Analysis: plot $p_n(\infty)$ |
| `asymptotic_distribution` | bool | Analysis: plot $dP_\alpha/dE$. Requires `tise.continuum.enabled: true`. |
| `expectation_values` | bool | Analysis: plot $\langle\hat{O}\rangle(t)$ for each operator enabled under `analysis.expectation_values` |

---

## Design Decisions

| Decision | Choice | Rationale |
|---|---|---|
| Potential format | List of string-encoded dicts | Keeps math expressions opaque to YAML; avoids YAML escaping issues with operators |
| Domain notation | Standard interval notation `[a, b)` | Unambiguous, familiar to physicists |
| Function expressions | Literal numbers only (no variable refs) | Simpler parser; avoids coupling expression evaluation to config state |
| Centrifugal / angular momentum | Fully user-defined in `potential` | No special-casing for radial vs. Cartesian; general 1D problem space |
| `bspline` placement | Top-level | Shared infrastructure used by both TISE and TDSE; not solver-specific |
| `continuum` placement | Nested under `tise` | Continuum states are computed by the TISE solver; nesting expresses that dependency |
| `tise.n_states` | Removed | Bound state count is output, not input; solver always computes all sub-threshold states |
| `tdse.chebyshev_order` | Removed | TDSE uses eigenstate expansion (exact), not Chebyshev propagation |
| Split `analysis` → `analysis` + `visualization` | Two top-level blocks | `analysis` = quantities Analysis computes from TISE+TDSE output; `visualization` = what to plot (raw solver output or analysis output). Conflating them mislabeled TISE/TDSE raw outputs as "analysis". |
| `expectation_values` format | Per-operator boolean flags | Consistent with existing flag style (`run`, `tise.continuum`); explicit and self-documenting |
| `visualization` scope | Boolean toggles only | Plot parameters (ranges, subsets) deferred to analysis.py / later schema iteration |
