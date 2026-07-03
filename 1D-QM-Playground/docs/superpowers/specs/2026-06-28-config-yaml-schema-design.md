# Design Spec: `config.yaml` Schema
*Date: 2026-06-28*

## Purpose

`config.yaml` is the single source of truth for a simulation run in the 1D QM Playground. It is read by the Python controller (`controller.py`), which validates it and passes it to each solver subprocess. The schema must cover the full pipeline (TISE, TDSE, Analysis) and remain agnostic of any particular solver implementation.

The canonical example file lives at `BSpline-CPP/1D-QM-Playground/config.yaml`.

---

## Top-Level Structure

```
run       — orchestration flags and output directory
physics   — physical constants (mass, hbar)
bspline   — B-spline basis parameters (shared by all solvers)
potential — piecewise potential definition
tise      — TISE solver settings (includes continuum sub-block)
tdse      — TDSE solver settings
analysis  — post-processing and plot selection (see TODO below)
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

> **TODO:** Concretize this section's config format once decided. The set of computable quantities is now resolved (see `docs/planning/architecture-06-20.md` — "Analysis: what to compute", stakeholder feedback 2026-07-03):
>
> 1. **Bound-state populations as a function of time** — |⟨ψₙ|ψ(t)⟩|² for each bound eigenstate n
> 2. **Spectral distribution across channels at end of simulation** — energy-resolved probability distribution at t = t_final
> 3. **Expectation values of key observables** — ⟨x̂⟩(t), ⟨p̂⟩(t), ⟨T̂⟩(t), ⟨V̂⟩(t), ⟨Ĥ⟩(t)
> 4. **Interval probability** (optional) — P_{[xa,xb]}(t) = ∫_{xa}^{xb} |ψ(x,t)|² dx for user-specified [xa, xb]
>
> The YAML schema for selecting and parameterizing these quantities is TBD.

Current stub:

| Field | Type | Description |
|---|---|---|
| `plots.eigenstates` | bool | Plot bound-state probability densities |
| `plots.phase_shifts` | bool | Plot δ(E) and dδ/dE (requires `tise.continuum.enabled: true`) |
| `plots.time_evolution` | bool | Plot or animate |ψ(x,t)|² (requires `run_tdse: true`) |

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
| Analysis section | Stub + TODO | Quantities resolved (2026-07-03); config format TBD |
