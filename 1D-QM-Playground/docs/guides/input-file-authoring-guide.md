# Input-File Authoring Guide

This guide walks you, step by step, through writing your own `config.yaml`-shaped input file for the 1D-QM-Playground solver — from "I have a 1D potential in mind" to a working config and an understanding of what the output means.

**Audience:** anyone who wants to run their own potential through the solver. No C++/Python source-reading required.

**What this guide is not:** a quantum-mechanics tutorial (it assumes you already know what potential you want to solve), and not a TDSE authoring guide — the `tdse:`/`analysis:` blocks are fully specified in the schema but **not runnable yet** (see [§3.7](#37-tdse-and-analysis-specified-not-runnable)).

**Before you start:** build the solver and install the Python dependencies once, per the top-level [`README.md`](../../README.md#quick-start-the-real-pipeline) Quick Start. Everything below assumes `TISE/build/tise_solver` exists and `python3 controller.py --config <your-file>` works.

## Contents

1. [Characterize your physics problem](#1-characterize-your-physics-problem)
2. [Find your starting point](#2-find-your-starting-point)
3. [Build a config, field by field](#3-build-a-config-field-by-field)
   - [3.1 `run:`](#31-run)
   - [3.2 `physics:`](#32-physics)
   - [3.3 `bspline:`](#33-bspline)
   - [3.4 `potential:`](#34-potential)
   - [3.5 `tise:`](#35-tise)
   - [3.6 `tise.continuum:`](#36-tisecontinuum)
   - [3.7 `tdse:`/`analysis:` (specified, not runnable)](#37-tdse-and-analysis-specified-not-runnable)
   - [3.8 `visualization:`](#38-visualization)
4. [Run it](#4-run-it)
5. [Read your output](#5-read-your-output)
6. [Iterating on a config you already have](#6-iterating-on-a-config-you-already-have)
7. [Troubleshooting reference](#7-troubleshooting-reference)
8. [Known hard limits](#8-known-hard-limits)
9. [Further reading](#9-further-reading)

---

## 1. Characterize your physics problem

Before opening a YAML file, answer these on paper:

- **What is $V(x)$?** Write it as one or more pieces, each a closed-form expression over an interval. A single smooth potential is one piece; a step, a well, a barrier, or anything with a kink is multiple pieces.
- **What box contains the states you care about?** The solver always works on a finite domain `[x_min, x_max]` (a hard Dirichlet wall at each end, unless the boundary-condition machinery detects and matches an unbounded tail — see below). Pick it big enough that the bound states you want aren't touching the wall.
- **Do you need only bound states, or also continuum/scattering states?** Bound states ($E<0$-style, confined) always come out of a solve. Continuum (scattering) states are optional and cost extra config (§3.6).
- **What does $V(x)$ do beyond your box?** This determines whether you should declare your potential's domain as bounded (matching the box) or genuinely unbounded past it:
  - **Flat** (goes to 0, or a constant): declare it bounded, matching the box — this is the common case (e.g. a finite square well capped by hard walls).
  - **Coulomb-like** ($V \sim C/x$, unbounded): declare the piece's domain as **genuinely unbounded** (e.g. `(0, inf)`), even though your box is finite — this is what lets the solver detect the true tail shape and match continuum states against Coulomb wave functions instead of the wrong (flat-asymptote) formula. See [§3.6.1](#361-the-l-field-and-coulomb-tail-matching).
  - **Irregular** (some other unbounded power law): also declare it unbounded; the solver detects this as "Case 3" and warns rather than silently mismatching it.

## 2. Find your starting point

`tests/*.yaml` are real, runnable, physics-verified reference configs (checked against closed-form analytic solutions — see [§9](#9-further-reading)). Find the row closest to your problem and start from that file rather than from scratch.

| Your potential looks like... | Start from | Demonstrates |
|---|---|---|
| Nothing special, just want to see the pipeline work | [`tests/free_particle.yaml`](../../tests/free_particle.yaml) | The simplest possible config: $V=0$ in a box (an infinite square well, not a literal free particle — the walls are still hard). |
| A well or barrier with a genuine jump in $V$ | [`tests/finite_square_well.yaml`](../../tests/finite_square_well.yaml) | A two-piece potential with a real step discontinuity; energy-dependent phase shifts. See §3.4's worked diff against `free_particle.yaml`. |
| A smooth, confining potential (e.g. $x^2$) | [`tests/harmonic_oscillator.yaml`](../../tests/harmonic_oscillator.yaml) | A single smooth piece, no continuum (bound-only); the cleanest example for reasoning about `n_nodes`/`order` accuracy — see §3.3. |
| A Coulomb-tailed / hydrogenic potential | [`tests/hydrogen.yaml`](../../tests/hydrogen.yaml) | The `tise.continuum.l` field and unbounded-domain declaration for Coulomb-tail continuum matching — see §3.6.1. |
| A genuine singularity **inside** your domain | [`tests/interior_singularity.yaml`](../../tests/interior_singularity.yaml) | The double-open-interval idiom for excluding a point singularity; automatic strategic node placement handles the rest. Bound states are unaffected, but continuum construction is refused (same as the right-edge-singular row below) — an interior singular join splits your domain into two physically decoupled regions, so a "continuum state" matched at `bspline.domain`'s right edge wouldn't mean anything anyway. |
| A singularity **exactly at** your box wall | [`tests/right_edge_singularity.yaml`](../../tests/right_edge_singularity.yaml) | Continuum construction is refused outright (no `phase_shifts.dat`) rather than producing wrong output; bound states are unaffected. |
| An unbounded tail that's neither flat nor Coulomb | [`tests/case3_irregular_tail.yaml`](../../tests/case3_irregular_tail.yaml) | The "Case 3" irregular-asymptote taper/warning path; domain starting just short of a true origin singularity. |

**Note:** `interior_singularity.yaml`, `right_edge_singularity.yaml`, and `case3_irregular_tail.yaml` all ship with `run_analysis: false` and no `visualization:` block — running them produces `.dat` files only, no plots. That's deliberate (they're numerical-correctness checks), not a bug in the config.

## 3. Build a config, field by field

### 3.1 `run:`

```yaml
run:
  run_tise:     true
  run_tdse:     false
  run_analysis: true
  output_dir:   "./data/my_run"
```

- `run_tise` (required) — run the TISE solver.
- `run_tdse` (optional) — **must stay `false`.** Setting it `true` aborts the entire pipeline before anything runs: `controller.py: run.run_tdse is not yet supported (Phase 5); set run.run_tdse: false`. See [§3.7](#37-tdse-and-analysis-specified-not-runnable).
- `run_analysis` (optional) — run `analysis.py` afterward to produce plots from the TISE output.
- `output_dir` (required) — everything gets written under `<output_dir>/tise/` (and, once TDSE ships, `<output_dir>/tdse/`).

### 3.2 `physics:`

```yaml
physics:
  mass: 1.0
  hbar: 1.0
```

**Recommendation: omit this block entirely.** Despite looking like free parameters, the solver hardcodes atomic units (`mass=1`, `hbar=1`) internally everywhere — the kinetic-energy matrix element, the accuracy-ceiling calculation, the continuum wavenumber $k=\sqrt{2E}$. If you include this block with anything other than `1.0`, the solver refuses to run:

```
physics.mass is fixed at 1.0 internally (fillBandedMatrices' kinetic-energy term/computeEAcc/k=sqrt(2E)
all hardcode it); setting it to anything else is not yet supported -- remove the field or set it to 1.0.
```

(`hbar` gets the identical treatment.) If you leave the block out entirely, nothing checks it — that's the simplest option.

### 3.3 `bspline:`

```yaml
bspline:
  n_nodes: 51
  order:   12
  domain:  [0.0, 100.0]
```

- `domain` — your box `[x_min, x_max]`, from [§1](#1-characterize-your-physics-problem).
- `n_nodes` — number of grid breakpoints. More nodes = finer resolution = higher accuracy, at the cost of a larger eigenproblem.
- `order` — B-spline order $k$ (degree $k-1$, continuity $C^{k-2}$ between ordinary knots). Total basis functions = `n_nodes + order - 2`.

**How many nodes do you actually need?** Node/knot placement is otherwise fully automatic (see [§8](#8-known-hard-limits)) — `n_nodes`/`order` are your only real levers on resolution. `tests/harmonic_oscillator.yaml`'s reference verification (the cleanest possible test — no discontinuities, no singularities, no box-wall contamination) derived a concrete rule of thumb: **$10^{-6}$ eigenvalue accuracy needs about 3.3 nodes per local de Broglie wavelength**,
$$h \;\lesssim\; \frac{2\pi}{3.3\sqrt{2\,(E - V_\text{min})}}$$
for the highest-energy state you care about. The solver's own built-in accuracy-ceiling warning ([§3.6](#36-tisecontinuum)) is looser than this — it only guarantees the 2-nodes/wavelength Nyquist floor — so treat the reported `E_acc` as roughly `E_acc/2.7`, not as your real ceiling. See [`docs/tests/reports/f4e8359/harmonic_oscillator.md`](../tests/reports/f4e8359/harmonic_oscillator.md) for the full derivation.

**If `n_nodes`/`order` are unreasonable** (e.g. `n_nodes < 2` or `order < 1`), the failure surfaces late and generically as `tise_solver: BSpline::init failed with code -1` (bad `n_nodes`) or `-2` (bad `order`) — not a friendly named-field message.

### 3.4 `potential:`

This is the field you'll spend the most time on. It's a YAML list of strings, each a `{'domain': '...', 'function': '...'}` piece.

**A minimal example, and the smallest possible step up from it.** `tests/free_particle.yaml` is one piece:
```yaml
potential:
  - "{'domain': '[0, 100]', 'function': '0'}"
```
`tests/finite_square_well.yaml` is otherwise identical — same `bspline:`, same `tise:` block — except this becomes two pieces, adding a real discontinuity in $V$ at $x=10$:
```yaml
potential:
  - "{'domain': '[0, 10)',  'function': '-1.0'}"
  - "{'domain': '[10, 100]', 'function': '0'}"
```
That's the whole lesson for adding a piece: split the domain, give each sub-interval its own expression, make sure the pieces meet exactly (see tiling rules below).

#### Domain interval syntax

`[`/`]` = inclusive bound, `(`/`)` = exclusive bound, in any combination: `[a,b]`, `(a,b)`, `[a,b)`, `(a,b]`. `inf`/`infinity` (case-insensitive) is a valid bound, e.g. `(0, inf)`.

> **Gotcha — write bounds as plain decimals, not scientific notation.** `python3 controller.py`'s config validation uses Python's `float()` for interval bounds, which happily accepts `"1e5"`. The C++ solver's interval parser does not — it only accepts plain decimal notation. A domain like `'[0, 1e5]'` **passes Python-side validation and then crashes the actual solve**: `tise_solver: Invalid interval: [0, 1e5]`. Write `100000` instead.

#### Tiling: pieces must exactly cover your `bspline.domain`

The pieces together must tile `bspline.domain` with no gaps and no overlaps. `controller.py` checks this before ever invoking the solver, with specific messages:

- A piece doesn't reach the domain's lower/upper edge: `potential does not cover domain lower bound 0.0: lowest piece '[5, 100]' starts at 5.0` (and the analogous message for the upper bound).
- Two adjacent pieces don't meet: `gap or overlap between potential pieces '[0, 10)' and '[12, 100]': 10.0 does not meet 12.0`.
- Two adjacent pieces are both inclusive at the same point: `gap or overlap between potential pieces ...: both are inclusive at shared boundary ...` — ambiguous (which piece's `function` applies exactly at that point?), so it's rejected.

**Exception, and the idiom for excluding a singular point:** a shared boundary where **both** sides are *exclusive* (e.g. `[0,20)` next to `(20,40]`) is tolerated as a measure-zero gap — this is exactly how you exclude a single point your potential is singular at. See [`tests/interior_singularity.yaml`](../../tests/interior_singularity.yaml), which uses `[0,20)` + `(20,40]` to exclude $x=20$, where its potential has a genuine `1/(x-20)` singularity.

**A piece's domain may (and for Coulomb-tail matching, must) extend past `bspline.domain`.** `tests/hydrogen.yaml` declares its single piece on `(0, inf)` even though the box is `[0.0, 100.0]` — see [§3.6.1](#361-the-l-field-and-coulomb-tail-matching) for why this matters.

If the C++ solver is run directly (bypassing `controller.py`'s pre-check) and a gap exists, it surfaces lazily, the first time a grid point lands in it: `tise_solver: Function domain does not cover x = 10.625`.

#### The `function` expression

Evaluated by [muparser](https://beltoforion.de/en/muparser/), with **`x` as the only variable**. No references to other config fields are allowed — all constants must be literal numbers baked directly into the string (e.g. write `0.5 * 1.0 * x^2`, not `0.5 * mass * x^2` — see `tests/harmonic_oscillator.yaml`).

Supported operators: `+ - * / ^`, comparisons `< <= > >= == !=`, logical `&& ||`, ternary `?:`. Supported functions: `sin cos tan asin acos atan atan2 sinh cosh tanh asinh acosh atanh log2 log10 ln exp abs sqrt rint sign sum avg min max`; constants `_pi`, `_e`.

> **Gotcha — a malformed expression does not fail cleanly.** A syntax error in `function` (unbalanced parens, a stray operator, etc.) is not caught anywhere in the solver and crashes the whole process: `terminate called after throwing an instance of 'mu::ParserError'` followed by `Aborted`, instead of a clean `tise_solver: ...` message. **Test a new or edited expression in a tiny, throwaway config first** before dropping it into a real run — if `tise_solver` aborts instead of printing a normal error, the `function` string is where to look.

### 3.5 `tise:`

```yaml
tise:
  n_pts_eigenstate: 301
  error_threshold:  1.0e-10
```

- `n_pts_eigenstate` — spatial resolution (grid points) of each `eigenstate_NNN.dat` output file.
- `error_threshold` — required by the schema, but **not actually consumed by `tise_solver`** today (only the separate, older `H-BoundStates` demo binary reads it). Include it for schema-completeness; it has no effect on a real run.

### 3.6 `tise.continuum:`

Only needed if you want scattering/continuum states, not just bound states.

```yaml
tise:
  continuum:
    enabled:     true
    E_threshold: 0.0
    E_max:       0.5
    n_energies:  5
    n_pts:       500
```

- `enabled` — if omitted or `false`, none of the rest of this block matters.
- `E_threshold`/`E_max` — the energy range to sample. **Nothing checks that `E_max > E_threshold`** — get this backwards and you'll get a degenerate or empty energy grid with no diagnostic.
- `n_energies` — must be a positive integer, or the solver refuses to run: `tise.continuum.n_energies must be a positive integer`.
- `n_pts` — spatial resolution of each `continuum_state_NNN.dat` output file.

**If your `E_max` exceeds what your grid can resolve**, you'll see a warning (not a hard failure) in `warnings.json`:
```
Warning: requested continuum E_max=50 exceeds the basis accuracy ceiling E_acc=0.049348 (set by the
B-spline node spacing); results at energies above E_acc are unreliable.
```
Fix by raising `bspline.n_nodes` or lowering `E_max` — see [§3.3](#33-bspline)'s resolution discussion (and remember the practical ceiling is closer to `E_acc/2.7` than the raw reported value).

**A specific continuum energy can also land suspiciously close to one of your bound eigenvalues** — a finite-box discretization artifact, not real physics:
```
Warning: continuum energy grid point E=0.1 is within 0.00327788 of confined eigenvalue E_13=0.0967221;
this energy's continuum state is likely a finite-box discretization artifact (see ADR-0007), not a
physical feature -- treat with suspicion.
```
Nothing needs fixing — just distrust that one energy's row, or nudge `E_threshold`/`E_max`/`n_energies` so grid points don't land there.

#### 3.6.1 The `l` field and Coulomb-tail matching

```yaml
tise:
  continuum:
    l: 1   # angular momentum for Coulomb-tail continuum matching
```

> **This is the single biggest silent-wrong-answer trap in the whole schema.**

`l` is **only ever consulted** when the solver's own asymptote classifier decides your potential's right-edge tail is a genuine, unbounded Coulomb tail ($V \sim C/x$) — for an ordinary flat-asymptote potential (like a square well capped at the box wall), `l` is read but never used.

Two things both have to be true for Coulomb-tail matching (and therefore `l`) to engage at all:

1. **Your potential piece's declared domain must be genuinely unbounded** (e.g. `(0, inf)`), not capped exactly at the box edge — even if the "true" physics has a Coulomb tail, declaring the domain as `(0, 100]` instead of `(0, inf)` means the solver never even looks beyond the wall to detect it.
2. The tail actually has to fit a $\sim 1/x$ power law when the solver samples far beyond the box.

`tests/hydrogen.yaml` is the worked example — an $\ell=1$ hydrogenic radial potential, `-1/x + 1/x^2`, declared on `(0, inf)`, with `l: 1` (because the centrifugal term $\ell(\ell+1)/2x^2 = 1\cdot2/2 = 1$ matches the `+1/x^2` in `function`). Running it produces exactly this diagnostic, confirming the match engaged:
```
tise_solver: warning: potential's right-edge tail is Coulomb (V ~ -0.999996/x); continuum phase shifts
are matched against Coulomb wave functions F_l/G_l (l=1, from tise.continuum.l) instead of sin/cos --
see docs/planning/coulomb-tail-continuum-matching.md.
```

**`l` cannot be inferred automatically** — a centrifugal term decays faster than the $1/x$ Coulomb term and is asymptotically invisible to the tail-shape fit. If you omit `l`, it **silently defaults to `0`** (s-wave). If your `function` has a nonzero centrifugal term baked in but you forget to set a matching `l` (or set the wrong one), **you get no error and no warning** — just continuum results matched against the wrong Coulomb wave functions.

### 3.7 `tdse:` and `analysis:` (specified, not runnable)

`config.yaml`'s `tdse:` (initial state, gauge, driving field, time step) and `analysis:` (populations, expectation values) blocks are fully specified in the schema — you'll see them in `config.yaml`'s own shipped example — but **nothing consumes them yet**. There is no `tdse_solver` binary. Setting `run.run_tdse: true` aborts the pipeline immediately (§3.1). Leave these blocks out, or leave them as shipped with `run.run_tdse: false` — either way they have zero effect on a TISE-only run.

### 3.8 `visualization:`

```yaml
visualization:
  eigenstates:  true
  phase_shifts: true
```

Only three fields actually do anything today: `eigenstates` and `bound_states_squared` (gate `eigenstate_NNN.png` generation), and `phase_shifts` (gates `phase_shifts.png`, only meaningful if continuum is enabled). Everything else in this block (`time_evolution`, `bound_state_populations`, `asymptotic_populations`, `asymptotic_distribution`, `expectation_values`) is schema-only, awaiting TDSE — setting them `true` produces no plot and no error.

## 4. Run it

```bash
python3 controller.py --config <your-file>.yaml
```

`controller.py` validates your config (§3's tiling/type checks), runs `tise_solver` into `<output_dir>/tise/`, prints any `warnings.json` entries to the console, then — if `run_analysis: true` — runs `analysis.py` against that output to produce plots.

**Validation happens entirely before any output file is written.** If your config is rejected (or the solver hits a hard error), nothing stale is left behind under `<output_dir>/`, so a failed run's output directory is always either absent or complete from the *last successful* run — never a half-written mess from the failed one.

## 5. Read your output

Under `<output_dir>/tise/`:

| File | Contents |
|---|---|
| `eigenvalues.dat` | index, $E_n$ for **every** computed state — not pre-filtered to bound states. Negative $E$ = bound; see the informational count in `warnings.json` (below) for how many. |
| `eigenvectors.dat` | the full B-spline coefficient matrix. |
| `eigenstate_NNN.dat` (+ `.png` if `visualization.eigenstates: true`) | $x$, $\phi_n(x)$ — one pair per computed state. |
| `hamiltonian.dat` / `overlap.dat` | the banded H/S matrices the eigenproblem was built from. |
| `phase_shifts.dat` (+ `.png` if `visualization.phase_shifts: true`) | $\varepsilon_i$, $\delta(\varepsilon_i)$, $d\delta/dE$ — only present if continuum is enabled *and* not refused (§3.6.1, right-edge singularities). |
| `continuum_state_NNN.dat` (+ `continuum_NNN.png`) | $x$, $\psi_{\varepsilon_i}(x)$ per continuum energy — same presence condition as `phase_shifts.dat`. |
| `warnings.json` | **always present.** An array of `{"category": "physics", "message": "..."}` entries — always includes an informational bound-state count, plus anything from §3.3/§3.6/§7. |

## 6. Iterating on a config you already have

Quick task → section pointers, once you've already got a working file and want to change something:

- *"A bound state is colliding with the wall"* → enlarge `bspline.domain` ([§3.3](#33-bspline)).
- *"My `E_max` exceeds the accuracy ceiling"* → raise `n_nodes` or lower `E_max` ([§3.3](#33-bspline)/[§3.6](#36-tisecontinuum)).
- *"I need to add a piece to my potential"* → [§3.4](#34-potential) (tiling rules + the `free_particle`→`finite_square_well` diff).
- *"I want to turn on continuum for an existing bound-only config"* → [§3.6](#36-tisecontinuum), and check whether your tail is flat or Coulomb ([§3.6.1](#361-the-l-field-and-coulomb-tail-matching)) before you do.
- *"`tise_solver` crashed instead of printing an error"* → almost certainly a malformed `function` expression ([§3.4](#34-potential)).
- *"A message I don't recognize showed up"* → [§7](#7-troubleshooting-reference).

## 7. Troubleshooting reference

| Message (or its shape) | What it means | Fix / see |
|---|---|---|
| `potential does not cover domain lower/upper bound ...` | A potential piece doesn't reach the box edge. | [§3.4](#34-potential) tiling |
| `gap or overlap between potential pieces ...: ... does not meet ...` | Two adjacent pieces leave a gap. | [§3.4](#34-potential) tiling |
| `gap or overlap between potential pieces ...: both are inclusive at shared boundary ...` | Two adjacent pieces overlap at a shared point (both inclusive). | [§3.4](#34-potential) tiling — use one inclusive, one exclusive |
| `tise_solver: Function domain does not cover x = <x>` | A gap the Python-side check didn't catch (e.g. it slipped past because the C++ interval parser rejected a bound Python accepted — see the scientific-notation gotcha below). | [§3.4](#34-potential) |
| `tise_solver: Invalid interval: [0, 1e5]` | Domain bound used scientific notation; Python's validator accepts it, the C++ parser doesn't. | [§3.4](#34-potential) — write plain decimals |
| `overlapping potential piece domains: '<a>' and '<b>' both cover x=<x>` | The C++ solver's own independent overlap re-check caught something. | [§3.4](#34-potential) tiling |
| `tise.continuum.n_energies must be a positive integer` | `n_energies` is missing, zero, or negative. | [§3.6](#36-tisecontinuum) |
| `physics.mass/hbar is fixed at 1.0 internally ... remove the field or set it to 1.0` | `physics:` block present with a value other than `1.0`. | [§3.2](#32-physics) — omit the block |
| `terminate called after throwing an instance of 'mu::ParserError'` / `Aborted` | Malformed `function` expression — not caught cleanly. | [§3.4](#34-potential) — test the expression alone first |
| `BSpline::init failed with code -1` / `-2` | `n_nodes < 2` / `order < 1`. | [§3.3](#33-bspline) |
| `Warning: requested continuum E_max=... exceeds the basis accuracy ceiling E_acc=...` | Grid too coarse for the requested continuum energy range. | [§3.3](#33-bspline)/[§3.6](#36-tisecontinuum) — raise `n_nodes` or lower `E_max` |
| `Warning: continuum energy grid point E=... is within ... of confined eigenvalue ...` | A box-discretization artifact, not physical. | [§3.6](#36-tisecontinuum) — distrust that row, or shift the energy grid |
| `potential is singular at the right domain edge ... skipping continuum construction entirely` | Your potential is singular exactly at the box wall. | [§2](#2-find-your-starting-point)'s `right_edge_singularity.yaml` row — bound states are still valid |
| `potential has a singular join strictly inside the domain ... skipping continuum construction entirely` | Your potential has a genuine interior singularity, splitting the domain into two decoupled regions. | [§2](#2-find-your-starting-point)'s `interior_singularity.yaml` row — bound states are still valid |
| `bound state <j> ... appears to be colliding with the outer wall ...` | Box too small for that state. | [§3.3](#33-bspline) — enlarge `bspline.domain` |
| `"<n> of <m> computed states are below E=0.0"` | Informational only — tells you how many computed states are physically bound. | [§5](#5-read-your-output) |

## 8. Known hard limits

Things you cannot configure your way around today — see the linked ADRs for the full reasoning:

- **Node/knot placement is fully automatic** (uniform grid + auto-detected structural knots at potential discontinuities/singularities). There's no config field to manually place nodes, choose a placement formula, or pick a density scheme — `n_nodes`/`order` are your only levers ([ADR-0002](../adr/0002-defer-wkb-collocation.md), [ADR-0015](../adr/0015-defer-user-supplied-node-placement-formula.md)).
- **Delta-function potentials can't be written in the `potential` DSL at all** — there's no distributional-function support, only closed-form expressions. The closest approximation is a narrow, tall rectangular barrier via `Step`-style pieces, which is not equivalent physics ([ADR-0014](../adr/0014-defer-delta-potential-join-detection.md)).
- **No resonance, ionization-rate, or complex-eigenvalue calculations.** The solver is real-valued only; a confining box plus discretized continuum is the only option for anything resonance-adjacent ([ADR-0011](../adr/0011-defer-cap-outgoing-wave-ecs-boundary-conditions.md)).
- **`physics.mass`/`hbar` aren't really adjustable** — see [§3.2](#32-physics).

## 9. Further reading

- [`docs/SDD.md`](../SDD.md) §6.1 "Configuration Schema" and §6.4 "Data Validation Rules" — the formal, exhaustive field-by-field reference this guide is a task-oriented complement to.
- [`docs/diagrams.md`](../diagrams.md) — schema and potential-DSL diagrams.
- [`docs/superpowers/specs/2026-06-28-config-yaml-schema-design.md`](../superpowers/specs/2026-06-28-config-yaml-schema-design.md) — the original design rationale for why the schema is shaped the way it is.
- [`docs/adr/`](../adr/) — design decisions behind every gotcha named above.
- [`docs/tests/reports/f4e8359/`](../tests/reports/f4e8359/) — the physics-verification reports (all 7 `tests/*.yaml` configs, all PASS) this guide's worked examples and accuracy numbers are drawn from.
- [`TISE/README.md`](../../TISE/README.md) — build instructions and dependency list.
- [`README.md`](../../README.md) — pipeline overview and inter-component contract.
- [`docs/planning/tdse-task-breakdown.md`](../planning/tdse-task-breakdown.md) — TDSE implementation status, for anyone curious when §3.7 stops being "not yet."
