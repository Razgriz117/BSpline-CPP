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
- **Do you need only bound states, or also continuum/scattering states?** Bound states ($E<0$-style, confined) always come out of a solve. Continuum (scattering) states are optional and cost extra config ([§3.6](#36-tisecontinuum)).
- **What does $V(x)$ do beyond your box?** This determines whether you should declare your potential's domain as bounded (matching the box) or genuinely unbounded past it:
  - **Flat** (goes to 0, or a constant): declare it bounded, matching the box — this is the common case (e.g. a finite square well capped by hard walls).
  - **Coulomb-like** ($V \sim C/x$, unbounded): declare the piece's domain as **genuinely unbounded** (e.g. `(0, inf)`), even though your box is finite — this is what lets the solver detect the true tail shape and match continuum states against Coulomb wave functions instead of the wrong (flat-asymptote) formula. See [§3.6.1](#361-the-l-field-and-coulomb-tail-matching), which also quotes the confirming console/`warnings.json` message you'll see once the match engages.
  - **Irregular** (some other unbounded power law, e.g. $V\sim x^{-1.5}$): also declare it unbounded; the solver fits the tail's power-law exponent, and if it matches neither flat ($p\approx0$) nor Coulomb ($p\approx1$), classifies the tail "Case 3"/Irregular and warns rather than silently mismatching it — see [§3.6.2](#362-irregular-tails-case-3-and-tapering) for the actual message and what the solver does about it.

**What "warns" actually means.** Both cases above can make `tise_solver` print a *warning* instead of failing outright: it appends a `{"category": "physics", "message": "..."}` entry to `<output_dir>/tise/warnings.json` and prints an identical `tise_solver: warning: ...` line to the console the instant it fires; `controller.py` then reprints every entry from that file a second time once the solve finishes ([§5](#5-read-your-output)). A warning means the result was still computed and written to disk, but should be treated with suspicion — usually because your grid is too coarse, your box is too small for the energies you asked for, or your physical setup is ambiguous in some specific, named way. This is categorically different from a **hard error** (a malformed config, a bad `function` expression, a domain that doesn't tile — the gotchas in [§3.2](#32-physics)/[§3.4](#34-potential)): errors abort the pipeline *before* `tise_solver` produces any output at all, so nothing under `<output_dir>/` gets written. Every warning message the solver can emit, with what it means and how to react, is catalogued in [§7](#7-troubleshooting-reference).

## 2. Find your starting point

`tests/*.yaml` are real, runnable, physics-verified reference configs (checked against closed-form analytic solutions — see [§9](#9-further-reading)). Find the row closest to your problem and start from that file rather than from scratch. **The "Demonstrates" column names the specific physics mechanism or solver edge case each file is the reference implementation for** — what it's checked against, and which part of the solver it exercises — so you can jump straight to a config that's structurally like your problem.

| Your potential looks like... | Start from | Demonstrates |
|---|---|---|
| Nothing special, just want to see the pipeline work | [`tests/free_particle.yaml`](../../tests/free_particle.yaml) | **Baseline: the trivial flat-tail case.** $V=0$ inside a hard-walled box (an infinite square well, not a literal free particle — the walls are still hard); checked against $E_n=n^2\pi^2/2L^2$, $\phi_n=\sqrt{2/L}\sin(n\pi x/L)$. With continuum enabled, the exact phase shift is $\delta(E)\equiv0$ at every energy — the simplest possible correctness check for flat-tail continuum matching ([§3.6](#36-tisecontinuum)), and the file that `finite_square_well.yaml` (below) diffs against for a real discontinuity. |
| Same, but with a non-default `physics.mass`/`hbar` | [`tests/free_particle_general_units.yaml`](../../tests/free_particle_general_units.yaml) | **Tests: `physics.mass`/`physics.hbar` at non-default values** ([§3.2](#32-physics)). Identical $V=0$ box to `free_particle.yaml`, but `mass: 2.0, hbar: 0.5` — the first reference config to set a non-default `physics:` block. Checked against the general closed forms $E_n=\hbar^2n^2\pi^2/(2\,\text{mass}\,L^2)$ and $k=\sqrt{2\,\text{mass}\,E}/\hbar$; `E_max` is rescaled by $\hbar^2/\text{mass}$ relative to `free_particle.yaml` to keep the same physical $k$-range. $\delta\equiv0$ still holds exactly — a flat potential scatters nothing regardless of units. |
| A well or barrier with a genuine jump in $V$ | [`tests/finite_square_well.yaml`](../../tests/finite_square_well.yaml) | **Tests: a real step discontinuity → nonzero, energy-dependent phase shifts.** Two pieces — $V=-1$ on $[0,10)$, $V=0$ on $[10,100]$ — otherwise identical `bspline:`/`tise:` blocks to `free_particle.yaml` so [§3.4](#34-potential) can diff the two directly. Checked against the transcendental bound-state condition $K\cot(Ka)=-\kappa$ (4 bound states) and the closed-form phase shift $\delta(E)=\arctan(\tfrac{k}{K}\tan Ka)-ka$ — unlike `free_particle`'s trivial $\delta\equiv0$, this is genuinely energy-dependent. |
| A smooth, confining potential (e.g. $x^2$) | [`tests/harmonic_oscillator.yaml`](../../tests/harmonic_oscillator.yaml) | **Tests: pure grid-resolution accuracy, no confounds.** A single smooth piece, no continuum (bound-only); box walls at $\pm20$ stay well outside every relevant turning point through $n\approx200$, so there's no box-edge, discontinuity, or singularity to muddy the error — checked against $E_n=n+\tfrac12$ (Hermite-function eigenstates). This is the cleanest file for reasoning about `n_nodes`/`order` accuracy, and is the source of [§3.3](#33-bspline)'s own "~3.3 nodes per de Broglie wavelength for $10^{-6}$ accuracy" rule. |
| A Coulomb-tailed / hydrogenic potential | [`tests/hydrogen.yaml`](../../tests/hydrogen.yaml) | **Tests: Coulomb-tail continuum matching, and the `tise.continuum.l` field.** Single piece on $(0,\infty)$, $V=-1/x+1/x^2$ — Coulomb plus an $\ell(\ell+1)/x^2$ centrifugal term baked into the expression for $\ell=1$, matched by `l: 1`; the only one of these files that uses `l`, and a mismatch between it and the baked-in centrifugal term is a silent-wrong-answer trap ([§3.6.1](#361-the-l-field-and-coulomb-tail-matching)). Checked against $E_n=-1/2n^2$ ($n\ge2$; $\ell=1$ excludes 1s); because the tail is declared genuinely unbounded, continuum states are matched against real Coulomb functions $F_1,G_1$ instead of plane waves, so — like `free_particle` — the target is $\delta\equiv0$, but this time testing the Coulomb-matching machinery specifically. |
| A genuine singularity **inside** your domain | [`tests/interior_singularity.yaml`](../../tests/interior_singularity.yaml) | **Tests: excluding an interior singularity, and why continuum can't be matched across a split domain.** Two half-open pieces around $x=20$ — the "double-open-interval idiom" — split $[0,40]$ into two physically decoupled regions (a field-free box, and a repulsive-Coulomb-in-a-box); automatic strategic node placement at the join handles the rest. Bound states on both sides check out against box-state energies / zeros of $F_0$, but continuum construction is refused — for a structurally different reason than the row below: here the domain itself splits into two regions that share no physical continuum, not a singularity sitting exactly at the matching edge. |
| A singularity **exactly at** your box wall | [`tests/right_edge_singularity.yaml`](../../tests/right_edge_singularity.yaml) | **Tests: a singularity exactly at the box wall → continuum refused outright.** Single piece with a repulsive-Coulomb-like singularity placed exactly at `bspline.domain`'s right edge, where continuum phase-shift matching would need to evaluate. Bound states check out against zeros of $F_0(1/k,100k)$ and are completely unaffected; continuum refuses outright instead of degrading — no `phase_shifts.dat`/`continuum_state_*.dat` at all ([§7](#7-troubleshooting-reference)). |
| An unbounded tail that's neither flat nor Coulomb | [`tests/case3_irregular_tail.yaml`](../../tests/case3_irregular_tail.yaml) | **Tests: the Irregular/"Case 3" tail classification and taper warning** ([§3.6.2](#362-irregular-tails-case-3-and-tapering)). Domain genuinely unbounded but starting at `bspline.domain: [0.1, 50.0]` — just short of the true origin singularity in $V=1/x^{1.5}$ — a different idiom from `interior_singularity.yaml`'s double-open-interval (this one keeps the box away from an edge singularity rather than excluding an interior point). The tail's fitted exponent ($p\approx1.5$) is neither flat nor Coulomb; continuum is disabled here, so the taper is correctly skipped for the bound-state solve. No closed form — checked against Richardson-extrapolated finite-difference eigenvalues. |

**Note on plotting these three.** `interior_singularity.yaml` and `right_edge_singularity.yaml` set `visualization: {eigenstates: true, bound_states_squared: true}` — every computed state in both is a genuine, verified bound state of a fully confining system (the interior/edge singularity forces $\psi\to0$ there exactly like a hard wall, per the two rows above), so $|\phi_n(x)|^2$ is physically meaningful. `case3_irregular_tail.yaml` sets `eigenstates: true` only, deliberately *without* `bound_states_squared`: its potential is purely repulsive and decays to $0$ at infinity, so — unlike the other two, where the box edge is a real wall — it has no true bound states at all ($V\ge0$ everywhere forbids $E<0$); every computed state there is a box-discretized stand-in for what would be a continuum state if its irregular tail had a closed form to match against. Squaring those into "probability density" plots would misrepresent them as confined, so they stay raw — matching `free_particle.yaml`'s same choice for the same reason. None of the three ever produces `phase_shifts.png`/`continuum_NNN.png`: continuum construction is refused for the first two and disabled by config for the third, so that data never exists to plot.

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

Both fields default to `1.0` (atomic units) if the block or either field is omitted — this remains the simplest option for the common case, and every worked example elsewhere in this guide leaves it out. If you *do* need a different particle mass or a rescaled $\hbar$, any **positive** value now works end to end: the kinetic-energy matrix element, the accuracy-ceiling calculation (`E_acc`, [§3.3](#33-bspline)), and every continuum formula in [§3.6](#36-tisecontinuum) (the wavenumber $k=\sqrt{2\,\text{mass}\,E}/\hbar$, the Sommerfeld parameter, and the continuum-state normalization) all honor whatever you set here — see [ADR-0017](../adr/0017-mass-hbar-generalization.md) for the derivation.

The only requirement is physical sanity — zero or negative values are rejected:

```
physics.mass must be positive (got -1.000000); zero or negative mass is not physically meaningful.
```

(`hbar` gets the identical treatment.) One practical gotcha worth knowing before you reach for a non-default mass: a potential's `function` expression can't reference `physics.mass` (per [§3.4](#34-potential)'s "no cross-field references" rule), so any coefficient baked into your potential string (e.g. a harmonic oscillator's spring constant, or a Coulomb tail's baked-in centrifugal term) keeps its literal value regardless of what you set here — only the *kinetic* side of the equation, and the continuum matching built on top of it, actually scales with `mass`/`hbar`. Changing `physics.mass` alone on an existing config changes the physical problem you're solving (e.g. a differently-shaped potential well relative to the kinetic term), not just a label — plan your expected closed-form comparison accordingly, the way `tests/free_particle_general_units.yaml` does for the flat-potential case (§2).

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

**If `n_nodes`/`order` are unreasonable** (e.g. `n_nodes < 2` or `order < 1`), the failure surfaces as `tise_solver: BSpline::init failed: n_nodes must be >= 2 (got ...)` or `order must be >= 1 (got ...)`, naming the bad field and its actual value (the numeric `BSpline::init` code, e.g. `-1`/`-2`, is kept alongside it for anyone tracing the underlying return contract).

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

`[`/`]` = inclusive bound, `(`/`)` = exclusive bound, in any combination: `[a,b]`, `(a,b)`, `[a,b)`, `(a,b]`. `inf`/`infinity` (case-insensitive) is a valid bound, e.g. `(0, inf)`. Numeric bounds accept scientific notation (e.g. `1e5`, `-1.5E-3`), matching what Python's `float()` accepts on the config-validation side.

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
- `E_threshold`/`E_max` — the energy range to sample. `E_max` must be strictly greater than `E_threshold`, or the solver refuses to run: `tise.continuum.E_max must be greater than E_threshold (got E_threshold=..., E_max=...)`.
- `n_energies` — must be a positive integer, or the solver refuses to run: `tise.continuum.n_energies must be a positive integer`.
- `n_pts` — spatial resolution of each `continuum_state_NNN.dat` output file.

**If your `E_max` exceeds what your grid can resolve**, you'll see a warning (not a hard failure) in `warnings.json`:
```
Warning: requested continuum E_max=50 exceeds the basis accuracy ceiling E_acc=0.049348 (set by the
B-spline node spacing); results at energies above E_acc are unreliable.
```
Fix by raising `bspline.n_nodes` or lowering `E_max` — see [§3.3](#33-bspline)'s resolution discussion (and remember the practical ceiling is closer to `E_acc/2.7` than the raw reported value).

**A flat right-edge tail also needs to be small *at* the wall, not just far past it.** `classifyAsymptote` only checks the tail's functional *shape* well beyond the box edge — not whether $V$ has actually decayed to zero exactly at `x=rMax`. If it hasn't (the box is too small for the continuum energies you asked for), you'll see:
```
Warning: potential at the right domain edge x=100 is V(rMax)=-0.02, not negligible compared to the
smallest requested continuum energy E=0.1 (ratio 0.2 exceeds 0.01); matchAsymptotic's flat-asymptote
matching assumes V(rMax) is approximately zero, so the resulting phase shifts may be inaccurate --
consider enlarging bspline.domain.
```
(Illustrative numbers, same convention as the warning above — none of the 7 reference configs are misconfigured enough to trigger this one.) Only fires for a **flat** right asymptote with continuum enabled, and checks against the *smallest* requested energy, not `E_max` (the relative distortion from a fixed leftover $V(r_\text{max})$ is worst at the smallest energy). A genuine Coulomb tail is exempt — its nonzero $V(r_\text{max})$ is expected and handled by its own matching path instead ([§3.6.1](#361-the-l-field-and-coulomb-tail-matching)). Fix by enlarging `bspline.domain`.

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

#### 3.6.2 Irregular tails ("Case 3") and tapering

If a piece's declared domain is genuinely unbounded ([§1](#1-characterize-your-physics-problem)) and the solver's tail fit comes out as neither flat ($p\approx0$) nor Coulomb ($p\approx1$), it's classified **Irregular / "Case 3"** and produces two chained `warnings.json` entries — regardless of whether `tise.continuum.enabled` is set (this is a tail-shape classification, not a continuum-only concern).

First, the classifier's own detection message, e.g. from `tests/case3_irregular_tail.yaml` (`1/x^1.5` on a genuinely unbounded domain):
```
Warning: potential asymptote on the right side is irregular (fitted power-law exponent p=1.5);
the potential will be smoothly tapered to zero over a transition width delta=4.99 approaching
the box boundary (see docs/planning/boundary-condition-case-3-smoothing.md), avoiding an abrupt
truncation. This remains an approximation -- the true asymptotic tail is not analytically known
-- so continuum normalization will be approximate.
```
Second, a summary that actually decides what happens to your bound-state solve, gated on `continuum.enabled`:
- **Continuum disabled** (this is what `case3_irregular_tail.yaml` ships with): `potential's right-edge tail is Irregular (Case 3); continuum is disabled, so the raw (untapered) potential is integrated to the wall -- the bound-state spectrum is unaffected.`
- **Continuum enabled**: `potential's right-edge tail is Irregular (Case 3); tapering it to 0 over a transition width of <delta> before x=<rMax> rather than integrating the raw tail up to the wall (see evaluateWindowedPotential).`

**Read the second message, not just the first.** The first message always says the potential "will be smoothly tapered" — that's the classifier's own recommendation, made before it even knows whether continuum is enabled. Whether a taper is actually applied to your solve is decided by the second message: if continuum is disabled, no taper touches your bound-state solve and the first message's wording is stale in that context, not a bug in your config.

### 3.7 `tdse:` and `analysis:` (specified, not runnable)

`config.yaml`'s `tdse:` (initial state, gauge, driving field, time step) and `analysis:` (populations, expectation values) blocks are fully specified in the schema — you'll see them in `config.yaml`'s own shipped example — but **nothing consumes them yet**. There is no `tdse_solver` binary. Setting `run.run_tdse: true` aborts the pipeline immediately ([§3.1](#31-run)). Leave these blocks out, or leave them as shipped with `run.run_tdse: false` — either way they have zero effect on a TISE-only run.

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

`controller.py` validates your config ([§3](#3-build-a-config-field-by-field)'s tiling/type checks), runs `tise_solver` into `<output_dir>/tise/`, prints any `warnings.json` entries to the console, then — if `run_analysis: true` — runs `analysis.py` against that output to produce plots.

**Validation happens entirely before any output file is written.** If your config is rejected (or the solver hits a hard error), nothing stale is left behind under `<output_dir>/`, so a failed run's output directory is always either absent or complete from the *last successful* run — never a half-written mess from the failed one.

## 5. Read your output

Under `<output_dir>/tise/`:

| File | Contents |
|---|---|
| `eigenvalues.dat` | index, $E_n$ for **every** computed state — not pre-filtered to bound states. Negative $E$ = bound; see the informational count in `warnings.json` (below) for how many. |
| `eigenvectors.dat` | the full B-spline coefficient matrix. |
| `eigenstate_NNN.dat` (+ `.png` if `visualization.eigenstates: true`) | $x$, $\phi_n(x)$ — one pair per computed state. |
| `hamiltonian.dat` / `overlap.dat` | the banded H/S matrices the eigenproblem was built from. |
| `phase_shifts.dat` (+ `.png` if `visualization.phase_shifts: true`) | $\varepsilon_i$, $\delta(\varepsilon_i)$, $d\delta/dE$ — only present if continuum is enabled *and* not refused ([§7](#7-troubleshooting-reference) — right-edge or interior singularities). |
| `continuum_state_NNN.dat` (+ `continuum_NNN.png`) | $x$, $\psi_{\varepsilon_i}(x)$ per continuum energy — same presence condition as `phase_shifts.dat`. |
| `warnings.json` | **always present.** An array of `{"category": "physics", "message": "..."}` entries — see [§1](#1-characterize-your-physics-problem) for what a warning means; always includes an informational bound-state count, plus anything from [§3.3](#33-bspline)/[§3.6](#36-tisecontinuum)/[§3.6.2](#362-irregular-tails-case-3-and-tapering)/[§7](#7-troubleshooting-reference). |

## 6. Iterating on a config you already have

Quick task → section pointers, once you've already got a working file and want to change something:

- *"A bound state is colliding with the wall"* → enlarge `bspline.domain` ([§3.3](#33-bspline)).
- *"My `E_max` exceeds the accuracy ceiling"* → raise `n_nodes` or lower `E_max` ([§3.3](#33-bspline)/[§3.6](#36-tisecontinuum)).
- *"I need to add a piece to my potential"* → [§3.4](#34-potential) (tiling rules + the `free_particle`→`finite_square_well` diff).
- *"I want to turn on continuum for an existing bound-only config"* → [§3.6](#36-tisecontinuum), and check whether your tail is flat or Coulomb ([§3.6.1](#361-the-l-field-and-coulomb-tail-matching)) before you do.
- *"An 'Irregular'/'Case 3' warning showed up"* → [§3.6.2](#362-irregular-tails-case-3-and-tapering).
- *"A 'V(rMax) not negligible' warning showed up"* → [§3.6](#36-tisecontinuum).
- *"`tise_solver` crashed instead of printing an error"* → almost certainly a malformed `function` expression ([§3.4](#34-potential)).
- *"A message I don't recognize showed up"* → [§7](#7-troubleshooting-reference).

## 7. Troubleshooting reference

| Message (or its shape) | What it means | Fix / see |
|---|---|---|
| `potential does not cover domain lower/upper bound ...` | A potential piece doesn't reach the box edge. | [§3.4](#34-potential) tiling |
| `gap or overlap between potential pieces ...: ... does not meet ...` | Two adjacent pieces leave a gap. | [§3.4](#34-potential) tiling |
| `gap or overlap between potential pieces ...: both are inclusive at shared boundary ...` | Two adjacent pieces overlap at a shared point (both inclusive). | [§3.4](#34-potential) tiling — use one inclusive, one exclusive |
| `tise_solver: Function domain does not cover x = <x>` | A gap the Python-side check didn't catch. | [§3.4](#34-potential) |
| `overlapping potential piece domains: '<a>' and '<b>' both cover x=<x>` | The C++ solver's own independent overlap re-check caught something. | [§3.4](#34-potential) tiling |
| `tise.continuum.n_energies must be a positive integer` | `n_energies` is missing, zero, or negative. | [§3.6](#36-tisecontinuum) |
| `tise.continuum.E_max must be greater than E_threshold (got ...)` | `E_max` is equal to or less than `E_threshold`. | [§3.6](#36-tisecontinuum) |
| `physics.mass/physics.hbar must be positive (got ...)` | `physics:` block present with a zero or negative value. | [§3.2](#32-physics) — use a positive value |
| `terminate called after throwing an instance of 'mu::ParserError'` / `Aborted` | Malformed `function` expression — not caught cleanly. | [§3.4](#34-potential) — test the expression alone first |
| `BSpline::init failed: n_nodes must be >= 2 ...` / `order must be >= 1 ...` | `n_nodes < 2` / `order < 1`. | [§3.3](#33-bspline) |
| `Warning: requested continuum E_max=... exceeds the basis accuracy ceiling E_acc=...` | Grid too coarse for the requested continuum energy range. | [§3.3](#33-bspline)/[§3.6](#36-tisecontinuum) — raise `n_nodes` or lower `E_max` |
| `potential at the right domain edge x=... is V(rMax)=..., not negligible compared to ...` | Box isn't large enough for the smallest requested continuum energy — flat-asymptote matching assumes V≈0 at the wall. | [§3.6](#36-tisecontinuum) — enlarge `bspline.domain` |
| `Warning: continuum energy grid point E=... is within ... of confined eigenvalue ...` | A box-discretization artifact, not physical. | [§3.6](#36-tisecontinuum) — distrust that row, or shift the energy grid |
| `potential's right-edge tail is Irregular (Case 3); ...` (preceded by `Warning: potential asymptote on the ... side is irregular ...`) | Tail is unbounded but fits neither flat nor Coulomb; solver tapers it near the wall if continuum is enabled, otherwise leaves it untapered (bound states unaffected either way). | [§1](#1-characterize-your-physics-problem)/[§3.6.2](#362-irregular-tails-case-3-and-tapering) |
| `potential is singular at the right domain edge ... skipping continuum construction entirely` | Your potential is singular exactly at the box wall. | [§2](#2-find-your-starting-point)'s `right_edge_singularity.yaml` row — bound states are still valid |
| `potential has a singular join strictly inside the domain ... skipping continuum construction entirely` | Your potential has a genuine interior singularity, splitting the domain into two decoupled regions. | [§2](#2-find-your-starting-point)'s `interior_singularity.yaml` row — bound states are still valid |
| `bound state <j> ... appears to be colliding with the outer wall ...` | Box too small for that state. | [§3.3](#33-bspline) — enlarge `bspline.domain` |
| `"<n> of <m> computed states are below E=0.0"` | Informational only — tells you how many computed states are physically bound. | [§5](#5-read-your-output) |

## 8. Known hard limits

Things you cannot configure your way around today — see the linked ADRs for the full reasoning:

- **Node/knot placement is fully automatic** (uniform grid + auto-detected structural knots at potential discontinuities/singularities). There's no config field to manually place nodes, choose a placement formula, or pick a density scheme — `n_nodes`/`order` are your only levers ([ADR-0002](../adr/0002-defer-wkb-collocation.md), [ADR-0015](../adr/0015-defer-user-supplied-node-placement-formula.md)).
- **Delta-function potentials can't be written in the `potential` DSL at all** — there's no distributional-function support, only closed-form expressions. The closest approximation is a narrow, tall rectangular barrier via `Step`-style pieces, which is not equivalent physics ([ADR-0014](../adr/0014-defer-delta-potential-join-detection.md)).
- **No resonance, ionization-rate, or complex-eigenvalue calculations.** The solver is real-valued only; a confining box plus discretized continuum is the only option for anything resonance-adjacent ([ADR-0011](../adr/0011-defer-cap-outgoing-wave-ecs-boundary-conditions.md)).

## 9. Further reading

- [`docs/SDD.md`](../SDD.md) §6.1 "Configuration Schema" and §6.4 "Data Validation Rules" — the formal, exhaustive field-by-field reference this guide is a task-oriented complement to.
- [`docs/diagrams.md`](../diagrams.md) — schema and potential-DSL diagrams.
- [`docs/superpowers/specs/2026-06-28-config-yaml-schema-design.md`](../superpowers/specs/2026-06-28-config-yaml-schema-design.md) — the original design rationale for why the schema is shaped the way it is.
- [`docs/adr/`](../adr/) — design decisions behind every gotcha named above.
- [`docs/tests/reports/f4e8359/`](../tests/reports/f4e8359/) — the physics-verification reports (all 7 `tests/*.yaml` configs, all PASS) this guide's worked examples and accuracy numbers are drawn from.
- [`TISE/README.md`](../../TISE/README.md) — build instructions and dependency list.
- [`README.md`](../../README.md) — pipeline overview and inter-component contract.
- [`docs/planning/tdse-task-breakdown.md`](../planning/tdse-task-breakdown.md) — TDSE implementation status, for anyone curious when [§3.7](#37-tdse-and-analysis-specified-not-runnable) stops being "not yet."
