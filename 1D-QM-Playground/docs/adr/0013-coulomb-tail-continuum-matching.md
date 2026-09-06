# ADR-0013: Coulomb-Tail Continuum Matching (Supersedes ADR-0010)

- **Status:** Accepted (implemented)
- **Date:** 2026-09-06

## Context

ADR-0010 deferred Coulomb-tail continuum matching: `tise::matchAsymptotic`
implemented only the flat-asymptote formula ($\sin/\cos$ matching against a
free-particle plane wave), which is not the correct asymptotic form for a
genuine $-Z/r$ Coulomb tail — the true exterior solution is a Coulomb wave
function $F_\ell(\eta,\rho)$/$G_\ell(\eta,\rho)$, whose phase grows
logarithmically with $r$, not linearly. Running continuum construction on a
hydrogen-like config with an unbounded domain silently produced a wrong
phase shift and normalization (confirmed before this work: ~0.02–1 rad
error depending on energy, vs. the true Coulomb-matched value).

This ADR reverses that deferral: the Coulomb-matching branch is implemented
and wired end-to-end.

## Decision

1. **New config field `tise.continuum.l`** (default 0). Angular momentum
   cannot be inferred from `classifyAsymptote`'s tail-shape fit — a
   centrifugal term $\ell(\ell+1)/2x^2$ decays faster than the $1/x$ Coulomb
   term and is asymptotically invisible to it — so it must be supplied
   explicitly, matching whatever centrifugal term is baked into the
   potential expression string.

2. **`AsymptoteClassification::fittedPowerLawCoefficient`** — a new field
   extracting the leading coefficient $C$ in $V\sim V_\infty + C/x^p$ (for
   Coulomb, $C=-Z$) from `classifyAsymptote`'s existing power-law tail fit.

3. **Hand-rolled Coulomb wave function evaluator**, not GSL. No
   special-function library is linked into this build
   (`TISE/CMakeLists.txt`: LAPACK/BLAS/muparser/nlohmann-json/yaml-cpp/GTest
   only); GSL's local availability was uncertain in the development sandbox
   (`pkg-config --exists gsl` failed; package-manager queries hung), and a
   hand-rolled evaluator matches this codebase's established pattern of
   self-contained numerics. The Coulomb phase shift
   $\sigma_\ell(\eta)=\arg\Gamma(\ell+1+i\eta)$ uses a Lanczos approximation
   for the complex Gamma function (validated to ~1e-15 against
   `mpmath.gamma`) plus the exact recurrence
   $\sigma_\ell=\sigma_0+\sum_{k=1}^{\ell}\arctan(\eta/k)$. $F_\ell$/$G_\ell$
   and their derivatives are obtained by starting from a WKB-amplitude-
   corrected asymptotic approximation far out and integrating the exact
   ODE $u''+(1-2\eta/\rho-\ell(\ell+1)/\rho^2)u=0$ inward to the requested
   $\rho$.

4. **Numerov's method, not Runge-Kutta, for the inward integration.** A
   first implementation used 4th-order Runge-Kutta. It was numerically
   correct but had an accuracy/cost tradeoff that failed on the project's
   real energy range (`E_max=10`, requiring $\rho$ into the hundreds and a
   `farMultiplier` of ~100 to keep the WKB starting error small): the
   accumulated phase error over that long oscillatory integration forced an
   impractically fine step (`stepSize=0.01`), and the real `config.yaml`
   (`E_max=10, n_energies=100`) failed to complete within a 60-second
   bound. Coarsening RK4's step directly made accuracy *worse*, not just
   cheaper — increasing `farMultiplier` at a fixed coarse step increased
   accumulated phase drift over the longer path faster than it reduced the
   WKB starting-point error. Numerov's method — the standard technique for
   a second-order ODE with no first-derivative term, exactly this ODE's
   form — has much better phase behavior per step for oscillatory
   equations. Switching to it allowed `farMultiplier=50, stepSize=0.1`
   (from RK4's `farMultiplier=100, stepSize=0.01`) at roughly 1/17th the
   integration steps.

5. **`matchAsymptotic` gains a Coulomb branch**, dispatched from
   `tise_solver_main.cpp` when `classifyAsymptote` reports
   `AsymptoteSubType::Coulomb`. The energy-independent coefficient
   `C = fittedPowerLawCoefficient` and `l = tise.continuum.l` are passed
   through as `std::optional<std::pair<int, Real>> coulombLC`; the
   energy-dependent Sommerfeld parameter $\eta=C/k$ ($k=\sqrt{2E}$) is
   recomputed fresh for every energy in the grid inside `matchAsymptotic`,
   not once per call — an early version computed it once outside the
   per-energy loop, which is wrong since $k$ varies across the grid; this
   was caught and fixed before shipping (see `TISE/tise.cpp`'s
   `amplitudeAndDelta` lambda). The matching itself mirrors the flat case's
   structure: match $\psi(R)$/$\psi'(R)$ against
   $A_E[\cos\delta\cdot F_\ell(\eta,kR)+\sin\delta\cdot G_\ell(\eta,kR)]$
   and its derivative, using the Wronskian $W=FG'-GF'$ (computed directly
   from the evaluated functions, never assumed to be exactly $\pm1$) in a
   2×2 linear solve.

## Consequences

- Hydrogen-like continuum configs (`tests/hydrogen.yaml`, `config.yaml`)
  now get physically correct Coulomb-matched phase shifts. Validated
  end-to-end (recovering a known synthetic phase shift through the full
  matching formula, cross-checked against `mpmath.coulombf`/`coulombg`)
  to ~6e-3 rad worst-case across $\ell=0..2$ and $E\in[0.05,10]$ — the full
  energy range this project's configs use — versus the flat formula's
  ~0.02–1 rad error on the same cases. This is not a special-function-
  library-grade arbitrary-precision result; `farMultiplier`/`stepSize`
  trade accuracy for speed if a different point on that curve is needed.
- `config.yaml` (`E_max=10, n_energies=100`) now completes in ~7 seconds
  (measured), down from failing to complete within 60 seconds under the
  initial RK4-based implementation.
- A potential's domain must be declared unbounded (`'(0, inf)'`, not
  capped exactly at the box edge like `'(0, 100]'`) for
  `classifyAsymptote`'s own `coveredBeyondDomain` precondition to detect a
  Coulomb tail at all — a config that declares the domain as ending
  exactly at `rMax` will never engage this feature, even if the true
  physics is Coulomb-tailed. This is existing, unchanged behavior of
  `classifyAsymptote`, just newly load-bearing for this feature.
- `tise.continuum.l` defaults to 0 (s-wave) if omitted, silently — a config
  with a centrifugal term baked into its potential expression but no
  matching `l` value will be matched against the wrong Coulomb functions
  with no diagnostic. This mirrors the same class of silent-mismatch risk
  ADR-0010 already flagged as unavoidable without an explicit `l` field.
- Full C++ (`ctest`) and Python (`pytest tests/`) regression suites pass
  with these changes; see
  `docs/planning/tise-coulomb-tail-and-boundary-validation-plan.md`'s
  implementation-status section for before/after numbers.

## Source

`docs/planning/coulomb-tail-continuum-matching.md` (derivation sketch and
integration plan, followed closely); `docs/tests/reports/9d0f04b/verify_known_solutions.py`
(the pre-existing `mpmath.coulombf`-based Python reference this
implementation was validated against); ADR-0010 (the deferral this
supersedes); `TISE/tise.hpp`/`TISE/tise.cpp` (module-level comment above
`evaluateCoulombFunctions`'s declaration has the full validation record and
numerical-method history).
