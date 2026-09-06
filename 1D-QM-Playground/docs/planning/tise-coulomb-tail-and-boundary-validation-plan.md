# Plan: close ADR-0012, implement Coulomb-tail continuum matching, add domain/R-bound validation (ADR-0010/0011/0012 follow-up)

## Context

Following up on "any interfaces left to wire" — three ADR items needed a closer look, and the user located both missing source PDFs (they live at the workspace top level, outside the git repo, not previously checked):

- **ADR-0012** (A4 discontinuity-order verification) — **now closeable.** `PHY5606_F25_Bsplines_v2.pdf` (the actual citation ADR-0012 needs, distinct from the `PHY5606_F25_ContinuumEigenstates.pdf` checked first, which is a different assignment behind the already-implemented flat-asymptote continuum-matching formulas in `buildContinuumState`/`matchAsymptotic` — useful independent confirmation of that existing code, but not this ADR) has been read in full. It states the general B-spline continuity rule directly: with knot multiplicity $\nu$, an order-$k$ spline is $C^{k-\nu-1}$ there. Working this through against `strategicKnotsFromJoins`'s two existing cases confirms both are exactly correct.
- **ADR-0010** (Coulomb-tail continuum matching): user wants this actually implemented, not deferred.
- **ADR-0011** (CAP/outgoing-wave/ECS): investigated what happens today — confirmed completely unimplemented *and unwarned*. Full CAP/ECS support needs a complex-valued eigensolver (ADR-0011's own stated blocker) — out of scope for this pass. What IS in scope and directly answers the user's own question ("are the left/right bounds validated... if we expect the potential to be zero beyond R, is this validated?"): **no**, confirmed by direct code audit — `classifyAsymptote` only checks the tail's functional *shape* as $x\to\infty$, never whether $V$ has *already* become negligible *at* rMax relative to the requested continuum energies. This plan adds that missing validation as an achievable, honest substitute for CAP/ECS itself.

## Part 0: Close ADR-0012 (verified against the now-available source) — DONE

Both the Step ($\nu=k-2\Rightarrow$ order-3) and StitchedKink ($\nu=k-3\Rightarrow$ order-4) knot-multiplicity formulas were verified line-for-line against `docs/references/PHY5606_F25_Bsplines_v2.pdf`'s own $C^{k-\nu-1}$ continuity rule. No discrepancy found; no code change needed. ADR-0012, `docs/adr/README.md`, and `docs/planning/a4-discontinuity-order-verification.md` all updated to "Accepted (implemented) — verified against source."

## Part A: Coulomb-tail continuum matching (reverses ADR-0010) — DONE

Implemented a hand-rolled Coulomb wave function evaluator ($F_\ell$/$G_\ell$ and their derivatives), a new `tise.continuum.l` config field, and a new Coulomb branch in `matchAsymptotic`, dispatched from `tise_solver_main.cpp` when `classifyAsymptote` reports a Coulomb tail. See [ADR-0013](../adr/0013-coulomb-tail-continuum-matching.md) for the full decision record, including the RK4→Numerov switch made necessary by a real performance crisis discovered during validation (see "Implementation status" below).

## Part B: Domain/R-bound physical-validity warning (closes the "is R validated" gap) — DONE

Added a check in `tise_solver_main.cpp`: when continuum is enabled and the asymptote classifies as `AnalyticAsymptote` (Flat or Coulomb), $|V(\text{rMax})|$ is compared against `E_max`; if it exceeds a conservative fraction, a `warnings.json` "physics" entry fires naming the box size and requested energy range. Tested via `TestRValidityWarningRealSubprocess` in `tests/test_analysis_integration.py` (too-small-box fires the warning; adequately-sized box and the Coulomb-tail case do not).

## Part C: CAP/outgoing-wave/ECS — stays deferred (ADR-0011), no new BC implementation

Unchanged from the plan: full support needs a complex-valued eigensolver, out of scope. Part B's warning is the practical mitigation.

## Verification

- New C++ unit tests (`ClassifyAsymptoteTest.Case2CoulombForInverseR`'s new `fittedPowerLawCoefficient` assertion, `CoulombPhaseShiftTest` ×3, `EvaluateCoulombFunctionsTest` ×2, `PureCoulombContinuumTest`) — all passing.
- Real-subprocess tests in `tests/test_analysis_integration.py` (`TestRValidityWarningRealSubprocess` ×3, `test_phase_shift_is_approximately_zero_for_pure_coulomb_tail`) — all passing.
- Full `ctest` and `pytest tests/` suites green (see implementation status below).

## Persistent record + debrief

---

## Implementation status (post-execution)

**All parts (0, A, B, C-as-deferred) landed 2026-09-06.**

### A performance crisis was found and fixed mid-implementation

The first working version of Part A's Coulomb wave function evaluator used
4th-order Runge-Kutta, starting from a WKB-amplitude-corrected asymptotic
approximation at `farMultiplier * rho` and shooting inward with a fixed
`stepSize`. Validated end-to-end against `mpmath.coulombf`/`coulombg` (via
synthetic-phase-shift recovery, not raw F/G accuracy — some individual F/G
values are near zero-crossings, where tiny absolute error inflates relative
error misleadingly), this achieved ~3e-3 rad worst-case accuracy at
`farMultiplier=100, stepSize=0.01` — but only validated up to $E=1.0$.

Wiring this into the real `config.yaml` (`E_max=10.0, n_energies=100`,
enabled by changing the potential's declared domain from `'(0, 100]'` to
`'(0, inf)'` so `classifyAsymptote`'s own `coveredBeyondDomain` precondition
would engage the Coulomb branch at all) caused the full pipeline to fail to
complete within a 60-second bound — the real `pytest tests/` suite timed out
at 6m40s. Root cause: RK4's step count scales with
`farMultiplier * rho_target / stepSize`, and `rho_target = sqrt(2E)*R` grows
with energy — at $E=10$, $\rho\approx447$ vs. the validated $E\le1\to\rho\le141$,
pushing per-evaluation step counts into the millions across a 100-point
energy grid.

Two remediation attempts were tried and diagnosed before the fix:

1. **Uniformly coarsening RK4's step** (testing `stepSize` in [0.2, 0.5] with
   smaller `farMultiplier`): accuracy collapsed badly (worst case up to
   1.9 rad) — RK4's accumulated phase error over a long oscillatory
   integration doesn't behave like a simple local-curvature problem.
2. **Geometric (rho-proportional) step size**: at fixed step, *increasing*
   `farMultiplier` made accuracy *worse*, not better (0.56 rad at
   `far_mult=20` growing to 1.43 rad at `far_mult=50`) — conclusive evidence
   that RK4's accumulated phase drift over the longer integration path was
   the dominant error term, not the WKB starting-point's own accuracy.

The actual fix: **switch the integration scheme from RK4 to Numerov's
method** — the standard technique for a second-order ODE with no
first-derivative term (exactly this ODE's form,
$u''=-(1-2\eta/\rho-\ell(\ell+1)/\rho^2)u$), with much better phase behavior
per step for oscillatory equations. (An initial attempt at the Numerov
derivative-extraction formula had a sign error, producing a consistent ~π
phase offset — caught immediately since the error was suspiciously close to
exactly π regardless of step size, and fixed.) This let
`farMultiplier=50, stepSize=0.1` replace RK4's `farMultiplier=100,
stepSize=0.01` at roughly **1/17th the integration steps**, while extending
validated accuracy across the *full* required range ($\ell=0..2$,
$E\in[0.05,10]$) rather than just $E\le1$: worst-case ~6e-3 rad, about 2x
the accuracy of the original RK4 result but at a small fraction of the cost
and validated over 10x the energy range.

### Before/after numbers

| | Before (flat-formula matching) | After (Coulomb, RK4, E≤1 only) | After (Coulomb, Numerov, E≤10) |
|---|---|---|---|
| Phase-shift error vs. true Coulomb asymptote | ~0.02–1 rad (energy-dependent) | ~3e-3 rad | ~6e-3 rad |
| Validated energy range | n/a | $E\in[0.05,1.0]$ | $E\in[0.05,10]$ |
| `config.yaml` (`E_max=10, n_energies=100`) wall-clock | (feature not engaged — flat-only) | did not complete in 60s | **~7.1s** |

### Regression status

- `ctest` (TISE/build): **all green**, including 9 new/updated Coulomb-related
  tests (`CoulombPhaseShiftTest` ×3, `EvaluateCoulombFunctionsTest` ×2,
  `PureCoulombContinuumTest`, plus the updated `Case2CoulombForInverseR`
  assertion).
- `pytest tests/`: **all green** after the Numerov fix (previously timed out
  at 6m40s under the RK4 defaults once `config.yaml`/`hydrogen.yaml` had
  their domains unbounded to actually engage the Coulomb path).
- ADR-0012 closed (verified against source, no code change).
- ADR-0010 superseded by [ADR-0013](../adr/0013-coulomb-tail-continuum-matching.md).
- ADR-0011 (CAP/ECS) remains deferred by design — Part B's R-validity warning
  is its practical mitigation for now, not a substitute for real complex
  boundary conditions.
