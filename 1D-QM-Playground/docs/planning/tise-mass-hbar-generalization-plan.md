# `physics.mass`/`physics.hbar` generalization

> Companion doc to `tise-release-readiness-plan.md`/`tise-known-solution-followup-plan.md`, same convention: plan first, "Implementation status (post-execution)" appended once the work actually lands.

## Context

`docs/planning/tise-release-readiness-plan.md`'s gap **G12** ("`physics.mass`/`physics.hbar` documented, never consumed, no error if set") was closed by Part C as a **guard-rail, not full generalization**: `tise_solver_main.cpp` was made to reject any `physics.mass`/`physics.hbar` value other than exactly `1.0`, turning a previously silent wrong-physics bug into an honest config error. The formulas themselves — the kinetic-energy matrix element, `computeEAcc`, and the continuum wavenumber/Sommerfeld parameter/energy-normalization amplitude inside `matchAsymptotic` — all still hardcoded atomic units (`mass=1`, `hbar=1`) internally.

This plan finishes that deferred half: `mass`/`hbar` become real, honored parameters throughout the TISE solver, entirely within the TISE module (confirmed no dependency on the not-yet-implemented TDSE module — `TISE/time_evolution.cpp`'s `tevol::` namespace already takes `mass`/`hbar` as explicit parameters, and the TDSE task-breakdown doc commits new TDSE code to the same convention; this work brings TISE's own core formulas in line with that).

## Gap inventory

| # | Gap | Source | Disposition |
|---|---|---|---|
| H1 | Kinetic-energy matrix element (`fillBandedMatrices`, `tise.cpp:616`) hardcodes the atomic-units `1/2` factor instead of `hbar^2/(2*mass)`. | This session's audit | **Fix** |
| H2 | `computeEAcc` already takes `mass` but not `hbar`; its own doc comment derives the formula "in atomic units (hbar=1)." Its one caller (`tise_solver_main.cpp`) passed a hardcoded `mass=1.0` regardless of config. | This session's audit | **Fix** |
| H3 | `matchAsymptotic`'s continuum wavenumber (`k=sqrt(2E)`, two sites) and Sommerfeld parameter (`eta=C/k`) hardcode atomic units. | This session's audit | **Fix** |
| H4 | `matchAsymptotic`'s energy-normalization amplitude `A_E` (two branches: flat-tail, Coulomb) is the atomic-units-only formula from `docs/SDD.md`'s documented target — requires re-deriving the general form (not a mechanical substitution), via the `dk/dE=mass/(hbar^2*k)` Jacobian converting `delta(k-k')`- to `delta(E-E')`-normalization. | This session's audit | **Fix** |
| H5 | `analyticHydrogenEnergy`/`eigenvalueError` (test/legacy-driver helpers) hardcode atomic units; needed generalized for a mass!=1 hydrogen analytic-comparison test. | This session's audit | **Fix** |
| H6 | `solveTISE` (the all-in-one entry point used by unit tests and legacy `TISE/main.cpp`) has no `mass`/`hbar` parameters at all. | This session's audit | **Fix** |
| H7 | `TISE/main.cpp`'s `PARTICLE_MASS=3.75`/`HBAR=1.0` constants were used for the wavepacket shape but never threaded into `solveTISE`/`eigenvalueError` — a latent inconsistency (eigenbasis silently built at `mass=1` regardless). | This session's audit | **Fix** (forced by H6's signature change) |
| H8 | Guard-rail (`tise_solver_main.cpp:159-170`) rejects any non-`1.0` value; needs replacing with a positivity check once every formula above honors the real value. | Part C, `tise-release-readiness-plan.md` | **Fix — last, after H1-H7** |
| — | `radialPotential` (`tise.cpp:55-58`) — dead code kept for reference, not on the production `tise_solver_main.cpp` path. | This session's audit | **Explicitly deferred, no action** |

## The physics

Derived and cross-checked against `docs/SDD.md:555`'s documented atomic-units target formula (confirmed the general form collapses to it exactly at `mass=hbar=1`):

- **Kinetic term**: `(hbar^2/(2*mass)) * bs.integral(fUni, iBs1, iBs2, 1, 1)`.
- **`computeEAcc`**: `hbar^2 * pi^2 / (2*mass*nodeSpacing^2)`.
- **Wavenumber**: `k = sqrt(2*mass*E) / hbar`.
- **Sommerfeld parameter**: `eta = mass*C / (hbar^2*k)` — derived from the radial Schrödinger equation with a `C/r` tail, confirmed the reduced Coulomb ODE (`localK`/phase/WKB-start/Numerov) needs no change, only correct `eta`/`rho` inputs.
- **Energy-normalization amplitude**: general `psi_E = sqrt(2*mass/(pi*hbar^2*k)) * sin(kx+delta)`, from converting the `delta(k-k')`-normalized `sqrt(2/pi)*sin(kx+delta)` via the Jacobian `sqrt(dk/dE) = sqrt(mass/(hbar^2*k))`. Both `A_E` branches (flat, Coulomb) get their leading `2.0/M_PI`-flavored factor replaced with `2.0*mass/(M_PI*hbar*hbar)`.

**Confirmed unaffected** (no formula change): `classifyBoundStates`, `checkWellContainment`, `classifyAsymptote`, `case3WindowFunction`/`evaluateWindowedPotential`, the R-bound check, `buildContinuumState`/`poleTolFraction`, all of `BSpline.cpp`/`.hpp`, grid/knot-placement code, `coulombPhaseShift`/`evaluateCoulombFunctions` (already take `eta`/`rho` as inputs), and all of `time_evolution.cpp`/`.hpp`.

## Part A: Signature changes (closes H1-H6)

New `mass`/`hbar` parameters added as **trailing, defaulted** arguments (`=1.0`), matching this codebase's existing extension convention (`dropSet`, `case3RightR/Delta`, `fineDE`, `coulombLC`, `coulombFarMultiplier/StepSize` were all added to `fillBandedMatrices`/`matchAsymptotic` this same way) — every existing call site keeps compiling unchanged; only new mass!=1 tests need to pass them explicitly:

- `computeEAcc(Real nodeSpacing, Real mass, Real hbar = 1.0)`.
- `fillBandedMatrices(..., Real mass = 1.0, Real hbar = 1.0)` (after `case3RightDelta`).
- `matchAsymptotic(..., Real mass = 1.0, Real hbar = 1.0)` (after `coulombStepSize`).
- `solveTISE(..., Real mass = 1.0, Real hbar = 1.0)` (after `continuumOutputPoints`); threads into its own internal `fillBandedMatrices`/`matchAsymptotic` calls.
- `analyticHydrogenEnergy(int n, int L, Real mass = 1.0, Real hbar = 1.0)` / `eigenvalueError(..., Real mass = 1.0, Real hbar = 1.0)`.

## Part B: Guard-rail replacement (closes H8)

Read `mass`/`hbar` from config (default `1.0` if absent — preserves current behavior for every existing `tests/*.yaml`, none of which set `physics:`), reject non-positive values instead of non-unity ones, thread into the three `tise_solver_main.cpp` call sites (`fillBandedMatrices`, `computeEAcc`, `matchAsymptotic`). Landed **last**, once every formula above had independent unit-test proof of correctness at `mass!=1`/`hbar!=1` — never a window where non-1.0 values were accepted but silently wrong.

## Part C: `TISE/main.cpp` side-fix (closes H7)

Thread `PARTICLE_MASS`/`HBAR` into `solveTISE` and `eigenvalueError` — forced by H6's signature change, and closes a real pre-existing inconsistency (the wavepacket shape assumed `mass=3.75` while the eigenbasis it was evolving silently assumed `mass=1`).

## Part D: Tests

**C++ (`TISE/tests/test_tise.cpp`, GoogleTest):**
- `ComputeEAccTest.ScalesWithHbarSquared` — extends the pre-existing `ScalesInverselyWithMassAndSpacingSquared` mass case.
- `FillBandedMatricesTest.KineticTermScalesWithMassAndHbarOverlapUnchanged` — cross-checks the kinetic block by hand at `mass=2.0, hbar=0.5`; confirms `Smat` (overlap) is unaffected.
- `FreeParticleContinuumGeneralUnitsTest` (3 cases: `EigenvaluesMatchGeneralBoxFormula`, `PhaseShiftMatchesZeroScatteringAtGeneralUnits`, `WrittenWavefunctionMatchesGeneralAnalyticSine`) — a self-contained free-particle system built at `mass=2.0, hbar=0.5`, checked against the general closed forms directly (not by comparison against a second atomic-units run).
- `AnalyticHydrogenEnergyTest.ScalesWithMassAndInverseHbarSquared`, `EigenvalueErrorTest.RespectsMassAndHbar`.
- `SolveTISETest.GroundStateMatchesAnalyticHydrogenAtNonUnitMassAndHbar` — confirms `solveTISE` threads `mass`/`hbar` through end-to-end, at `mass=0.5, hbar=1.5` (chosen so the hydrogenic Bohr radius `a0'=hbar^2/mass=4.5` stays >= the atomic-units `a0=1`, keeping the existing 41-node grid at least as well-resolved — see the "hard finding" below).

**Python (`tests/test_analysis_integration.py`, real-subprocess):**
- New reference config `tests/free_particle_general_units.yaml` — the first `tests/*.yaml` file to ever set a non-default `physics:` block (`mass=2.0, hbar=0.5`). Built on `free_particle.yaml` rather than `harmonic_oscillator.yaml`/`hydrogen.yaml` because both of those bake a coefficient (`0.5*1.0*x^2`'s spring constant, or the centrifugal term) directly into the potential's muparser literal — `V≡0` has nothing to reinterpret, and exercises both the bound-state (kinetic term) and continuum (`k`/`computeEAcc`/`A_E`) halves of the change in one file. `E_max=0.0625` (`= 0.5 * hbar^2/mass`) keeps the same physical `k`-range as the already-characterized `docs/tests/reports/f4e8359/free_particle.md`.
- `TestFreeParticleGeneralUnitsRealSubprocess` (3 cases), modeled on `TestFreeParticleContinuumPhysics`.
- `test_controller_integration.py`'s `test_non_unity_mass_raises_and_leaves_no_partial_files` — premise now false, replaced with `test_non_unity_mass_now_succeeds`, `test_non_positive_mass_raises_and_leaves_no_partial_files`, `test_non_positive_hbar_raises_and_leaves_no_partial_files`.

## Verification

1. `cmake --build TISE/build && ctest --test-dir TISE/build --output-on-failure` — all pass.
2. `pytest tests/` — all pass.
3. Manually run `tests/free_particle_general_units.yaml` via `controller.py`, spot-check output against the closed forms above.

## Implementation status (post-execution)

All of Parts A-D landed, TDD (RED-then-GREEN) followed for every new capability. Final regression state: all 4 C++ suites green (`ctest --test-dir TISE/build`: BSplineTests, UtilsTests, TISETests — 194 individual `TEST`/`TEST_F` cases, including the 8 new cases listed in Part D — TimeEvolutionTests), and all 217 Python tests green (`pytest tests/`). Net Python test-function diff for this phase: `tests/test_controller_integration.py` +3/-1 (the retired `test_non_unity_mass_raises_and_leaves_no_partial_files` replaced by `test_non_unity_mass_now_succeeds`/`test_non_positive_mass_raises_and_leaves_no_partial_files`/`test_non_positive_hbar_raises_and_leaves_no_partial_files`), `tests/test_analysis_integration.py` +3 (`TestFreeParticleGeneralUnitsRealSubprocess`).

**Part A — signature changes: fully done**, exactly as scoped. Every function's `mass`/`hbar` landed as trailing defaulted params; zero existing call sites needed updating.

**Part B — guard-rail replacement: fully done.** One consequence not fully anticipated by the plan: `matchAsymptotic`'s new `mass`/`hbar` trailing params sit after `coulombFarMultiplier`/`coulombStepSize`, which the one production call site (`tise_solver_main.cpp`) previously omitted (relying on their defaults) — reaching the new params meant spelling out `50.0, 0.1` explicitly there. Cosmetic, not a behavior change, but worth noting for anyone diffing that call site.

**Part C — `TISE/main.cpp` side-fix: fully done.**

**Part D — tests: fully done, with one real mid-implementation correction.** The first version of the `A_E`-scaling C++ test asserted `ar.A_E[0]` itself equals the closed-form target amplitude directly — this is wrong: `A_E` is a matching-derived scale factor that also absorbs the resolvent (bound-eigenstate-expansion) construction's own basis-dependent internal normalization (fixed by `B_N`'s coefficient being pinned to `1`, not by a unit-normalization step), so it is not simply predictable from `k`/`mass`/`hbar` alone. The test was rewritten to check the *written* wavefunction (`A_E * psi_num(x)`, via `writeContinuumInfo`) against the closed form instead — mirroring the file's own pre-existing `WrittenWavefunctionMatchesAnalyticSine` pattern — which passed immediately once corrected, confirming the `A_E` formula derivation itself was right; the test's premise, not the implementation, was the bug.

**Hard finding, worth recording for future test design**: the first version of `SolveTISETest.GroundStateMatchesAnalyticHydrogenAtNonUnitMassAndHbar` used `mass=2.0, hbar=0.5` (matching the free-particle test's own choice) and failed by 18% against the analytic formula — not a formula bug, but a genuine test-design flaw: `mass=2.0, hbar=0.5` shrinks hydrogen's Bohr radius `a0'=hbar^2/mass` to `0.125` (an 8x contraction vs. the atomic-units `a0=1`), badly under-resolving the true wavefunction on the same 41-node grid tuned for the atomic-units case. Unlike the free particle's purely geometric box spectrum (mass/hbar-independent by construction), a bound state's required grid resolution depends on the physical length scale set by mass/hbar — a real subtlety this generalization introduces that any future mass/hbar-varying bound-state test needs to account for. Fixed by choosing `mass=0.5, hbar=1.5` instead (`a0'=4.5`, at least as well-resolved as atomic units), after which the test passed at the same tight `1e-4` tolerance as its atomic-units sibling.

**Verification checklist**: items 1-2 done and confirmed (both suites fully green, counts above). Item 3 (manual `controller.py` run) was subsumed by the real-subprocess Python integration tests, which exercise the identical code path.

**Process note**: no git commits were made during this phase — every change is an uncommitted working-tree modification, left for the user to review/stage/commit on their own terms.
