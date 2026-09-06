# BUG-0001: `writeContinuumInfo` Evaluated Continuum States in the Wrong Coefficient Basis

- **Status:** Fixed
- **Date introduced:** 2026-08-12 (`9c318c1`, "Added energy grid ..., added function to write to file as per B.3.4, added unit tests")
- **Date fixed:** 2026-08-30 (`9d0f04b`, "test: e2e known solution config files")
- **Report date:** 2026-09-04

## Symptom

Continuum-state wavefunction plots produced by `analysis.py` (from
`data/tise/continuum_state_*.dat`) looked physically wrong at every energy,
regardless of the potential: a spike near the wall and near-zero amplitude
at both edges, instead of an oscillatory scattering wave. For the
free-particle case the numeric wavefunction should overlay the analytic
`A·sin(kx)` solution; instead it was off by roughly a factor of 3 in
amplitude and the wrong shape entirely.

## Root Cause

`buildContinuumState` (`TISE/tise.cpp`) returns each continuum state's
coefficients in the **confined-eigenstate basis** `{φₙ} ∪ {B_N}` — i.e.
expansion coefficients over the bound-problem's eigenvectors, plus the
boundary B-spline `B_N` — not raw B-spline coefficients. `bspline::BSpline::eval`
requires raw B-spline coefficients.

`writeContinuumInfo`, added in commit `9c318c1` (2026-08-12), called
`bs.eval(x, states[i].data(), n)` directly on the confined-basis
coefficients, silently treating them as B-spline coefficients.

`matchAsymptotic` had the identical bug at the time, and was fixed
separately in commit `76d0d8c` (2026-08-21, "fixed matchAsymptotic so that
we can replicate square well phase shift") by inserting a basis-transform
loop (project the confined-basis coefficients back through
`eigen.vectors` into B-spline coefficient space) before its own `bs.eval`
calls. That fix was never propagated to `writeContinuumInfo` — the two
functions needed the identical transform but only one got it. When branch
`4ff1005` ("merge gracedev's continuum-state construction (B1-B5)",
2026-08-23) merged into `main`, it carried `matchAsymptotic`'s fix and
`writeContinuumInfo`'s still-broken, untransformed version together, so the
bug shipped to `main` in that merge and remained there for the plotting
path.

## Fix

Commit `9d0f04b` (2026-08-30):

- Factored the basis-transform out of `matchAsymptotic` into a shared
  helper, `continuumStateToBSplineCoeffs` (`TISE/tise.cpp`), specifically
  so the transform exists in exactly one place and can't diverge between
  callers again.
- Applied that helper inside `writeContinuumInfo` before each `bs.eval`
  call, replacing the direct `bs.eval(x, states[i].data(), n)` call on
  unconverted coefficients.
- Threaded `eigen` and `dropSet` through `writeContinuumInfo`'s signature
  (both trailing parameters; `dropSet` defaulted to the classic `{1}`
  drop-set) and updated both call sites (`solveTISE` in `tise.cpp`, and
  `tise_solver_main.cpp`'s own orchestration) to pass them through.
- Added a pole-proximity warning in `buildContinuumState`: its
  `c_n(E) = (...)/(E - Eₙ)` sum has no guard against a requested continuum
  energy landing close to a confined eigenvalue (bound, or an
  ADR-0007 unfiltered box-discretization artifact). New trailing, defaulted
  parameters `poleTolFraction=0.1, warnOut=std::cerr`; warns when a grid
  point is within `poleTolFraction * (grid[1] - grid[0])` of any
  `eigen.values[i]`. Wired into `tise_solver_main.cpp`'s existing
  `addWarning`/`warnings.json` path (same pattern as the pre-existing
  E_acc warning). Diagnostic only — does not change any computed value.

## Verification

- C++ (`TISE/tests/test_tise.cpp`): new `FreeParticleContinuumTest` fixture
  checks eigenvalues against the infinite-square-well analytic values,
  `matchAsymptotic`'s phase shift δ≈0, and — the core regression —
  `writeContinuumInfo`'s output against the analytic `A·sin(kx)` solution.
  A self-consistency addition to `SquareWellPhaseShiftTest`, and two new
  `BuildContinuumStateTest` cases cover the pole warning.
- Python (`tests/test_analysis_integration.py`): `TestFreeParticleContinuumPhysics`
  runs the real `tise_solver` subprocess end-to-end and checks eigenvalues,
  phase shifts, and continuum wavefunction shape against analytic values;
  `TestHarmonicOscillatorBoundStates` checks bound eigenvalues against
  `(n+½)ω`.
- Visual confirmation: `data/tise/continuum_*.png` went from the
  wall-adjacent-spike/near-zero-at-both-edges shape (present at every
  energy, for every potential tested) to a clean physical scattering wave;
  the free-particle numeric wavefunction now overlays the analytic
  `A·sin(kx)` solution essentially exactly.

A continuum state's overall sign is not physically meaningful (the same
freedom as an eigenvector's sign, and dependent on which sign the
underlying LAPACK solve happens to return at a given energy), so both test
suites compute an aggregate correlation/dot-product sign against the
analytic solution first and compare under that one consistent sign, rather
than asserting an exact `psi == +A·sin(kx)` match.

## Source

`TISE/tise.cpp` (`continuumStateToBSplineCoeffs`, `writeContinuumInfo`,
`buildContinuumState`), `TISE/tise.hpp`, `TISE/tise_solver_main.cpp`;
regression tests in `TISE/tests/test_tise.cpp` and
`tests/test_analysis_integration.py`. Introduced by `9c318c1`; sibling fix
in `matchAsymptotic` by `76d0d8c`; merged into `main` unresolved by
`4ff1005`; fixed by `9d0f04b`.
