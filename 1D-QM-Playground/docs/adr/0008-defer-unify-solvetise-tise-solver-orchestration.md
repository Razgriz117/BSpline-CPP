# ADR-0008: Defer Unifying `solveTISE`/`tise_solver` Orchestration

- **Status:** Accepted (deferred)
- **Date:** 2026-08-28 (identified while wiring strategic node placement (A4/A4b)
  into `solveTISE`, SDD §5.2.3)

## Context

Two separate C++ entry points independently orchestrate the same
fill-matrices → diagonalize → (optionally) build-continuum sequence:
`tise::solveTISE` (`TISE/tise.cpp`, called by `main.cpp`/`H-BoundStates`) and
`tise_solver_main.cpp`'s `main()` (the config-driven `tise_solver` CLI). They
exist separately because their output destinations and gating differ:
`solveTISE` always builds continuum output and writes `phase_shifts.dat`/
`continuum_state_NNN.dat` relative to the current working directory;
`tise_solver_main.cpp` gates continuum construction on
`config.yaml`'s `tise.continuum.enabled` and writes everything under an
explicit `--output-dir`, plus the `eigenvalues.dat`/`eigenvectors.dat`/
`hamiltonian.dat`/`overlap.dat` files `solveTISE` never produces.

Wiring strategic grid construction and singular-join B-spline removal into
`solveTISE` (this session's task) adds a *third* independently-duplicated
concern `tise_solver_main.cpp` would need to hand-sync if it ever needed the
same capability: it would have to re-derive the same
`detectPotentialStructure`/`strategicKnotsFromJoins`/`buildStrategicRadialGrid`/
`bSplinesTouchingX` orchestration `solveTISE` now contains, rather than
reusing it.

The natural alternative — extract a shared, file-write-free "core solve"
(returning something like `solveTISE`'s new `SolveTISEResult`) that both
`solveTISE` and `tise_solver_main.cpp` call, with each entry point handling
only its own output/gating on top — was considered as part of resolving
`solveTISE`'s Gap 3 (how it exposes the grid/basis/drop-set it used to
callers). It is the same shape as a rejected "caller-constructed
`SolveContext`" alternative to `SolveTISEResult`.

## Decision

Do not unify `solveTISE`/`tise_solver_main.cpp`'s orchestration now.
`tise_solver_main.cpp` is left untouched by the A4/A4b wiring work; it
continues to build its own uniform grid via its own direct calls to
`fillBandedMatrices`/`solveGeneralizedEigenproblem`/`buildContinuumState`/
`matchAsymptotic`, independently of whatever `solveTISE` now does
internally.

## Consequences

- `tise_solver`'s config-driven output does **not** benefit from strategic
  node placement or singular-point B-spline removal — it remains on the
  plain uniform grid / classic-adjacent `{1}` drop-set convention it already
  used. This is a real capability gap for general (non-hydrogen,
  structure-bearing) potentials run through `tise_solver`, not merely a
  cosmetic duplication.
- Any future change to `fillBandedMatrices`/`buildContinuumState`/
  `matchAsymptotic`'s calling convention (as this session's continuum
  boundary-coupling generalization already required) must be applied
  consistently to both `solveTISE` and `tise_solver_main.cpp` by hand; there
  is no shared code path enforcing that they stay in sync.
- **Revisit trigger:** if `tise_solver` ever needs strategic placement too
  (a real, foreseeable need — it's the general-purpose, config-driven
  entry point the project is otherwise growing), or if the two copies drift
  and cause a bug from exactly this kind of hand-sync failure — extract the
  shared core-solve function (returning `SolveTISEResult` or equivalent)
  instead of hand-syncing a fourth duplicated concern.

## Source

`TISE/tise.cpp` (`solveTISE`), `TISE/tise_solver_main.cpp`; identified while
wiring A4/A4b into `solveTISE`, 2026-08-28. See
`docs/planning/engineer-a-plan-A4-wiring-design.md` for the wiring design
this ADR's context arose from.
