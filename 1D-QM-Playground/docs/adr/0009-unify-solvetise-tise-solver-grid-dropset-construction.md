# ADR-0009: Unify `solveTISE`/`tise_solver` Grid + Drop-Set Construction

- **Status:** Accepted (implemented) — **supersedes ADR-0008**
- **Date:** 2026-08-28 (TISE release-readiness pass, ahead of TDSE work)

## Context

ADR-0008 accepted, as a real and known gap, that `tise_solver_main.cpp` (the
config-driven CLI production binary) always built a plain uniform grid and
hardcoded the classic `{1}` B-spline drop-set, unlike `solveTISE`
(`H-BoundStates`' entry point), which had strategic node placement
(REQ-F-050) and A4b interior-singular-B-spline removal. ADR-0008 named its
own revisit trigger explicitly: *"if `tise_solver` ever needs strategic
placement too... or if the two copies drift and cause a bug."*

A release-readiness audit (`docs/planning/tise-release-readiness-plan.md`,
Part A) confirmed this gap is not cosmetic: a config.yaml with an interior
Step/StitchedKink/Singular potential join fed to the real `tise_solver`
binary produces **silently degraded eigenvalues with `EXIT_SUCCESS` and zero
warning** — the exact failure mode ADR-0008 flagged as a real capability
gap, not merely reduced accuracy. Given the goal of shipping TISE
end-to-end before TDSE work begins, this was escalated from "accepted,
deferred" to "fix now."

## Decision

Extract `solveTISE`'s grid/dropset-construction logic (previously inlined
only there: `detectPotentialStructure` → `strategicKnotsFromJoins` →
`buildStrategicRadialGrid` → per-join `bSplinesTouchingX` removal →
`fillDropSet`/`nEnBound` computation) into a new shared function,
`tise::buildStrategicGridAndDropSet` (`TISE/tise.hpp`/`tise.cpp`), returning
a `StrategicGridResult{grid, bs, nBSplines, fillDropSet, nEnBound,
rightEdgeSingular}`.

- `solveTISE` now calls this shared function instead of inlining the same
  logic; everything after (continuum construction, its own file-writing) is
  unchanged.
- `tise_solver_main.cpp` now calls the same function instead of
  `buildUniformRadialGrid` + a hardcoded `{1}` drop-set, and threads the
  resulting `fillDropSet`/`nBSplines` through every downstream call that
  previously relied on the `nullopt`→classic-`{1}` default
  (`fillBandedMatrices`, `buildContinuumState`, `matchAsymptotic`,
  `writeContinuumInfo`, the `eigenstate_NNN.dat` loop's
  `eigenstateCoefficients` calls, and `writeEigenvectors`, which gained a
  `dropSet` parameter it didn't have before). It also now surfaces
  `rightEdgeSingular` through its `warnings.json` path — previously only
  `solveTISE` could detect this (via a bare `cerr` print `tise_solver_main.cpp`
  had no access to at all), and neither production binary's users could see
  it in a machine-readable form.
- `computeEAcc`'s `nodeSpacing` input in both orchestration paths now uses
  `minInterNodeGap(grid)` instead of a flat uniform-spacing formula — correct
  once a genuinely non-uniform strategic grid is possible in either path
  (previously `minInterNodeGap` had zero production callers anywhere,
  including inside `solveTISE` itself, despite being the one function that
  actually builds non-uniform grids).

A potential with no detectable Step/StitchedKink/Singular structure produces
an unchanged uniform grid and the classic `{1}` drop-set in both entry
points, byte-identical to pre-ADR-0009 behavior — verified directly against
the real `config.yaml`'s hydrogen potential (edge-singular, not interior):
identical eigenvalues and identical warnings before/after this change.

## Consequences

- `tise_solver`'s config-driven output now genuinely benefits from strategic
  node placement and interior singular-B-spline removal, matching
  `solveTISE`/`H-BoundStates` and REQ-F-050's literal claim ("strategic node
  placement... as the default collocation scheme").
- Only one copy of the grid/dropset construction logic exists now; the
  hand-sync risk ADR-0008 named is resolved by construction, not by
  discipline.
- `solveTISE`'s own file-writing/continuum-construction responsibilities and
  `tise_solver_main.cpp`'s own config-driven output/gating remain
  deliberately separate — this ADR only unifies the grid/dropset
  construction half, not the two entry points' full orchestration (they
  still differ in output destination, continuum gating, and warning
  surfacing, by design).
- Real-subprocess regression coverage added:
  `tests/interior_singularity.yaml` (genuine interior `1/(x-20)`-style
  singularity) and `tests/right_edge_singularity.yaml`
  (`1/(rMax-x)`-style right-edge singularity), both exercised through the
  real `tise_solver` binary in `tests/test_analysis_integration.py` — these
  are what would catch a regression to the pre-ADR-0009 behavior or a missed
  drop-set-threading site.

## Source

`TISE/tise.hpp`/`tise.cpp` (`buildStrategicGridAndDropSet`, `solveTISE`),
`TISE/tise_solver_main.cpp`; `docs/planning/tise-release-readiness-plan.md`
Part A. Supersedes `docs/adr/0008-defer-unify-solvetise-tise-solver-orchestration.md`.
