# ADR-0007: Defer Bound-State Filtering of TISE Eigenvalue/Eigenvector Output

- **Status:** Accepted (deferred)
- **Date:** 2026-08-27 (identified and decided while wiring `tise_solver` to real TISE
  physics, SDD §7.2.1/§6.3, closing the Controller↔TISE interface gap)

## Context

SDD §6.3 documents `eigenvalues.dat`/`eigenvectors.dat`'s file format but not a
selection criterion for which (or how many) computed states to include. The only
existing precedent, `main.cpp`'s `EigenState_NNN`-writing loop, breaks at the first
eigenvalue whose error against `tise::analyticHydrogenEnergy` exceeds
`config.yaml`'s `tise.error_threshold` — hardcoded to the hydrogen analytic
spectrum, meaningless for the arbitrary piecewise potentials this generalized
solver now supports.

`tise::classifyBoundStates(EigenResult, threshold)` (REQ-F-020) is the library's
general-purpose tool for a related but distinct question: it classifies computed
states as bound vs. continuum by a plain energy threshold (every existing test
uses `0.0`), returning a proven-prefix `isBound`/`nBound`. Two options were
considered for `tise_solver`'s output: (a) write only the first `nBound` states
per `classifyBoundStates(er, 0.0)`, or (b) write all `nEn` computed states
regardless of sign.

## Decision

`tise_solver` writes **all `nEn` computed eigenvalues/eigenvectors** to
`eigenvalues.dat`/`eigenvectors.dat`. No bound-state filtering is applied at the
`tise_solver` level. `config.yaml`'s `tise.error_threshold` field is not consumed
by `tise_solver` (it remains meaningful only to the separate, hydrogen-specific
`H-BoundStates` demo executable).

## Consequences

- `eigenvalues.dat`/`eigenvectors.dat` include `E >= 0` "states" that are really
  finite-box discretization artifacts (their energy/wavefunction shape depends on
  `bspline.domain`/`n_nodes`, not the physical problem) alongside genuine bound
  states. A consumer wanting bound-states-only must filter itself — the same
  `classifyBoundStates` function remains available as a C++ library call for this;
  no Python equivalent is exposed to `analysis.py` today.
- `tise_solver` makes no bound/continuum physics judgment call internally, keeping
  it a thinner, more mechanical translation of "solve and report" — consistent
  with this pass also not wiring per-state `checkWellContainment` warnings (that
  would require the same classification this ADR defers).
- File size is `nEn` rows rather than a filtered `nBound` subset — a modest
  difference for the current `config.yaml` (`nEn=59`), potentially more relevant
  for much larger bases.
- **Revisit trigger:** if a downstream consumer (Analysis, once REQ-F-060 quantities
  land in a later phase, or a human inspecting output) needs pre-filtered
  bound-only data and re-deriving it downstream proves inconvenient or
  threshold-convention-inconsistent, OR if output size/transfer becomes a real
  concern for large bases — revisit whether `tise_solver` itself should filter via
  `classifyBoundStates(er, 0.0)`, or add a separate, explicitly-named output
  (e.g. `bound_eigenvalues.dat`) alongside the full one rather than replacing it.

## Source

SDD §6.3 (`data/tise/` persistent storage format), §7.2.1 (Controller→TISE
contract); identified while wiring `tise_solver` to real physics, 2026-08-27.
Cross-reference: `tise::classifyBoundStates` (`TISE/tise.hpp`) — the
not-invoked-here tool for the deferred alternative.
