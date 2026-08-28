# Architecture Decision Records

This directory records decisions to **defer** a design alternative rather than adopt it now. Each ADR captures the context that was considered, the decision, and — critically — the condition under which the decision should be revisited. Decisions that were **adopted** live as numbered requirements (REQ-F-/REQ-NF-) in the [SDD](../SDD.md) [§3](../SDD.md#3-requirements), not here.

| ID | Title | Status | Summary |
|---|---|---|---|
| [0001](0001-defer-fedvr-basis.md) | Defer FEDVR as an Alternative Basis | Accepted (deferred) | Finish and validate the B-spline solvers before adding FEDVR as a configurable basis option. |
| [0002](0002-defer-wkb-collocation.md) | Defer WKB/Phase-Angle-Proportional Collocation | Accepted (deferred) | Use uniform + strategic-node knot placement by default; WKB-proportional spacing is a future upgrade. |
| [0003](0003-defer-multiparticle-3d-extension.md) | Defer Multi-Particle / 3D Extension | Accepted (deferred) | Stay strictly 1D, single-particle for this project; multi-particle/3D is future pedagogical work. |
| [0004](0004-defer-visualization-plot-parameters.md) | Defer Visualization Plot-Parameter Schema | Accepted (deferred) | Keep `visualization` boolean-only until the Analysis module's detailed design (SDD [§5.4](../SDD.md#54-analysis-module)) is complete. |
| [0005](0005-defer-analysis-output-artifact-format.md) | Defer Analysis Output-Artifact Format | Accepted (deferred) | Analysis produces no output artifact (no files, no plots, no `data/analysis/`) until REQ-F-060 is computable in SDD [§10.2](../SDD.md#102-phased-implementation-sequence) Phase 7/8; distinct from ADR-0004's narrower plot-*parameter* scope. |
| [0006](0006-tdse-propagator-choice.md) | Defer Crank-Nicolson-Family TDSE Propagators | Accepted (deferred) | Use a symmetric split-operator propagator for the TDSE driven case (SDD [§5.3.3](../SDD.md#533-internal-design) Figure 8); Crank-Nicolson-family and hybrid alternatives are compared and deferred/rejected. |
| [0007](0007-defer-bound-state-filtering-tise-eigenvalue-output.md) | Defer Bound-State Filtering of TISE Eigenvalue/Eigenvector Output | Accepted (deferred) | `tise_solver` writes all computed eigenvalues/eigenvectors to `data/tise/`; filtering to bound-only states via `classifyBoundStates` is left to downstream consumers. |
| [0008](0008-defer-unify-solvetise-tise-solver-orchestration.md) | Defer Unifying `solveTISE`/`tise_solver` Orchestration | Accepted (deferred) | `solveTISE` gains strategic node placement (A4/A4b); `tise_solver_main.cpp` keeps its own independent, unwired orchestration rather than extracting a shared core-solve function now. |

See also [SDD Appendix B](../SDD.md#b-open-design-questions) for design questions that are still genuinely unresolved (no decision, not even a decision to defer).

## Format

Each ADR uses the same five fields: **Status**, **Date**, **Context**, **Decision**, **Consequences** (including an explicit revisit trigger), and **Source**.
