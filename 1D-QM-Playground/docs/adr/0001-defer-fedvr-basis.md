# ADR-0001: Defer FEDVR as an Alternative Basis

- **Status:** Accepted (deferred)
- **Date:** 2026-07-03

## Context

Finite-Element Discrete Variable Representation (FEDVR) is mathematically a limiting case of the B-spline basis this project already uses: it is what you get when every interior knot is pushed to multiplicity $k-1$, splitting the domain into independent finite elements, and the local B-spline basis within each element is replaced by its Discrete Variable Representation (a basis that diagonalizes the position operator $\hat{x}$ within the element).

FEDVR offers two properties the current B-spline basis does not:

- The potential energy matrix $V$ becomes **diagonal** (versus banded for B-splines).
- The overlap matrix becomes the **identity** (versus banded for B-splines), turning the generalized eigenvalue problem $\mathbf{H c} = E\mathbf{S c}$ into a standard eigenvalue problem $\mathbf{Hc} = E\mathbf{c}$.

Together these make each application of $H$ to a state vector $O(N)$, which is what makes Krylov/Chebyshev-type explicit TDSE propagators attractive for large systems — an efficiency class not accessible to the current banded B-spline representation without first diagonalizing.

Because FEDVR is constructed *within* B-spline space (same knot vector and assembly routines; only the local per-element basis changes via the DVR transform), the stakeholder noted that FEDVR support could emerge as a comparatively low-effort side product once the B-spline infrastructure is solid, rather than requiring a parallel implementation effort now.

Full exploration, including the matrix-structure comparison and the migration path, is retained in `docs/planning/fedvr-exploration.md`.

## Decision

Defer FEDVR implementation. Complete and validate the B-spline TISE and TDSE solvers first. Revisit FEDVR afterward as a configurable basis option, e.g. a `basis.type: fedvr` flag in `config.yaml`, alongside the existing B-spline basis — not a replacement for it.

## Consequences

- The current and planned `config.yaml` schema (SDD [§6.1](../SDD.md#61-configuration-schema)) has no `basis.type` field; B-splines are the only basis for the lifetime of this decision.
- No architectural provision is made now for a diagonal-potential/identity-overlap code path, DVR quadrature points, or Krylov/Chebyshev propagators.
- **Revisit trigger:** the B-spline TISE and TDSE solvers are validated against the physics benchmarks in SDD [§9.3](../SDD.md#93-verification-and-validation). At that point, FEDVR should be re-evaluated as a `basis.type` extension per the migration path already documented in `docs/planning/fedvr-exploration.md` §6.

## Source

`docs/planning/fedvr-exploration.md` (full exploration); referenced from SDD [§4.1](../SDD.md#41-architectural-style-and-rationale) and [§5.2.3](../SDD.md#523-internal-design); corresponds to the excluded Figure 11 in `docs/diagrams.md`.
