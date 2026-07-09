# ADR-0003: Defer Multi-Particle / 3D Extension

- **Status:** Accepted (deferred)
- **Date:** 2026-07-03

## Context

The dynamics-operator scope decision (REQ-F-010) restricts this project to a single particle in 1D, with the Hamiltonian fully determined by a user-supplied $V(x)$ and the exposed operator set limited to $\hat{H}$, $\hat{x}$, and $\hat{p}_x$. This keeps the physics and the code simple: one discrete symmetry (parity), at most two escape channels, no angular momentum or centrifugal terms, no coupling between particles.

A natural pedagogical extension, raised during design discussion, would be to add a second interacting particle and/or extend to 3D. This opens substantially more physics (entanglement, correlation, exchange symmetry, decoherence) but also substantially more implementation complexity (multi-particle Hilbert space construction, antisymmetrization, higher-dimensional bases) and more demanding simulations.

## Decision

Explicitly out of scope for this project. The priority is a clean, self-contained 1D single-particle solver carried through to publication. The multi-particle/3D extension may be pursued afterward depending on where the project stands at completion — it is future work, not a near-term requirement.

## Consequences

- No architectural provision is made for multi-particle Hilbert spaces, antisymmetrization, or bases beyond 1D single-particle $\hat{x}$.
- REQ-F-010 is scoped explicitly to 1D, single-particle, and should not be reinterpreted to implicitly allow multi-particle extensions without revisiting this ADR first.
- **Revisit trigger:** the core 1D single-particle solver (TISE + TDSE) reaches publication-ready completion, and there is appetite and time remaining to pursue the pedagogical extension.

## Source

`docs/planning/architecture-06-20.md`, "Operators" → "Future directions (deferred)" (stakeholder feedback, 2026-07-03); referenced from SDD [§3.1](../SDD.md#31-functional-requirements) (REQ-F-010) and [§5.1](../SDD.md#51-controller) note.
