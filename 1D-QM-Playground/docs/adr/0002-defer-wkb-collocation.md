# ADR-0002: Defer WKB/Phase-Angle-Proportional Collocation

- **Status:** Accepted (deferred)
- **Date:** 2026-07-03

## Context

Node (knot) placement for the B-spline basis breaks into two independent concerns: **strategic placement**, driven by structural features of the potential (degeneracies at deltas/steps, removals at singularities — always applied automatically regardless of this decision), and **density distribution**, driven by accuracy requirements across the domain.

Two density-distribution schemes were considered:

- **Uniform** — equal spacing everywhere. Simplest to implement; potential-agnostic.
- **WKB-proportional** — node density proportional to the local classical momentum at a reference energy $E$, so that resolution concentrates wherever the wavefunction is expected to oscillate fastest:

$$n(x) = \alpha\sqrt{2m(E - V(x))}$$

WKB-proportional spacing is the more principled general-purpose scheme (it adapts to any potential without Coulomb-specific assumptions), but it requires choosing a reference energy $E$ and adds implementation complexity beyond what a first working solver needs.

## Decision

Use a uniform B-spline basis as the default node-density scheme, augmented only by the strategic nodes already dictated by potential structure (REQ-F-050). WKB-proportional spacing is not implemented in the initial version.

## Consequences

- The initial `bspline` configuration (SDD [§6.1](../SDD.md#61-configuration-schema)) exposes no reference-energy or WKB-mode field; density is uniform by construction, strategic nodes are automatic.
- Continuum and bound-state accuracy in regions of rapid oscillation is bounded by uniform-grid resolution until this is revisited (see the basis accuracy ceiling $E_\text{acc}$ in REQ-F-040).
- **Revisit trigger:** the uniform + strategic-node scheme is found to be an accuracy or efficiency bottleneck in production runs (e.g., continuum phase shifts inaccurate near $E_\text{acc}$, or excessive node counts needed to compensate). At that point, implement WKB-proportional density as a selectable `bspline` density mode with a user- or auto-selected reference energy.

## Source

`docs/planning/architecture-06-20.md`, "Collocation scheme" (stakeholder feedback, 2026-07-03); referenced from SDD [§5.2.3](../SDD.md#523-internal-design) and REQ-F-050.
