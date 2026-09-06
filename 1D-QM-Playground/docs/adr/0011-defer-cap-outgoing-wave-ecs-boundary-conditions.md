# ADR-0011: Defer CAP / Outgoing-Wave (Siegert) / Exterior-Complex-Scaling Boundary Conditions

- **Status:** Accepted (deferred)
- **Date:** 2026-09-04 (formalized during a post-implementation cleanup pass;
  the question itself has been open, unformalized, since the SDD's initial
  population — see Source)

## Context

SDD [§5.2.3](../SDD.md#523-internal-design)'s "Initial analysis" (`docs/SDD.md:1325`–1329)
identifies three approaches for resonance/ionization boundary conditions
beyond the plain Dirichlet wall REQ-F-030 formalizes:

- **Outgoing-wave (Siegert) boundary conditions**: replace the Dirichlet
  condition at $r_\text{max}$ with the requirement that $\psi$ be a purely
  outgoing wave, $\psi(r_\text{max}) \sim e^{ikr}$. This yields complex
  eigenvalues $E = E_r - i\Gamma/2$ directly (resonance position and width).
  Rigorous implementation requires exterior complex scaling (ECS).
- **Exterior complex scaling (ECS)**: rotate the radial coordinate into the
  complex plane beyond some radius, $r \to r_0 + (r-r_0)e^{i\theta}$ for
  $r > r_0$. This makes outgoing (divergent) Siegert states
  square-integrable in the rotated coordinate, turning the resonance search
  into an ordinary (complex-symmetric, not Hermitian) eigenvalue problem.
- **Complex absorbing potential (CAP)**: add a smooth imaginary potential
  $-iW(r)$ near the outer wall that damps outgoing flux before it reflects
  off the Dirichlet boundary. `docs/SDD.md:1327` calls this "the more
  tractable near-term option," and `docs/SDD.md:1329`'s recommendation
  names it "the recommended next addition for continuum and ionization
  calculations, since it does not require changes to the real-valued
  eigensolver."

`docs/SDD.md:1376` records that none of the three were revisited during the
2026-07-03 stakeholder feedback pass, and states explicitly: *"unlike
FEDVR, WKB collocation, the multi-particle extension, or the visualization
schema (ADR-0001 through ADR-0004), no decision has been made to defer it
either. The next time this is discussed with the stakeholder, it should be
promoted to either a new REQ (if adopted) or a new ADR (if a conscious
decision to defer is made)."* SDD [§12.B](../SDD.md#b-open-design-questions)
(Open Design Questions) has carried this item as genuinely unresolved ever
since. `docs/planning/tise-task-breakdown.md` §6 independently lists it as
out of scope for the Engineer A/B implementation phases, citing the same
SDD passage. That promotion — to REQ or ADR — has still not happened; this
ADR is that promotion, on the "defer" side.

All three options share a common prerequisite the current implementation
does not have: none of `DSBGV` (LAPACK's real-symmetric banded generalized
eigensolver, used throughout `TISE/tise.cpp`), the banded
`H`/`S`-matrix-assembly code, or `bspline::BSpline` support complex-valued
matrix elements or a complex-valued (or complex-coordinate) eigenproblem.
Adopting any of the three is therefore not a small, additive change — it
touches the Hamiltonian-assembly and eigensolver layer the rest of the TISE
implementation is built on.

## Decision

Do not adopt CAP, outgoing-wave/Siegert boundary conditions, or ECS now.
The solver continues to support only real-valued Hamiltonians via `DSBGV`,
with Dirichlet walls (or Case 1–3 asymptote matching, per REQ-F-030 and
ADR-0010) at the domain's unbounded edges. This ADR formalizes, as a
conscious deferral, the open question SDD §12.B has carried since the SDD
was first populated — it does not change any code.

## Consequences

- Resonance calculations (complex eigenvalues $E_r - i\Gamma/2$) and
  ionization-flux studies that need absorbing/outgoing boundary conditions
  are not supported. Any such work must currently either accept the
  finite-box pseudostate approximation the existing real-eigensolver
  pipeline already produces (bound + discretized continuum, per SDD
  [§5.2.3](../SDD.md#523-internal-design)'s Figure 7 tree), or be done
  outside this codebase.
- The companion planning doc,
  `docs/planning/complex-boundary-conditions.md`, gives a comparative
  overview of the three options and a recommendation for which to
  prototype first, so a future implementer is not starting from the SDD's
  brief paragraph-level discussion alone.
- SDD [§12.B](../SDD.md#b-open-design-questions)'s entry for this question
  should be annotated, per its own stated convention ("Entries... are never
  deleted as they're resolved — they stay verbatim, annotated with a Status
  line and a pointer to whichever REQ or ADR resolved them"), to point at
  this ADR. That annotation is a small, targeted SDD pointer addition made
  alongside this ADR (see `docs/SDD.md:1376`/`:1329`) — a full SDD sync of
  §12.B's phrasing is out of scope for this change.
- **Revisit trigger:** if a future project phase needs resonance widths,
  ionization rates via outgoing flux, or any other genuinely complex-valued
  spectral quantity — the natural next step given the project's existing
  TDSE roadmap, which already anticipates ionization-probability analysis
  (SDD [§2.2](../SDD.md#22-goals-and-objectives)) — revisit this decision.
  The companion planning doc's recommendation is to prototype CAP first
  (smallest change: an added imaginary diagonal term, still solvable via a
  complex-Hermitian or perturbative-width extension, without touching the
  real-coordinate B-spline basis or knot grid), before attempting ECS
  (which requires complexifying the radial coordinate and grid itself).

## Source

`docs/SDD.md:1325`–1329 (Initial analysis: outgoing-wave, CAP), `:1376`
(Status: genuinely unresolved, "should be promoted to... a new ADR"),
[§12.B](../SDD.md#b-open-design-questions) (Open Design Questions);
`docs/planning/tise-task-breakdown.md` §6 ("Explicitly Out of Scope").
See `docs/planning/complex-boundary-conditions.md` for the comparative
overview and prototyping recommendation.
