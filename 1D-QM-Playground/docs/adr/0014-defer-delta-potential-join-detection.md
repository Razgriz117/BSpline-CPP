# ADR-0014: Defer Delta-Potential Join Detection in Strategic Node Placement

- **Status:** Accepted (deferred)
- **Date:** 2026-09-11 (formalized during a node-placement capability review
  against the original design docs, prompted by
  `docs/planning/engineer-a-plan-cleanup.md` item 8, "Verify node placement
  implementation is at least as capable as what we discussed in original
  design meetings"; the gap itself has existed, unformalized, since A4 was
  implemented)

## Context

Both design-time sources for REQ-F-050's required-treatment table list
**four** potential-join types, not three:

- `docs/planning/architecture-06-20.md`, "Strategic node placement"
  (2026-07-03 stakeholder feedback): *"Delta potential $\delta(x-x_0)$ —
  Pile up degenerate knots at $x_0$ to capture the discontinuity in the
  first derivative of $\psi$."*
- `docs/planning/tise-task-breakdown.md`'s A4 task spec reproduces the same
  row: *"Delta potential $\delta(x-x_0)$ — Pile up degenerate knots at
  $x_0$ (discontinuity in $\psi'$)."*

The knot-multiplicity formula for this case was even derived and verified:
`docs/planning/engineer-a-plan-A4.md:381` uses the delta case as a
**self-consistency cross-check** of the Step/StitchedKink formulas — "a
delta potential forces a ψ' jump (`n=1`) → `extra = order-2`... matching
the standard delta-matching condition" — and ADR-0012 confirms the
underlying continuity formula ($C^{k-\nu-1}$) against
`PHY5606_F25_Bsplines_v2.pdf`. The math this case needs is therefore
already worked out and independently verified.

What was never built is the *detection* side. `TISE/tise.hpp:521-527`'s
`JoinType` enum has exactly three actionable cases plus `Continuous`:

```cpp
enum class JoinType
{
    Continuous,
    Step,
    StitchedKink,
    Singular
};
```

There is no `Delta` case, and `detectPotentialStructure`
(`TISE/tise.cpp`) has no code path that could produce one — it classifies
joins by probing the value/slope of the piecewise potential expression on
either side of a boundary (`isSingularApproaching` plus one-sided finite
differences), which can only ever see `Continuous`/`Step`/`StitchedKink`/
`Singular` behavior in an ordinary function, never a distributional
(delta-function) potential.

This is not purely an oversight, however: `config.yaml`'s `potential`
field is a list of `{domain, function}` pieces evaluated by `muparser`
(`docs/SDD.md:1123`, `tise::evaluateFunction`) — an ordinary closed-form
expression parser with no representation for a Dirac delta at all. A
literal delta potential cannot be expressed in the current config schema
in the first place, so the detector has nothing to detect even if it were
extended. No ADR or planning document connects these two facts explicitly;
the gap was found only during this capability review, not during A4's
original implementation or any later audit.

## Decision

Do not implement delta-potential join detection now. `JoinType` stays
limited to `Continuous`/`Step`/`StitchedKink`/`Singular`; the derived
`extra = order-2` formula remains unused in production, recorded only as
the A4 cross-check. This ADR formalizes, as a conscious deferral, a gap
between the original design's four-row table and the three-case
implementation that was never previously written down.

## Consequences

- A config that expresses a delta-like potential only approximately (e.g.
  a narrow, tall `Step`-bounded rectangular barrier, the closest
  expressible analogue) gets `Step` treatment (`extra = order-3`) at each
  of its two edges, not the dedicated `extra = order-2` single-knot
  treatment a true delta would call for. This is a pre-existing behavior,
  unchanged by this ADR — no regression results from deferring rather than
  implementing.
- If the config schema ever gains a way to express distributional
  potentials (e.g. a dedicated `delta` piece type alongside `function`),
  `JoinType::Delta` and its `extra = order-2` treatment can be added to
  `detectPotentialStructure`/`strategicKnotsFromJoins` directly — the
  formula and its independent verification (ADR-0012, the A4 cross-check)
  already exist and would not need to be re-derived.
- **Revisit trigger:** if the config schema is extended to represent delta
  (or other distributional) potentials, or if a stakeholder specifically
  requests delta-potential support for a teaching example — revisit this
  decision and wire the already-derived formula into `JoinType`.

## Source

`docs/planning/architecture-06-20.md` "Strategic node placement"
(2026-07-03 stakeholder feedback table); `docs/planning/tise-task-breakdown.md`
A4 task spec (same table); `docs/planning/engineer-a-plan-A4.md:381` (the
`extra = order-2` cross-check derivation); `docs/adr/0012-defer-a4-discontinuity-order-verification.md`
(literature verification of the underlying continuity formula);
`TISE/tise.hpp:521-527` (`JoinType` enum, no `Delta` case);
`docs/SDD.md:1123` (muparser-based `potential.function` expression parser,
no distributional-potential support); `docs/planning/engineer-a-plan-cleanup.md`
item 8 (the capability review that surfaced this gap).
