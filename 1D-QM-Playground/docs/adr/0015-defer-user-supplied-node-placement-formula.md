# ADR-0015: Defer User-Supplied Arbitrary Node-Placement Formula

- **Status:** Accepted (deferred)
- **Date:** 2026-09-11 (formalized during a node-placement capability review
  against the original design docs, prompted by
  `docs/planning/engineer-a-plan-cleanup.md` item 8, "Verify node placement
  implementation is at least as capable as what we discussed in original
  design meetings"; the request itself dates to the 2026-07-03 stakeholder
  feedback pass and has been open, unformalized, ever since)

## Context

`docs/planning/architecture-06-20.md`'s "Stakeholder feedback (2026-07-03)"
section opens with two sentences before the strategic-placement/density
split that the rest of the section (and REQ-F-050) develops in detail:

> "The program should select a node distribution heuristically based on
> the potential, but leave open the option for the user to supply an
> arbitrary sequence of points via a formula $x(n)$."

The second half of that sentence — a user-supplied placement formula — is
never mentioned again. `docs/planning/tise-task-breakdown.md`'s A4 task
spec (REQ-F-050) carries forward only the first half (the strategic
placement table) and explicitly scopes *out* WKB-proportional density
("Do not implement a reference-energy-driven density scheme"); it does not
mention the arbitrary-formula option at all, in scope or out. `docs/SDD.md`
REQ-F-050's requirement text ("uniform spacing augmented by strategic
potential-driven nodes as the default collocation scheme") and
`docs/adr/0002-defer-wkb-collocation.md` both address only the
uniform-vs-WKB density question. No planning document, requirement, or ADR
carries the user-formula half of the stakeholder request forward, and none
formally defers it either — unlike WKB (ADR-0002) or CAP/outgoing-wave
boundary conditions (ADR-0011, which explicitly promoted a similarly
long-open SDD §12.B item to a decision), this request simply stops
appearing in any downstream document after the 2026-07-03 notes.

`config.yaml`'s `bspline` block (`config.yaml:17-21`, `docs/SDD.md`
[§6.1](../SDD.md#61-configuration-schema) Figure 4) confirms the current
state: it exposes only `n_nodes`, `order`, and `domain`. There is no field
through which a user could supply a placement formula, sequence, or even
select between grid strategies (e.g. `uniform` vs. a hypothetical
`custom`/`formula` option) — consistent with node placement being fully
automatic today (`docs/SDD.md:595`: "no config flag... this is determined
automatically from the potential, not user-set"), but that automatism was
only ever intended to cover the *heuristic* half of the stakeholder
sentence, not to imply the *override* half was rejected.

## Decision

Do not implement a user-supplied node-placement formula now. Node
placement remains fully automatic (uniform base + strategic knots per
REQ-F-050), with no config surface for a user-specified $x(n)$ sequence.
This ADR formalizes, as a conscious deferral, a stakeholder request that
has been open and unformalized since the 2026-07-03 design-feedback pass.

## Consequences

- Users who need a specific, hand-chosen node distribution the automatic
  heuristic does not produce (e.g. extra resolution around a feature the
  potential-structure detector cannot see, such as a region of rapid but
  continuous variation) have no way to request it through `config.yaml`
  today.
- Any future implementation would need to decide where a user sequence
  composes with strategic knots (REQ-F-050's automatic degenerate-knot
  insertion/B-spline removal) — e.g. whether user-supplied points are
  merged with strategic knots the way `buildStrategicRadialGrid` already
  merges strategic knots onto a uniform base, or whether supplying a
  custom sequence opts out of strategic placement entirely. This design
  question is itself unresolved and would need to be worked out at
  implementation time, not just the config-schema surface.
- **Revisit trigger:** if a stakeholder or user reasserts this request
  (e.g. for a teaching example needing a hand-tuned grid, or if the
  automatic heuristic is found insufficient for a specific potential in
  practice), revisit this decision. At that point, resolve the composition
  question above and extend `config.yaml`'s `bspline` block with an
  explicit field (e.g. an optional list of knot locations or a formula
  string parsed the same way `potential.function` is via `muparser`).

## Source

`docs/planning/architecture-06-20.md` "Stakeholder feedback (2026-07-03)"
(opening sentence, never revisited); `docs/planning/tise-task-breakdown.md`
A4 task spec (REQ-F-050, carries forward only the strategic-placement half);
`docs/SDD.md:205` (REQ-F-050 requirement text), `:595` (automatic,
no-config-flag placement); `config.yaml:17-21` (`bspline` block, no
formula/sequence field); `docs/adr/0002-defer-wkb-collocation.md` (the
sibling density-scheme deferral, which does not cover this request);
`docs/adr/0011-defer-cap-outgoing-wave-ecs-boundary-conditions.md` (the
precedent for formally promoting a long-open, unformalized design-meeting
item to an ADR); `docs/planning/engineer-a-plan-cleanup.md` item 8 (the
capability review that surfaced this gap).
