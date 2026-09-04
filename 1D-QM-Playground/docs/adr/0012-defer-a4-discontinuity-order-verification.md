# ADR-0012: Defer Formal Verification of the A4 Derivative-Discontinuity-Order Derivation

- **Status:** Accepted (deferred)
- **Date:** 2026-09-04 (formalized during a post-implementation cleanup pass;
  the derivation itself dates to the Engineer A A4 work, already
  implemented and tested — this ADR defers only its *literature
  cross-check*, not the derivation's use)

## ACTION NEEDED FROM A HUMAN

**This ADR cannot be closed by an agent.** Formal verification requires a
specific source PDF that does not currently exist anywhere in this
repository:

- **File needed:** `PHY5606_F25_Bsplines_v2.pdf` (L. Argenti's PHY5606
  B-splines course notes), already listed as an external reference in
  `docs/SDD.md` [§1.5](../SDD.md#15-references) — but, like every other
  external reference there, cited by filename only, with no actual file
  present.
- **Confirmed:** `find . -iname "*.pdf"` returns nothing anywhere in this
  repository — no PDF of any kind is checked in, not just this one.
- **Needed action:** a human locates this file and adds it to the repo.
  Suggested location: a new `docs/references/` directory (none exists yet;
  none of the other external references SDD §1.5 lists have one either, so
  creating it is itself a small new convention, not a fix for an existing
  broken path).
- Once added, the derivation in `docs/planning/engineer-a-plan-A4.md:381`
  (reproduced below) should be checked line-for-line against the course
  notes' treatment of knot multiplicity vs. spline continuity order, and
  optionally cross-checked against C. de Boor's *A Practical Guide to
  Splines* (also SDD §1.5, and considerably more available as a published
  textbook than the course notes).

This is flagged, not blocking: the derivation is already implemented and
covered by hand-traced test cases (see Context below) and is not gated on
this verification. See the companion planning doc,
`docs/planning/a4-discontinuity-order-verification.md`, for the fallback
numerical-verification approach if the literature cross-check turns out to
be inconclusive once the PDF is available.

## Context

`strategicKnotsFromJoins` (`TISE/tise.cpp`, derived in
`docs/planning/engineer-a-plan-A4.md`) converts a detected potential-join
type (`Step`, `StitchedKink`, `Singular`, `Continuous`) into an "extra knot
multiplicity" to insert at that join point, so that the resulting B-spline
basis has exactly the continuity the physics requires there — no more (which
would waste basis flexibility) and no less (which would fail to represent a
genuine kink/discontinuity in $\psi$ or its derivatives).

The formula, from `docs/planning/engineer-a-plan-A4.md:381`:

> For a B-spline basis of order `k`, a knot of *total* multiplicity `m`
> gives continuity `C^{k-1-m}`. To make ψ's derivative of order `n` the
> first discontinuous one, need continuity `C^{n-1}`, i.e. `m = k-n`. An
> ordinary grid point already has base multiplicity 1, so *extra*
> multiplicity is `m-1 = k-n-1`. Step forces a ψ'' jump (`n=2`, from
> `ψ'' = 2m(V-E)ψ`) → `extra = order-3`. StitchedKink forces a ψ''' jump
> (`n=3`, from differentiating once more) → `extra = order-4`.

The plan document itself calls this "the team's own re-derivation,
independent of the (missing-from-repo) source PDF," and includes a
self-consistency cross-check (a delta-function potential forcing a $\psi'$
jump, $n=1$, giving `extra = order-2` — matching the standard delta-matching
condition independently of the missing PDF).

**Citation chain traced during this cleanup pass:** the required-treatment
table this formula feeds (Step→order-3, StitchedKink→order-4) is stated, in
the plan documents, to originate from `docs/planning/architecture-06-20.md`'s
"Collocation scheme" stakeholder-feedback section — but that section is a
design-meeting note recording a stakeholder decision about *what treatment
each join type gets*, not a literature derivation of *why* a given knot
multiplicity produces a given continuity order. The underlying B-spline
theory (knot multiplicity ↔ continuity order — a standard, well-established
result, not something specific to this project) is exactly what
`PHY5606_F25_Bsplines_v2.pdf` and de Boor's textbook (both already listed in
`docs/SDD.md` §1.5) are cited *for* — but, per the ACTION NEEDED section
above, the PDF is not actually present in the repo to check against.

The derivation is already implemented (`strategicKnotsFromJoins`) and tested
(hand-traced against all four canonical `detectPotentialStructure` cases,
`docs/planning/engineer-a-plan-A4.md:383`–388, plus the GTest cases the same
document adds). This ADR is about the *formal verification* of the
derivation against its literature source, not about the derivation's
correctness in a working sense — the existing tests already exercise it as
implemented, they just cannot confirm it matches the course-notes' treatment
line-for-line, since that source is not available to check against.

## Decision

Do not attempt to formally verify the A4 discontinuity-order derivation
against `PHY5606_F25_Bsplines_v2.pdf` now — the file is not available in
this repository, and this ADR does not block on obtaining it (per the
ACTION NEEDED section, that is a human's task, not an agent's). The
derivation, as implemented in `strategicKnotsFromJoins` and validated by its
existing self-consistency cross-check and hand-traced test cases, continues
to be used as-is.

## Consequences

- The knot-multiplicity formula used in production (`strategicKnotsFromJoins`)
  rests on an internal re-derivation plus one self-consistency check (the
  delta-potential cross-check), not on a literature-verified derivation.
  This is a real, if likely low-risk, gap: the internal re-derivation is
  standard B-spline theory and the self-consistency check is a genuine
  independent test, but neither is the same as confirming against the
  specific course-notes treatment this project already committed to citing.
- No functional change results from this ADR — `strategicKnotsFromJoins`'s
  behavior, and all tests that currently pass against it, are unaffected.
- The companion planning doc,
  `docs/planning/a4-discontinuity-order-verification.md`, records what the
  derivation assumes, what a literature cross-check (once the PDF is
  available) should confirm, and a fallback numerical-verification approach
  (a convergence-order test against known analytic potentials per join
  type) usable even if the literature cross-check turns out inconclusive.
- **Revisit trigger:** once a human adds `PHY5606_F25_Bsplines_v2.pdf` to
  the repository (suggested location: `docs/references/`), perform the
  line-for-line cross-check described in the companion planning doc and
  update this ADR's status accordingly (either closing the gap, or, if a
  discrepancy is found, opening a follow-up correcting
  `strategicKnotsFromJoins`). Independently, if a bug is ever observed in
  strategic-node placement's continuity behavior (e.g. a Step or
  StitchedKink join not producing the expected discontinuity order in a
  real run), treat that as an urgent trigger to revisit this derivation
  regardless of whether the PDF has been added yet.

## Source

`docs/planning/engineer-a-plan-A4.md:381` (the derivation and its
self-consistency cross-check), `:383`–388 (hand-traced test verification);
`docs/planning/architecture-06-20.md` "Collocation scheme" stakeholder
section (source of the required-treatment table this formula implements,
not of the multiplicity↔continuity theory itself); `docs/SDD.md`
[§1.5](../SDD.md#15-references) (external references: `PHY5606_F25_Bsplines_v2.pdf`,
C. de Boor, *A Practical Guide to Splines*). See
`docs/planning/a4-discontinuity-order-verification.md` for the fallback
numerical-verification plan.
