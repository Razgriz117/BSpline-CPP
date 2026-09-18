# ADR-0012: A4 Derivative-Discontinuity-Order Derivation — Verified Against Source

- **Status:** Accepted (implemented) — verified against source
- **Date:** 2026-09-04 (formalized during a post-implementation cleanup
  pass, then verified 2026-09-06 once the missing source PDF was located);
  the derivation itself dates to the Engineer A A4 work, already
  implemented and tested — this ADR originally deferred only its
  *literature cross-check*, which is now complete.

## Context

**Resolved 2026-09-06.** The blocking source PDF, `PHY5606_F25_Bsplines_v2.pdf`
(L. Argenti's PHY5606 B-splines course notes), has been located and added to
the repository at `docs/references/PHY5606_F25_Bsplines_v2.pdf` (alongside
`PHY5606_F25_ContinuumEigenstates.pdf`, `PHY5606_F25_Projects.pdf`, and the
Bachau et al. 2001 review — all previously cited by filename only in
`docs/SDD.md` [§1.5](../SDD.md#15-references) with no actual file present).

The course notes state the general B-spline continuity rule directly (p.1):

> "If a knot is distinct from its nearest neighbors, then the splines are at
> least $C^{k-2}$ there. If $\nu$ consecutive knots coincide, the regularity
> of the point decreases, and the splines are at least $C^{k-\nu-1}$ there."

This is the exact rule `strategicKnotsFromJoins`'s own internal re-derivation
already assumed (`docs/planning/engineer-a-plan-A4.md:381`: "a knot of total
multiplicity `m` gives continuity `C^{k-1-m}`" — identical formula, $m=\nu$).
Working the two production cases through it line-for-line:

- **Step** (V itself discontinuous): physically requires $\phi,\phi'$
  continuous but allows $\phi''$ to jump directly (since
  $\phi''=2(V-E)\phi$ inherits V's jump) — i.e. exactly $C^1$. Solving
  $k-\nu-1=1$ gives $\nu=k-2$; the base grid already contributes one
  occurrence of the knot, so extra multiplicity is $\nu-1=k-3=$
  `order-3` — **matches `TISE/tise.cpp`'s
  `case JoinType::Step: extra = std::max(0, order - 3)` exactly.**
- **StitchedKink** (V continuous, V' discontinuous): $\phi''$ stays
  continuous but $\phi'''$ jumps (inherits V's jump via
  $\frac{d}{dx}[2(V-E)\phi]$) — exactly $C^2$. Solving $k-\nu-1=2$ gives
  $\nu=k-3$, so extra multiplicity is $k-4=$ `order-4` — **matches
  `case JoinType::StitchedKink: extra = std::max(0, order - 4)` exactly.**

Both formulas used in production are now independently verified against the
actual cited source, not just the internal re-derivation and self-consistency
cross-check described below (which are both still valid and now corroborated,
not superseded).

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

**Citation chain traced during the original cleanup pass:** the required-treatment
table this formula feeds (Step→order-3, StitchedKink→order-4) is stated, in
the plan documents, to originate from `docs/planning/architecture-06-20.md`'s
"Collocation scheme" stakeholder-feedback section — but that section is a
design-meeting note recording a stakeholder decision about *what treatment
each join type gets*, not a literature derivation of *why* a given knot
multiplicity produces a given continuity order. The underlying B-spline
theory (knot multiplicity ↔ continuity order — a standard, well-established
result, not something specific to this project) is exactly what
`PHY5606_F25_Bsplines_v2.pdf` and de Boor's textbook (both already listed in
`docs/SDD.md` §1.5) are cited *for* — and, per the Context section above, the
course notes now confirm the formula line-for-line.

The derivation is implemented (`strategicKnotsFromJoins`) and tested
(hand-traced against all four canonical `detectPotentialStructure` cases,
`docs/planning/engineer-a-plan-A4.md:383`–388, plus the GTest cases the same
document adds), and is now additionally confirmed against its literature
source.

## Decision

The A4 discontinuity-order derivation, as implemented in
`strategicKnotsFromJoins`, is verified correct against `PHY5606_F25_Bsplines_v2.pdf`
(course notes' own $C^{k-\nu-1}$ continuity rule, worked through for both the
Step and StitchedKink cases in the Context section above) — no discrepancy
was found. The derivation continues to be used as-is; no code change results
from this verification.

## Consequences

- The knot-multiplicity formula used in production (`strategicKnotsFromJoins`)
  now rests on three independent confirmations: the internal re-derivation,
  the self-consistency check (the delta-potential cross-check), and this
  literature verification against the course notes it was always meant to
  match. No gap remains.
- No functional change results from this ADR — `strategicKnotsFromJoins`'s
  behavior, and all tests that currently pass against it, are unaffected.
- The companion planning doc, `docs/planning/a4-discontinuity-order-verification.md`,
  is updated to record the completed cross-check (its fallback numerical-
  verification approach — a convergence-order test against known analytic
  potentials per join type — remains available but is no longer needed to
  close this gap).
- **No further revisit trigger** — closed. (Independently, if a bug is ever
  observed in strategic-node placement's continuity behavior in a real run,
  that would still warrant revisiting this derivation, as with any
  production code.)

## Source

`docs/references/PHY5606_F25_Bsplines_v2.pdf` p.1 (the course notes' own
knot-multiplicity/continuity-order rule, $C^{k-\nu-1}$, that this ADR
verifies the implementation against); `docs/planning/engineer-a-plan-A4.md:381`
(the internal derivation and its self-consistency cross-check), `:383`–388
(hand-traced test verification); `docs/planning/architecture-06-20.md`
"Collocation scheme" stakeholder section (source of the required-treatment
table this formula implements, not of the multiplicity↔continuity theory
itself); `docs/SDD.md` [§1.5](../SDD.md#15-references) (external references:
`PHY5606_F25_Bsplines_v2.pdf`, C. de Boor, *A Practical Guide to Splines*).
See `docs/planning/a4-discontinuity-order-verification.md` for the completed
cross-check record and its (no-longer-needed) fallback numerical-verification
plan.
