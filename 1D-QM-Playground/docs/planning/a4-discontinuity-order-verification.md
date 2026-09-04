# A4 Derivative-Discontinuity-Order Derivation: Verification Plan

**Status:** Deferred — see [ADR-0012](../adr/0012-defer-a4-discontinuity-order-verification.md).
This document is the companion planning doc for that ADR.

## The one thing to know before reading further

**Verifying this derivation against the literature requires
`PHY5606_F25_Bsplines_v2.pdf` (L. Argenti's PHY5606 B-splines course notes),
and that file does not currently exist anywhere in this repository.**
`find . -iname "*.pdf"` returns nothing — no PDF of any kind is checked in.
A human needs to locate this file and add it to the repo (suggested
location: a new `docs/references/` directory) before the line-for-line
cross-check described below can actually be performed. This document
describes the verification plan for when that happens, plus a
numerical fallback that does not require the PDF at all — it does not
perform the verification itself.

## What the derivation claims

`strategicKnotsFromJoins` (`TISE/tise.cpp`, derived in
`docs/planning/engineer-a-plan-A4.md:381`) needs to decide, for each detected
join point in a piecewise potential, how many *extra* B-spline knot copies
to insert there so the resulting basis has exactly the right continuity —
enough to represent a genuine kink or discontinuity, not so much that basis
flexibility is wasted where the wavefunction is actually smooth.

The claimed relationship, in the derivation's own words:

> For a B-spline basis of order `k`, a knot of *total* multiplicity `m`
> gives continuity `C^{k-1-m}`.

This is a standard fact from B-spline theory (independent of anything
quantum-mechanical): each unit of knot multiplicity beyond 1 reduces the
number of continuous derivatives at that knot by one, starting from a
maximally-smooth simple knot's $C^{k-2}$ continuity for an order-`k`
(degree-$(k-1)$) B-spline basis. The derivation then applies this to the
physics: the Schrödinger equation $\psi'' = 2m(V-E)\psi$ ties $\psi$'s
second derivative directly to $V$, so wherever $V$ itself is discontinuous
(a **Step** join), $\psi''$ inherits that discontinuity, and wherever $V$'s
*first derivative* is discontinuous but $V$ itself is continuous (a
**StitchedKink** join — two smoothly-varying pieces meeting at a corner,
like $|x|$), $\psi'''$ (one derivative further, since Schrödinger relates
$\psi''$ to $V$, not $\psi'''$ to $V$) inherits it instead. Working backward
from "which derivative order `n` should be the first discontinuous one" to
"how much extra knot multiplicity does that require" gives:

| Join type | Discontinuous derivative order `n` | Continuity required, `C^{n-1}` | Extra multiplicity, `order - n - 1` |
|---|---|---|---|
| Step (V jumps) | 2 (ψ'') | C¹ | `order - 3` |
| StitchedKink (V' jumps, V continuous) | 3 (ψ''') | C² | `order - 4` |
| Singular / Continuous | — (no join-induced discontinuity) | — | 0 |

The derivation includes one independent self-consistency check: a
delta-function potential term forces a $\psi'$ jump ($n=1$), giving
`extra = order - 2` — exactly one less than the Step case, and matching the
standard, well-known delta-function matching condition
($\psi'(x_0^+) - \psi'(x_0^-) = 2m\cdot(\text{delta strength})\cdot\psi(x_0)$,
a boundary-matching result independent of B-splines entirely). This
cross-check is genuine evidence the multiplicity↔continuity mapping is
applied correctly, but it is not the same as confirming the base B-spline
fact itself (`m` extra multiplicity ⟹ `C^{k-1-m}` continuity) against the
course notes' own treatment, which may use different indexing conventions,
define "order" differently (order vs. degree is a common source of
off-by-one errors across different B-spline references), or otherwise phrase
the result in a way worth double-checking explicitly.

## What the PDF cross-check should confirm, once available

1. **Order vs. degree convention.** Confirm `PHY5606_F25_Bsplines_v2.pdf`
   uses the same "order `k`" convention as this codebase (order = degree +
   1; a cubic B-spline is order 4) — this is the single most common source
   of off-by-one bugs when porting a formula from one B-spline reference to
   another. `docs/planning/bsplines.md` (already in the repo, cited in SDD
   §1.5) should also be checked for which convention *it* uses internally,
   as a second cross-reference independent of the missing PDF.
2. **The base multiplicity↔continuity relationship itself**
   (`m` multiplicity ⟹ `C^{k-1-m}` continuity at a simple B-spline knot),
   confirmed against whatever form the course notes state it in.
3. **Any project-specific caveats** the course notes might raise for the
   radial/physical setting this project uses (e.g., behavior of the
   relationship right at a domain edge or at $r=0$, which is already
   handled specially elsewhere in this codebase via the classic `{1}`
   drop-set convention and doesn't obviously interact with this derivation,
   but is worth explicitly ruling out during the cross-check rather than
   assuming).
4. **De Boor's textbook** (*A Practical Guide to Splines*, also cited in SDD
   §1.5) as a second, independently-available cross-check — unlike the
   course notes, this is a published book and may be easier to obtain in
   the meantime if only a partial or preliminary check is wanted before the
   PDF itself is located.

## Fallback: numerical verification (usable even without the PDF)

If the literature cross-check is inconclusive, or as an independent check
regardless of the literature outcome, a **convergence-order test** can
verify the derivation's practical consequences directly, without needing
any external reference document:

1. Construct a potential with a known, exactly-analytic solution across a
   deliberately-placed Step join (e.g. a finite square well/barrier, whose
   piecewise-exponential/piecewise-sinusoidal solution is known in closed
   form on each side, with a known, computable derivative-continuity
   structure at the join).
2. Run `strategicKnotsFromJoins`'s current formula (extra = order - 3 for
   Step) and confirm the B-spline solution's $\psi''$ at the join point
   converges to the analytic solution's $\psi''$ discontinuity at the
   correct rate as basis order/grid density increase, while $\psi$ and
   $\psi'$ remain continuous (as they should).
3. Repeat with the extra multiplicity deliberately set one *too low* (i.e.
   testing what happens if the formula were `order - 4` instead of
   `order - 3` for a Step join): confirm this produces a *detectably worse*
   fit at the join — i.e. that the current formula is not just "a value
   that works" but "the value the physics actually requires," by showing a
   plausible off-by-one alternative measurably underperforms.
4. Repeat the same two checks for a StitchedKink join, using a potential
   with a known corner (e.g. a piecewise-linear potential, like the
   `10-x`/`x` example already used in `engineer-a-plan-A4.md`'s own
   hand-traced test cases), checking $\psi'''$'s discontinuity instead of
   $\psi''$'s.
5. Report convergence order (how the join-point derivative error scales
   with basis order or grid refinement) rather than just a pass/fail
   threshold, since that is what actually distinguishes "the right
   multiplicity" from "a multiplicity that happens to look adequate at one
   resolution."

This numerical approach cannot substitute for confirming the underlying
B-spline theory statement itself is quoted correctly (that specifically
needs the literature source), but it can independently confirm the
*physics application* of the derivation is correct — which, combined with
the already-existing hand-traced test cases
(`docs/planning/engineer-a-plan-A4.md:383`–388), is meaningful evidence on
its own even before the PDF is located.

## Non-goals of this document

This document does not perform either the literature cross-check (blocked
on the PDF) or the numerical convergence-order test (not attempted here) —
it records the plan for both so that whoever picks this up next (a human
adding the PDF, or an implementer running the convergence test) has a
concrete starting point rather than an open-ended "go verify this."
