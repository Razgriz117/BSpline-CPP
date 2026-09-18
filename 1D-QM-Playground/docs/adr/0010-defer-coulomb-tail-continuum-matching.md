# ADR-0010: Defer Coulomb-Tail Continuum Matching (Figure 7 Case 2, Non-Flat Sub-Branch)

- **Status:** Superseded by [ADR-0013](0013-coulomb-tail-continuum-matching.md)
- **Date:** 2026-09-04 (formalized during a post-implementation cleanup pass;
  the gap itself was identified earlier, during the original TISE
  implementation phase — see Source); superseded 2026-09-06 once Coulomb-tail
  matching was actually implemented (ADR-0013).

## Context

REQ-F-030 requires the solver to determine boundary conditions automatically
from the potential's asymptote, including "Coulomb-sine-matched" continuum
normalization for a domain whose unbounded edge has a $-Z/r$ tail (SDD
[Figure 7](../SDD.md#523-internal-design) Case 2). `docs/SDD.md:562` states
this outright: the Coulomb-tail sub-branch of Case 2 "is required by
REQ-F-030 but not yet worked out at this level of detail by any source
document." The same gap is flagged independently in three other places —
`docs/planning/tise-task-breakdown.md`'s A1 note ("not yet derived by any
source document... do not attempt it as part of this task") and B5 note
("explicitly flagged as a follow-up, not part of this task"), and
`docs/planning/tise-release-readiness-plan.md`'s gap-audit table (listed
under "Explicitly deferred, no action").

What exists today: `tise::matchAsymptotic` (`TISE/tise.cpp:705`) implements
only the **flat**-asymptote matching formula — for each continuum energy $E$
it evaluates the confined-basis solution $\bar\psi_E(R)$ and
$\bar\psi_E'(R)$ at the box edge $R$, then extracts the amplitude and phase
shift by matching against a free-particle plane wave $\sin(kx+\delta)$:

$$
k = \sqrt{2E}, \qquad
A_E = \sqrt{\frac{2/\pi}{k\,\bar\psi_E(R)^2 + \bar\psi_E'(R)^2/k}}, \qquad
\delta(E) = \arctan\!\left(\frac{k\,\bar\psi_E(R)}{\bar\psi_E'(R)}\right) - kR
$$

This is the correct asymptotic form when $V(R)\to 0$ flat, or (per Case 3)
when the true asymptote is approximated as flat with a documented warning.
It is not the correct asymptotic form for a genuine $-Z/r$ Coulomb tail: the
long-range $1/r$ potential term never becomes negligible relative to the
kinetic term at any finite $R$, so the free-particle wave $\sin(kx+\delta)$
is not asymptotically similar to the true solution there — the true
solution instead approaches a **Coulomb wave function**, whose phase grows
logarithmically with $r$ rather than linearly (see the companion planning
doc for the derivation sketch). Matching a Coulomb-tailed problem against
`matchAsymptotic`'s flat-case formula (as `tise_solver_main.cpp`/`solveTISE`
would today, if asked to run continuum construction on a hydrogen-like
config with an unbounded domain) silently produces a wrong phase shift and
a wrong normalization constant — this is the same "runs to completion,
wrong answer" failure class ADR-0009 fixed for the grid/drop-set gap, not
yet fixed here.

## Decision

Do not implement Coulomb-tail continuum matching now. `matchAsymptotic`
continues to implement only the flat-asymptote formula; no Coulomb-specific
branch is added to `matchAsymptotic`, `buildContinuumState`, or the Case
1–3 asymptote classifier's dispatch logic. Any config that combines an
unbounded, genuinely-Coulomb-tailed domain (Case 2's non-flat sub-branch)
with continuum construction enabled continues to receive flat-asymptote
matching with no runtime warning — the same silent-approximation behavior
that exists today.

This decision only formalizes, as an ADR, a scope boundary that was already
set (informally, in planning docs) at the start of TISE implementation; it
does not change any code or existing behavior.

## Consequences

- Continuum-state normalization and phase shifts for a Coulomb-tailed
  unbounded domain remain numerically wrong (not merely approximate) if
  continuum construction is invoked on such a config. Bound-state
  diagonalization is unaffected — this gap is specific to the continuum
  branch of Case 2.
- Hydrogen bound-state runs (the project's current primary validation case,
  `H-BoundStates`) are unaffected: they either do not request continuum
  output, or the box is large enough that any resulting phase-shift error
  has not been part of the project's accuracy claims to date.
- The companion planning doc, `docs/planning/coulomb-tail-continuum-matching.md`,
  captures the derivation sketch (logarithmic Coulomb phase term vs. the
  flat case's linear $kR$ term) and the integration points
  (`matchAsymptotic`, `buildContinuumState`) so a future implementer does
  not have to re-derive the scope from scratch.
- **Revisit trigger:** if a config with a genuine unbounded Coulomb tail
  (not a Case-1 confining potential, and not a domain small enough that
  Case-3's flat approximation is defensible) needs continuum output —
  e.g. photoionization / scattering-phase-shift work that the current
  `tise_solver` production binary does not yet support for hydrogen-like
  targets — implement the Coulomb-matching branch following the companion
  planning doc and Bachau et al. (2001)'s treatment, and wire it into
  `matchAsymptotic`'s (or a new sibling function's) dispatch alongside the
  existing flat-case branch.

## Source

`docs/SDD.md:562` (REQ-F-030 gap statement, Figure 7 Case 2);
`docs/planning/tise-task-breakdown.md` (A1 and B5 follow-up notes, §6
"Explicitly Out of Scope"); `docs/planning/tise-release-readiness-plan.md`
(gap-audit table); `TISE/tise.cpp:705` (`matchAsymptotic`),
`TISE/tise.cpp:653` (`buildContinuumState`). See
`docs/planning/coulomb-tail-continuum-matching.md` for the derivation
sketch and integration plan.
