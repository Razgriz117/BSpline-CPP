# Coulomb-Tail Continuum Matching (Figure 7 Case 2, Non-Flat Sub-Branch)

**Status:** Implemented — see [ADR-0013](../adr/0013-coulomb-tail-continuum-matching.md)
(supersedes [ADR-0010](../adr/0010-defer-coulomb-tail-continuum-matching.md)).
This document's derivation sketch and integration plan were followed closely
during implementation; it is kept as the derivation reference (the "why"
behind `tise::evaluateCoulombFunctions`/`matchAsymptotic`'s Coulomb branch),
not as an open item.

## Problem statement

REQ-F-030 (SDD [§3](../SDD.md#3-requirements)) requires the solver to pick a
boundary-condition treatment automatically from the potential's behavior at
an unbounded domain edge. [Figure 7](../SDD.md#523-internal-design)'s Case 2
covers "known analytic asymptote," which SDD [§5.2.3](../SDD.md#523-internal-design)
splits into two sub-branches:

- **Flat asymptote** ($V(x) \to V_\infty$, a constant, as $x \to \infty$):
  implemented — this is what `tise::matchAsymptotic`
  (`TISE/tise.cpp:705`) actually does today.
- **Coulomb asymptote** ($V(x) \to -Z/x$ as $x \to \infty$, the physically
  important case for any hydrogen-like or ionic radial problem): **not
  implemented**. `matchAsymptotic` is unconditionally called with the
  flat-case formula regardless of which sub-branch actually applies; there
  is no dispatch on asymptote type inside the continuum-construction path.

This matters physically, not just as a code-organization gap. For a flat
asymptote, the radial Schrödinger equation reduces at large $x$ to
$-\tfrac12\psi'' \approx (E - V_\infty)\psi$, whose two independent
solutions are the free-particle plane waves $e^{\pm ikx}$ with
$k=\sqrt{2(E-V_\infty)}$ — a real, constant wavenumber, giving the familiar
$\sin(kx+\delta)$ asymptotic form with a **linear** phase $kx$. For a
Coulomb tail, the potential term $-Z/x$ never becomes negligible relative to
the kinetic term no matter how large $x$ is (it falls off at the same
$1/x^2$ rate, roughly, as the next correction to the plane wave) — the
equation never truly decouples into the free-particle form. The correct
asymptotic solutions are the **Coulomb wave functions** $F_l(\eta, kx)$ and
$G_l(\eta, kx)$ (regular and irregular, respectively), whose phase grows
like

$$
kx - \eta\ln(2kx) - \frac{l\pi}{2} + \sigma_l(\eta)
$$

where $\eta = -Z/k$ is the Sommerfeld parameter and
$\sigma_l(\eta) = \arg\Gamma(l+1+i\eta)$ is the Coulomb phase shift. The
$-\eta\ln(2kx)$ term is the qualitative difference from the flat case: it is
a **logarithmic** correction to the phase that never vanishes as
$x\to\infty$ (unlike, say, a short-range potential's correction, which does
vanish). Matching against $\sin(kx+\delta)$ instead of the true
$F_l/G_l$-combination silently discards this term, producing a phase shift
and normalization constant that are wrong by an $x$-dependent (i.e.,
box-size-dependent) amount — not a small, controlled approximation error
the way Case 3's flat-approximation-with-warning is.

## Where this plugs into the existing code

Two functions in `TISE/tise.cpp` would need a Coulomb-aware branch:

- **`buildContinuumState`** (`TISE/tise.cpp:653`) constructs, for each
  energy on a continuum energy grid, the confined-basis coefficients of the
  scattering state $\bar\psi_E$ via the boundary-coupling formula
  (`states[E_idx][i] = (coeffs1[i] - E*coeffs2[i]) / (E - eigen.values[i])`).
  This part is asymptote-agnostic — it only depends on the confined
  eigenbasis and the Hamiltonian/overlap matrices, not on what lies beyond
  $R$. **No change needed here.**
- **`matchAsymptotic`** (`TISE/tise.cpp:705`) is where the asymptote-specific
  physics lives: it evaluates $\bar\psi_E(R)$, $\bar\psi_E'(R)$ via
  `bs.eval`, then matches them against the assumed asymptotic form to
  extract `A_E` (normalization) and `delta` (phase shift). **This is the
  function that needs a Coulomb-tail branch.** The natural shape of the
  change is a boolean/enum parameter (or a wrapping caller that dispatches
  based on the Case 1/2/3 classification A1's asymptote classifier already
  computes) that selects between the existing flat-case formula and a new
  Coulomb-case formula.

## Sketch of the Coulomb-matching formula

Following Bachau et al. (2001, *Rep. Prog. Phys.* **64**, 1815, already
cited in `docs/SDD.md` [§1.5](../SDD.md#15-references) as the project's
primary numerical-methods reference) and the standard scattering-theory
treatment: instead of writing the exterior solution as
$A\sin(kx+\delta)$, write it as a linear combination of the regular and
irregular Coulomb functions,

$$
\bar\psi_E(x) \;\underset{x\to\infty}{\longrightarrow}\; A_E\bigl[\cos\delta(E)\, F_l(\eta, kx) + \sin\delta(E)\, G_l(\eta, kx)\bigr]
$$

and match both the value and derivative of the confined-basis solution at
$x=R$ against this form and its derivative (computed from the standard
recurrence/asymptotic-expansion relations for $F_l$, $G_l$ — these are
well-tabulated special functions with existing numerical libraries, e.g.
GSL's `gsl_sf_coulomb_wave_FG_e`, or Boost.Math if a suitable routine
exists there). The two matching conditions
($\bar\psi_E(R) = \dots$, $\bar\psi_E'(R) = \dots$) give two equations for
the two unknowns $A_E$ and $\delta(E)$, exactly analogous to how the
existing flat-case formula solves for `A_E`/`delta` from
$\bar\psi_E(R)$/$\bar\psi_E'(R)$ — the derivation strategy carries over
directly, only the target functions change from $\sin$/$\cos$ to $F_l$/$G_l$.

Concretely, the two building blocks that must be supplied are:

1. A numerically stable evaluation of $F_l(\eta, kR)$, $G_l(\eta, kR)$ and
   their derivatives at the matching radius $R$, for the $\eta$, $k$, $l$
   values in play (this is the part requiring an external special-function
   routine — B-spline evaluation itself does not change).
2. The 2×2 linear (or trigonometric) solve for $A_E$, $\delta(E)$ from the
   value/derivative match — structurally the same kind of step
   `matchAsymptotic`'s existing `atan`/`sqrt` formulas already perform for
   the flat case, just with $F_l$/$G_l$ (and their derivatives) in place of
   $\sin(kR)$/$\cos(kR)$.

## Validation approach, once implemented

- Hydrogen ($Z=1$, $l=0$) continuum states at low energy have tabulated
  Coulomb phase shifts (trivially zero for a pure Coulomb potential with no
  additional short-range term, since $\delta(E)$ measures the *deviation*
  from pure Coulomb scattering — a pure $-1/r$ potential run through this
  matching should recover $\delta \approx 0$ to numerical precision, a
  strong and cheap sanity check).
- A hydrogen-plus-short-range-correction potential (e.g. a screened or
  truncated Coulomb term added inside $R$) gives a nonzero, independently
  computable phase shift to compare against, mirroring how the flat-case
  formula was presumably validated against known square-well/free-particle
  phase shifts.
- Cross-check against Bachau et al.'s own worked examples where available.

## Non-goals of this document

This document does not propose an implementation timeline or claim the
sketch above is a complete, ready-to-code specification — it is deliberately
a "what shape does this take" note, written so that whoever eventually picks
this up (per ADR-0010's revisit trigger) starts from the physics and
integration points already scoped, rather than from zero.
