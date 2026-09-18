# Complex Boundary Conditions: CAP vs. Outgoing-Wave (Siegert) vs. ECS

**Status:** Deferred — see [ADR-0011](../adr/0011-defer-cap-outgoing-wave-ecs-boundary-conditions.md).
This document is the companion planning doc for that ADR: a comparative
overview of the three options, written for a reader who knows the current
TISE solver (real-valued, `DSBGV`-based, Dirichlet/Case-1–3 boundary
conditions per REQ-F-030) but not necessarily complex-scaling/CAP theory.

## Problem statement

The current solver produces, for any potential, a finite set of real
eigenvalues from a finite box $[0, r_\text{max}]$ with a Dirichlet wall
(or a Case-2/Case-3 asymptote-matched continuum, per REQ-F-030) at the
outer edge. This is enough to characterize **bound states** exactly and
**genuine scattering states** approximately (as finite-box pseudostates,
per SDD [§5.2.3](../SDD.md#523-internal-design)). It cannot directly produce
two physically important quantities:

- **Resonance energies and widths.** A resonance is a quasi-bound state
  that eventually decays — it has a complex energy
  $E = E_r - i\Gamma/2$, where $\Gamma$ (the width) is inversely
  proportional to the resonance's lifetime. A real-valued Hermitian
  eigenproblem (which is what `DSBGV` solves) can only ever produce real
  eigenvalues — resonance widths are fundamentally invisible to the current
  pipeline, not merely approximated.
- **Ionization flux / outgoing probability current** through the outer
  boundary. With a hard Dirichlet wall, all probability that reaches
  $r_\text{max}$ reflects back rather than escaping — appropriate for a
  bound-state calculation, but wrong for tracking how much of a
  wavefunction has "ionized" past the box in a driven (TDSE) calculation.

Three standard techniques address one or both of these. They are not
mutually exclusive in general (CAP and ECS are sometimes combined), but for
a first adoption, one must be chosen to prototype first.

## The three options

### 1. Complex absorbing potential (CAP)

**Idea.** Add a smooth, purely imaginary term $-iW(r)$ to the potential,
turned on only near the outer wall (e.g. for $r > r_0$, ramping from 0 to
some strength $W_0$). The Hamiltonian becomes
$H \to H - iW(r)$. Physically, $W(r)$ acts like a "sponge": it damps the
part of the wavefunction that reaches the outer region before it can
reflect off the Dirichlet wall, so the box no longer artificially confines
outgoing flux.

**What it touches in this codebase.** The Hamiltonian matrix `H` that
`fillBandedMatrices` (`TISE/tise.cpp`) assembles would need a
new, energy-independent, purely imaginary diagonal-ish contribution
(banded, following the same $\langle B_i | W | B_j\rangle$ overlap-integral
structure the real potential term already uses). The resulting matrix
$H - iW$ is complex, but only complex-symmetric (not general complex or
Hermitian) if $W$ is real and the basis is real — this is the mildest
possible extension of the existing real-symmetric-banded machinery: it
still keeps the real-space B-spline basis and grid completely unchanged, only
the eigensolver changes (from `DSBGV`, real-symmetric, to a complex or
complex-symmetric banded/generalized eigensolver — LAPACK has routines for
this, e.g. `ZGBSV`-family or a complex-symmetric variant, though banded
complex-symmetric generalized eigensolvers are less standard than `DSBGV`
and may require either a dense fallback or a different LAPACK routine
family entirely).

**Trade-offs.** CAP introduces two free parameters (onset radius $r_0$ and
strength $W_0$, sometimes also a ramp shape/exponent) that must be tuned:
too weak and reflection still occurs; too strong and the CAP itself
reflects (a sharp imaginary step is just as bad as a hard wall). There is
no first-principles way to pick these without some empirical tuning per
problem, though standard recipes exist (e.g. a low-order polynomial or
Manolopoulos-style CAP with tabulated near-optimal parameters). It does not
directly give resonance widths as cleanly as Siegert/ECS — extracting
$\Gamma$ from a CAP calculation typically requires either fitting the decay
of a wavepacket's norm over time (a TDSE-side calculation, not a TISE
one) or using CAP inside a complex-eigenvalue formalism, which is a step
beyond the "just damp outgoing flux" description above.

### 2. Outgoing-wave (Siegert) boundary conditions

**Idea.** Instead of $\psi(r_\text{max}) = 0$, require $\psi$ to look like
a purely outgoing wave at the boundary: $\psi(r) \sim e^{ikr}$ for
$r \gtrsim r_\text{max}$, with $k$ complex when $E$ is complex. Solving the
Schrödinger equation with this boundary condition instead of Dirichlet
directly produces complex eigenvalues $E = E_r - i\Gamma/2$ for resonance
states — this is the cleanest, most direct route to resonance
positions/widths, with no tunable damping parameters the way CAP has.

**What it touches.** In practice, imposing a genuinely outgoing-wave
condition on a real grid is numerically delicate (the outgoing solution
grows exponentially for complex $k$ in the wrong direction unless the
contour is rotated — see ECS below), so Siegert boundary conditions are
rarely implemented directly on the real axis; they are almost always
implemented *via* exterior complex scaling, which is why option 3 below is
usually presented as "how you actually get Siegert states" rather than as
a fully separate alternative.

### 3. Exterior complex scaling (ECS)

**Idea.** Rotate the radial coordinate into the complex plane beyond some
radius $r_0$:

$$
r \;\to\; \tilde r(r) = \begin{cases} r & r \le r_0 \\ r_0 + (r - r_0)\,e^{i\theta} & r > r_0 \end{cases}
$$

for some rotation angle $\theta$. Under this transformation, a Siegert
(outgoing, exponentially divergent on the real axis) state becomes
exponentially *decaying* and square-integrable — the "complex rotation"
turns an awkward boundary condition into an ordinary $L^2$ eigenproblem,
just on a complex-deformed contour. The resulting Hamiltonian is
complex-symmetric (not Hermitian) with genuinely complex eigenvalues for
resonances; bound states remain on the real axis, unaffected by the
rotation (by Cauchy's theorem / analyticity arguments, assuming $\theta$
is chosen appropriately and the potential is analytic enough).

**What it touches.** This is the most invasive of the three options: the
radial grid itself becomes complex-valued beyond $r_0$, so every place the
current code assumes `Real` (`TISE/tise.hpp`'s scalar type alias) grid
points, B-spline knots, or matrix elements would need a complex
counterpart, or the whole numeric type would need to become
`std::complex<double>`-based beyond $r_0$. `bspline::BSpline`'s evaluation
and integration routines, `fillBandedMatrices`, and the eigensolver all sit
downstream of that choice.

**Trade-offs.** ECS is the most rigorous and general of the three (no
tunable damping parameters, works for resonances of any width, and is the
standard approach in the atomic/molecular physics literature Bachau et al.
2001 surveys), but it is also the largest implementation lift by a wide
margin — it is not a localized addition the way CAP is.

## Comparison table

| | Resonance widths | New free parameters | Codebase impact | Numerical maturity here |
|---|---|---|---|---|
| CAP | Indirect (needs extra step) | Onset radius, strength, ramp shape | Small: one added complex diagonal-ish term; still real B-spline basis/grid | Compatible with existing real-valued `DSBGV`-adjacent infrastructure at the basis level; needs a new complex eigensolver call |
| Outgoing-wave (Siegert) | Direct | None (in principle) | Not practical to implement directly without ECS | N/A (see ECS) |
| ECS | Direct | Rotation angle $\theta$, rotation radius $r_0$ | Large: complex-valued grid/basis/matrix elements throughout | Requires a substantially new complex-arithmetic code path |

## Recommendation

**Prototype CAP first**, as SDD `docs/SDD.md:1329`'s own "Initial analysis"
recommendation already concluded ("CAP support at the outer boundary is the
recommended next addition... since it does not require changes to the
real-valued eigensolver" — noting that in practice a complex eigensolver
*is* eventually needed once $-iW(r)$ makes $H$ complex, but the basis/grid
side of the codebase is untouched, which is the more invasive half of the
alternatives). CAP is the smallest change relative to the existing
codebase, has the lowest implementation risk, and is directly useful for
ionization-flux studies the TDSE roadmap already anticipates — even without
a full resonance-width extraction, a working CAP unblocks "how much
probability left the box" questions immediately. ECS should be considered
later, if and when resonance widths specifically (not just absorbing
outgoing flux) become a project goal, since it is the more rigorous and
more expensive option. Outgoing-wave/Siegert boundary conditions are best
understood as "what ECS gives you," not as an independently simpler
alternative to implement.

## Non-goals of this document

This document does not specify a CAP functional form, tune parameters, or
provide pseudocode for a `TISE/tise.cpp` implementation — it is a
comparative decision-support document, written so that whichever option is
eventually picked up (per ADR-0011's revisit trigger) starts from an
informed choice rather than from the SDD's brief paragraph-level survey.
