# ADR-0006: Defer Crank-Nicolson-Family TDSE Propagators

- **Status:** Accepted (deferred)
- **Date:** 2026-08-21

## Context

SDD §5.3.3 Figure 8 explicitly leaves the driven-case ($H(t)=H_0+H_\text{int}(t)$) propagator undecided: *"numerically integrate forward in time (propagator scheme is an implementation detail, not user-configurable)."* The config-schema-design doc confirms this is deliberate — an earlier `tdse.chebyshev_order` config field was removed for exactly this reason. No prior ADR covers it.

Four candidate numerical schemes were considered for advancing the driven-case state vector.

Throughout, $H_\text{int}(t)=-\hat O\,\mathcal E(t)$, with $\hat O=\hat x$ (length gauge) or $\hat O=\hat p_x$ (velocity gauge) per `tdse.operator.gauge`, and $\mathcal E(t)$ the field envelope (`tdse.field.expression`). $C$ is the TISE eigenvector matrix (`tise::EigenResult::vectors`), satisfying $C^T\mathbf SC=\mathbb 1$ and $C^TH_0^\text{BS}C=\mathrm{diag}(E_n)$ by construction of the generalized eigenproblem. $\hat O_\text{eig}\equiv C^T\hat O_\text{BS}C$ is the coupling operator rotated into the eigenbasis (`docs/planning/bsplines.md`, "Eigenstates as Linear Transformations").

Four options were compared.

### Option 1a — Dense Crank-Nicolson (Cayley), eigenbasis

$$\left(\mathbb{1} + i\theta H_\text{eig}(t_{n+1/2})\right)\boldsymbol\alpha^{n+1} = \left(\mathbb{1} - i\theta H_\text{eig}(t_{n+1/2})\right)\boldsymbol\alpha^{n}, \qquad \theta=\frac{\Delta t}{2\hbar},\quad H_\text{eig}(t)=\mathrm{diag}(E_n)+H_{\text{int,eig}}(t)$$

Exactly unitary for any $\Delta t$: the Cayley transform of a Hermitian generator always has unit-modulus eigenvalues. $H_{\text{int,eig}}(t)=-\hat O_\text{eig}\mathcal E(t)$ is dense (rotating a banded B-spline-space operator through $C$ densifies it in general), so each step is a dense complex linear solve, and since $\mathcal E(t)$ changes every step, the Cayley matrix must be re-factored every step — no reuse between steps.

### Option 1b — Generalized (banded) Crank-Nicolson, B-spline basis

Projecting $i\hbar\,\partial_t\psi=\hat H(t)\psi$ directly onto $\{B_j\}$ (not the eigenbasis) gives the mass-matrix ODE $i\hbar\mathbf S\dot{\mathbf d}(t)=H_\text{BS}(t)\mathbf d(t)$ — the time-domain analog of $\mathbf{Hc}=E\mathbf{Sc}$. Trapezoidal discretization gives the generalized Cayley step:

$$\left(\mathbf S+i\theta H_\text{BS}(t_{n+1/2})\right)\mathbf d^{n+1}=\left(\mathbf S-i\theta H_\text{BS}(t_{n+1/2})\right)\mathbf d^{n}$$

**This is exactly equivalent to Option 1a — proven two independent ways:**

*Proof (i), Löwdin substitution.* $\mathbf S$ is real symmetric positive-definite (a Gram matrix of independent functions), so a unique $\mathbf S^{1/2}$ exists. Substituting $\mathbf d=\mathbf S^{-1/2}\tilde{\mathbf d}$ and left-multiplying by $\mathbf S^{-1/2}$ collapses the recursion to plain identity-form CN, $(\mathbb 1+i\theta\tilde H)\tilde{\mathbf d}^{n+1}=(\mathbb 1-i\theta\tilde H)\tilde{\mathbf d}^{n}$ with $\tilde H=\mathbf S^{-1/2}H_\text{BS}\mathbf S^{-1/2}$ Hermitian — exactly unitary, and $\|\tilde{\mathbf d}\|_2^2=\mathbf d^\dagger\mathbf S\mathbf d=\int|\psi|^2dx$ is exactly the physical norm. So 1b is exactly unitary w.r.t. the physical norm, without ever forming $\mathbf S^{\pm1/2}$ explicitly.

*Proof (ii), direct substitution of $C$ (stronger — shows identical trajectories, not just "both unitary").* Substituting $\mathbf d=C\boldsymbol\alpha$ and left-multiplying by $C^T$, using $C^T\mathbf SC=\mathbb1$:

$$\left(C^T\mathbf SC+i\theta\,C^TH_\text{BS}C\right)\boldsymbol\alpha^{n+1}=\left(C^T\mathbf SC-i\theta\,C^TH_\text{BS}C\right)\boldsymbol\alpha^{n} \;\Longrightarrow\; (\mathbb 1+i\theta H_\text{eig})\boldsymbol\alpha^{n+1}=(\mathbb 1-i\theta H_\text{eig})\boldsymbol\alpha^{n}$$

— **Option 1a's exact recursion.** Since $C$ is fixed and $\boldsymbol\alpha^0=C^{-1}\mathbf d^0$ for matching initial conditions, induction gives $\boldsymbol\alpha^n_{(1a)}=\boldsymbol\alpha^n_{(1b)}$ **exactly**, step for step — not merely two separately-unitary schemes, the identical trajectory. (Useful later: this is a free, strong cross-validation test if 1b is ever implemented — transform every 1b snapshot through $C^{-1}$ and assert near-machine-precision agreement with 1a.)

$\hat O_\text{BS}$ (length: $x_{ij}=\int B_ixB_jdx$; velocity: $p_{ij}=-i\hbar\int B_iB_j'dx$) vanishes for $|i-j|\ge k$ by the same compact-support argument that makes $\mathbf H_0^\text{BS}$/$\mathbf S$ banded (`TISE/BSpline.hpp`'s `integral()` — multiplication by $x$ or differentiation doesn't enlarge a B-spline's knot support). So $H_\text{BS}(t)$ stays banded, bandwidth $k-1$, for all $t$, both gauges — 1b is a **banded** complex linear solve, not dense. But $(\mathbf S+i\theta H_\text{BS})$ is banded and non-Hermitian ($\mathbf S$ real-symmetric plus $i\theta\times$Hermitian), so it needs a *general*-banded solver (LAPACK `zgbsv`-family, general-banded `AB` storage with pivot-fill rows) — a new binding, mirroring how `tise.cpp` already declares `dsbgv_` via `extern "C"`, but a different storage layout than `dsbgv`'s symmetric-banded one, and non-trivial new code, not a copy-paste.

### Option 2 — Symmetric split-operator (Strang splitting), eigenbasis

$$|\psi(t+\Delta t)\rangle = e^{-i\hat H_0\Delta t/(2\hbar)}\;e^{+i\,\Delta t\,\mathcal E(t+\Delta t/2)\,\hat O/\hbar}\;e^{-i\hat H_0\Delta t/(2\hbar)}\;|\psi(t)\rangle$$

Sign check: $i\hbar\,d|\psi\rangle/dt=H(t)|\psi\rangle \Rightarrow e^{-iH_\text{int}(t_{n+1/2})\Delta t/\hbar}=e^{-i(-\hat O\mathcal E(t_{n+1/2}))\Delta t/\hbar}=e^{+i\hat O\mathcal E(t_{n+1/2})\Delta t/\hbar}$ — matches. The two outer factors, in the eigenbasis, are exactly `tevol::timeEvolveState` called with $t=\Delta t/2$ — reused verbatim, zero new numerics, already covered by ~15 existing passing tests. $\hat O_\text{eig}$ is Hermitian in both gauges (length: real-symmetric, $x_\text{eig}=C^Tx_\text{BS}C$; velocity: $p_\text{eig}=-i\hbar(C^TKC)$ with $K$ explicitly antisymmetrized, and $i\times$real-antisymmetric is Hermitian) and **time-independent** — only the scalar $\mathcal E(t)$ varies. Diagonalizing $\hat O_\text{eig}=W\Lambda W^\dagger$ **once**, up front, via `Eigen::SelfAdjointEigenSolver` (already reachable through `time_evolution.hpp`'s existing `#include <Eigen/Dense>` — `Eigen/Eigenvalues` is pulled in transitively, zero new dependency), makes every subsequent step a cheap diagonal-phase update:

$$e^{i\theta\hat O_\text{eig}}v = W\,\mathrm{diag}\!\left(e^{i\theta\lambda_k}\right)\,W^\dagger v, \qquad \theta=\frac{\Delta t\,\mathcal E(t_{n+1/2})}{\hbar}$$

with no per-step factorization at all. Exactly unitary: a product of three exponentials of Hermitian generators is exactly unitary, for any $\Delta t$.

### Option 3 — Split-operator + Crank-Nicolson hybrid

Same sandwich structure as Option 2, but the interaction step uses the Cayley/CN approximation to $e^{i\theta\hat O_\text{eig}}$ instead of the exact exponential. **Considered and rejected.** The Cayley system $\mathbb 1+i\theta\hat O_\text{eig}/2$ changes every step (θ varies), so a naive implementation pays a fresh dense factorization every step for no offsetting benefit. Even a "smart" implementation — noticing $\hat O_\text{eig}$ itself never changes and diagonalizing it once anyway — collapses to $W\,\mathrm{diag}\!\left(\frac{1-i\theta\lambda_k/2}{1+i\theta\lambda_k/2}\right)\,W^\dagger$, computed from the *same* precomputed $\lambda_k$ Option 2 already needs, at the same cost, but now only a Padé[1/1] **approximation** to $e^{i\theta\lambda_k}$ instead of the exact value — accuracy given up for nothing. No regime was found where this wins over Option 2; its only coherent motivation is reusing a generic CN primitive that would otherwise need to exist for some unrelated reason.

### Why unitarity and accuracy order don't discriminate between these

All of 1a/1b/2 are **exactly** unitary for any $\Delta t$ — the Cayley transform of a Hermitian generator and the exponential of $i\times$Hermitian are both exactly unitary by construction. This had been the previous draft's stated justification for Crank-Nicolson; it doesn't survive as a discriminator once split-operator is in the comparison. Both families are standard $O(\Delta t^2)$-global-error schemes — a tie in asymptotic order (which has the smaller error *constant* for this specific $H_0,\hat O$ pair is an empirical question the Rabi-flopping benchmark in `tdse-task-breakdown.md` §5 would answer, not asserted here).

### What does discriminate

1. **Exact vs. approximate reduction to the already-built field-free path.** Under CN (1a/1b), zero-field DRIVEN matches the existing `STATIC` path only to $O(\Delta t^2)$. Under Option 2, zero-field DRIVEN reduces to $e^{-i\hat H_0\Delta t/\hbar}$ by the exponential group law — **bit-for-bit identical** (mod floating point) to what `timeEvolveState` already computes. A strictly stronger correctness property, not a cosmetic one.
2. **Cost, concretely, at this project's actual scale** ($n_\text{En}\approx59$ per `bsplines.md`'s stated parameters; `dt=0.3`, `t_final=300` → 1000 steps, per `config.yaml`):

   | Option | Per-step cost | × 1000 steps |
   |---|---|---|
   | 1a — dense CN, fresh factorization every step | $O(n^3)\approx2.05\times10^5$ | $\approx2\times10^8$ |
   | 1b — banded CN, fresh banded factorization every step | $O(nk^2)\approx8.5\times10^3$ | $\approx8.5\times10^6$ (~24× fewer than 1a) |
   | 2 — precompute $W,\Lambda$ once; per-step $O(n^2)$ mat-vecs + $O(n)$ phases | $\approx3.5\times10^3$ | $\approx3.5\times10^6$ (~59× fewer than 1a) |

   All three totals are sub-millisecond to a few milliseconds in absolute terms at this scale — not an urgent problem today. But Option 2 is strictly cheapest with **no accuracy tradeoff**, so there is no cost to taking the win now, and the margin only grows if basis size grows later.
3. **Dependency footprint.** Option 2 needs nothing beyond Eigen's already-linked core. Option 1b needs a new general-banded LAPACK binding (`zgbsv`-family, non-trivial new storage-layout code). Option 1a needs nothing new but pays 1a's cost.

## Decision

Use the symmetric split-operator / Strang-splitting propagator (Option 2 above — eigenbasis, coupling operator diagonalized once via `Eigen::SelfAdjointEigenSolver`) as the driven-case propagator. Crank-Nicolson-family propagators (Options 1a and 1b) are not implemented. The split-operator/CN hybrid (Option 3) is not implemented either.

## Consequences

- **`docs/planning/tdse-task-breakdown.md`'s Engineer A tasks build the split-operator method directly**, not a generic Crank-Nicolson primitive: precompute $\hat O_\text{eig}=W\Lambda W^\dagger$ once; apply the sandwich (`timeEvolveState` half-step → exact interaction kick via $W,\Lambda,\theta$ → `timeEvolveState` half-step) per snapshot interval. The `+i` sign in the interaction kick (not `-i`) is the propagator-side responsibility this decision assigns — get it wrong and the field pushes population the wrong way with no crash to flag it, so it is called out explicitly here and in the task breakdown.
- **The field-free `STATIC` path (`tevol::timeEvolveState`, already implemented and tested) is kept as an explicit fast path**, not replaced. It remains useful as a cheap oracle: DRIVEN-with-zero-field must now match it *exactly* (not approximately) per the split-operator method's exponential group-law property — a stronger regression check than the CN family could have offered. Collapsing the fast path into "just call the DRIVEN loop with a zero field" is a legitimate future simplification specifically *because* the equivalence is exact here; it would have been a silent correctness regression under Options 1a/1b, where the match is only $O(\Delta t^2)$. This asymmetry should not be forgotten if the code is later refactored by someone unaware of which propagator is active.
- **Engineer B's coupling-operator task builds only the fixed matrix `buildCouplingOperatorEigenbasis` (gauge dispatch → dense Hermitian $O_\text{eig}$)** — already exactly what the adopted method needs. The closure-producing `TimeDependentMatrix`/`makeInteractionHamiltonian` shape a CN-family implementation would have wanted is **not** built; Engineer A consumes B's coupling-operator matrix and field-evaluation function directly as two plain-function handoffs instead of one fused closure.
- No architectural provision is made for a general-banded complex LAPACK solver, a `TimeDependentMatrix`-style closure interface, or a reusable generic Crank-Nicolson primitive as production code — all deferred pending the revisit trigger below, except Option 3 (see below).
- **Option 1b is the well-specified target if this decision is ever revisited** — kept in full above because the equivalence proofs and banded-solver argument are exactly the artifacts a future implementer would need, and shouldn't have to be re-derived. **Revisit trigger:** a future coupling model breaks the "fixed operator × time-varying scalar envelope" structure the adopted method's efficiency depends on — e.g. multi-color/chirped fields with independently time-evolving operator character, or any extension beyond ADR-0003's current 1D/single-particle/single-dipole scope. At that point, prefer 1b within the CN family over 1a for its banded efficiency (per the flop-count table above), and use Proof (ii)'s equivalence to 1a as 1b's primary correctness cross-check.
- **Option 1a is retained only as 1b's would-be internal reference/cross-validation step**, never as a second, independently-maintained implementation — it has no standalone revisit trigger of its own.
- **Option 3 has no revisit trigger**, unlike 1a/1b: the analysis above found no regime in which it outperforms the adopted method on either accuracy or cost, so revisiting it would require a new argument this ADR doesn't anticipate, not just a changed circumstance. Recorded here so the reasoning isn't silently re-litigated later — mirroring how this project already records *why* the archived Chebyshev/van-Dijk design (`docs/TDSE-original-design/`) was dropped rather than just noting that it was.

## Source

`docs/SDD.md` §5.3.3 Figure 8 (propagator left open), §6.1 (`tdse` config block, `hbar` convention); `docs/superpowers/specs/2026-06-28-config-yaml-schema-design.md` (`tdse.chebyshev_order` removal rationale); `docs/planning/bsplines.md` ("Eigenstates as Linear Transformations," the $O_\text{eig}=C^TO_\text{BS}C$ rotation formula and round-trip diagram); `TISE/BSpline.hpp` (`integral()`'s compact-support behavior); `TISE/tise.cpp` (`fillBandedMatrices`, the `dsbgv_` `extern "C"` binding a future `zgbsv_` binding would mirror); `TISE/time_evolution.hpp`/`.cpp` (`timeEvolveState`, reused verbatim by Option 2); `docs/planning/tdse-task-breakdown.md` (the task assignments this decision drives); **ADR-0001** (`docs/adr/0001-defer-fedvr-basis.md` — its Krylov/Chebyshev $O(N)$-efficiency argument is the same family of reasoning this ADR's flop-count table applies concretely to the dense-vs-banded-vs-precomputed-diagonalization comparison); **ADR-0003** (`docs/adr/0003-defer-multiparticle-3d-extension.md` — the scope boundary this ADR's revisit trigger references).
