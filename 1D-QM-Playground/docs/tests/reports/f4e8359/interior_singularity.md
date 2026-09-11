# Test report: `tests/interior_singularity.yaml` — iteration 3

**Verified against:** the exact split-domain spectrum — particle in a box on $[0,20]$ ($E = n^2\pi^2/800$) ∪ repulsive-Coulomb box on $[20,40]$ (zeros of $F_0(1/k, 20k)$).
**Branch / commit:** `TISE-Generalization` @ `f4e8359`. **Previous:** `8236239` (FAIL, 41 %).
**Date:** 2026-09-11.
**Verdict:** **PASS — the iteration-2 defect is fixed exactly as recommended.** All twelve lowest eigenvalues agree with the exact spectrum to $10^{-13}$–$10^{-11}$ (previously 41 % high on the field-free side). The basis has 50 states, the number predicted for "multiplicity $k-1$ at the singular point plus one dropped B-spline". The dead zone is gone: eigenvectors are non-zero right up to $x=20$ on their own side. One item remains open: continuum states and phase shifts are still written for this split domain, where they have no meaning.

Companion files: `verify_new_tests.py`, `figures/interior_singularity.png`.

---

## 1. Input

Unchanged: 41 nodes, order 8, $[0,40]$, $V=0$ on $[0,20)$ and $1/(x-20)$ on $(20,40]$, continuum enabled at $E=0.1\ldots0.5$, `run_analysis: false`.

## 2. What changed in the solver (commit `0bf17be`)

`buildStrategicGridAndDropSet` now treats a Singular interior join like a domain edge: the knot at $x_s=20$ is given multiplicity $\text{order}-1 = 7$, which makes the basis merely $C^0$ there and leaves exactly one B-spline non-zero at $x_s$; that one is dropped. Basis: $47 + 6$ extra knots $- 1 - 1 - 1 = 50$ (was $47 - 9 - 2 = 36$). The e2e test now asserts $E_1 = \pi^2/800$ and $E_3 = 0.1005405216$ (first Coulomb-side root).

## 3. The known solution

`../8236239/interior_singularity.md` §3: the one-sided $1/s$ singularity is non-integrable, so $\phi(20)=0$ is forced and the two halves decouple. Left: $E_n = 0.0123370\,n^2$. Right: $F_0(1/k, 20k)=0 \Rightarrow 0.100541, 0.168654, 0.254096, 0.359297, 0.485622, \dots$

## 4. What the solver produced

| # | $E$ (f4e8359) | exact | side | rel. error | at 8236239 |
|---|---|---|---|---|---|
| 1 | 0.0123370055013608 | 0.0123370055013617 | box | 7e-14 | +41 % |
| 2 | 0.0493480220054468 | 0.0493480220054468 | box | 1e-15 | +41 % |
| 3 | 0.1005405215766150 | 0.1005405215766166 | Coulomb | 1.6e-14 | +0.09 % |
| 4 | 0.1110330495122504 | 0.1110330495122553 | box | 4e-14 | +41 % |
| 5 | 0.1686540851059307 | 0.1686540851059007 | Coulomb | 1.8e-13 | +1.0 % |
| 6 | 0.1973920880217923 | 0.1973920880217872 | box | 2.6e-14 | — |
| 7 | 0.2540957784525549 | 0.2540957784527245 | Coulomb | 6.7e-13 | — |
| 8 | 0.3084251375344702 | 0.3084251375340424 | box | 1.4e-12 | — |

The values coincide with the iteration-2 scipy demonstration (`proper_bspline` column there) to all printed digits, as they must — it is the same discretisation. Each eigenvector lives entirely on one side and now extends to $x=20.0$ (the `eigenstate_NNN.dat` grid point at 20 is the first zero), with no band of forced zeros. `figures/interior_singularity.png`: solver state 1 lies on top of $\sqrt{2/20}\sin(\pi x/20)$; state 3 on top of the normalised $F_0$; the spectrum panel shows solver, exact and the iteration-2 scipy check coinciding.

### 4.1 Continuum output

`phase_shifts.dat` is still written: $\delta = 0.92, -2.02, -4.39, -6.40, -8.24$ with $d\delta/dE \approx -17\ldots-30$. These numbers are the flat-asymptote matching at $x=40$ applied to a state built from two decoupled boxes, with $V(40) = 0.05$ not flat. They describe nothing physical; iteration-2 recommendation 4 (refuse or warn for split domains) has not been implemented.

## 5. Correctness assessment

| Quantity | 8236239 | f4e8359 |
|---|---|---|
| left-side (box) eigenvalues | +41 % | **≤ 1.4e-12** |
| right-side (Coulomb) eigenvalues | +0.1–1 % | **≤ 7e-13** |
| eigenfunctions near $x=20$ | zero on [19,21] | correct, vanish linearly at 20 |
| basis size | 36 | 50 (as designed) |
| continuum / phase shifts | meaningless, written | meaningless, still written |

**Overall: PASS on the spectrum and eigenfunctions.**

## 6. Recommendations

1. **Refuse or clearly flag continuum construction when an interior Singular join splits the domain.** The detection already exists (`detectPotentialStructure` returns the join); when any Singular join lies strictly inside $(x_{\min}, x_{\max})$ and the continuum is enabled, either skip it with the same message used for right-edge singularities, or build per-subdomain continuum states with the sub-box that touches $x_{\max}$ as the only one matched.
2. **Consider multiplicity $k$ (full)** rather than $k-1$ for one-sided singularities, which decouples the two halves in the basis itself and permits different node densities per side. Not needed for correctness — the present result is at round-off — but it would let the Coulomb side, whose wavefunctions are stiffer near $s=0$, be refined independently.
3. Keep the two analytic assertions in the e2e test; add $E_5 = 0.1686540851$ (second Coulomb-side root) so that both sides are guarded by a state above the first.
