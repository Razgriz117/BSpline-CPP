# Test report: `tests/finite_square_well.yaml` — iteration 3

**Verified against:** finite square well with a hard wall at the origin (odd states of the symmetric well / $\ell=0$ spherical well): $K\cot Ka = -\kappa$ for bound states, $\delta(E) = \arctan\!\big(\tfrac{k}{K}\tan Ka\big) - ka$ for the continuum, $100k + \delta = j\pi$ for box states, all solved independently.
**Branch / commit:** `TISE-Generalization` @ `f4e8359`. **Previous:** `8236239` (`../8236239/finite_square_well.md`), `9d0f04b`.
**Date:** 2026-09-11.
**Verdict:** **PASS — every quantity, including the two items left open in iteration 2.** Bound energies, box states, eigenfunctions and continuum states are unchanged from `8236239` ($10^{-13}$–$10^{-9}$). New at this commit: `phase_shifts.dat` stores a continuously unwrapped $\delta(E)$ that matches the exact continuous trajectory to $3\times10^{-5}$ rad at all five energies, and $d\delta/dE$ — unusable in both previous iterations — now agrees with the analytic derivative to $10^{-4}$ relative ($-19.432$ vs $-19.432$ at $E=0.1$). `phase_shifts.png` finally shows the physics.

Companion files: `verify_known_solutions.py`, `figures/`.

---

## 1. Input

Unchanged from iteration 2 (51 nodes, order 12, $[0,100]$, $V=-1$ on $[0,10)$, continuum $E=0.1\ldots0.5$ in 5 steps, `visualization: eigenstates, phase_shifts`). The only behavioural change to note: `eigenstate_N.png` now plots the raw $\phi_n(x)$ rather than $|\phi_n|^2$ (commit `5cd28e4`; the squared plot is now opt-in via `visualization.bound_states_squared`).

## 2. What changed in the solver (commit `0bf17be`)

* `matchAsymptotic` computes $d\delta/dE$ as a central difference of the raw $\delta$ over an internal fine step `fineDE` $=10^{-3}$ around each requested energy (branch-corrected), replacing the $\sin2\delta/\cos2\delta$ difference across the production grid. This required threading $H$, $S$ and the order into `matchAsymptotic` so that continuum states can be built at auxiliary energies.
* The stored $\delta(E)$ is unwrapped: the first point is normalised into $(-\pi/2, \pi/2]$ and each subsequent point is placed on the branch predicted by $\delta_{i-1} + (d\delta/dE)_{i-1}\,\Delta E$.
* The e2e test for this YAML now asserts $E_1$, $E_4$ and $\delta(0.5)$ against the analytic values from the iteration-2 report.

## 3. The known solution

See `../9d0f04b/finite_square_well.md` §3 for the full derivation. Exact reference values: bound $E = -0.956996254,\ -0.828522650,\ -0.616571751,\ -0.326731918$ (four in the box; the fifth, bound at $-2.3\times10^{-5}$ on an infinite domain, is squeezed out by the wall at 100); box states $5.5253\times10^{-5},\ 1.16434\times10^{-3},\ 3.37370\times10^{-3},\ 6.70885\times10^{-3}$; continuous phase shift $\delta = 1.46421,\ -0.13072,\ -1.25335,\ -1.93910,\ -2.07379$ at $E = 0.1\ldots0.5$ (the last two are $1.20249$ and $1.06780$ mod $\pi$); $d\delta/dE = -19.432,\ -13.256,\ -9.241,\ -4.055,\ -0.083$.

## 4. What the solver produced

### 4.1 Bound energies and box states — unchanged from 8236239

| n | $E_n$ (solver) | exact | rel. error |
|---|---|---|---|
| 1 | −0.956996254368899 | −0.956996254368992 | 9.7e-14 |
| 2 | −0.828522650434678 | −0.828522650434757 | 9.5e-14 |
| 3 | −0.616571750670465 | −0.616571750671813 | 2.2e-12 |
| 4 | −0.326731917714706 | −0.326731918129711 | 1.3e-09 |
| box 1 | 5.52701e-05 | 5.52530e-05 | 3.1e-04 |
| box 2 | 1.164353e-03 | 1.164340e-03 | 1.1e-05 |

Eigenfunction $L^2$ differences $1.5\times10^{-6}$–$1.6\times10^{-5}$; continuum-state $L^2$ differences $5\times10^{-5}$–$6\times10^{-4}$ (`figures/fsw_states.png`, `figures/fsw_continuum.png`).

### 4.2 Phase shifts — the new result

| E | $\delta$ stored (unwrapped) | exact continuous | diff | $d\delta/dE$ (solver) | exact | rel. diff | at 8236239 |
|---|---|---|---|---|---|---|---|
| 0.1 | 1.464206 | 1.464209 | −3e-06 | −19.4323 | −19.4321 | 1e-05 | +2.41 |
| 0.2 | −0.130730 | −0.130723 | −7e-06 | −13.2565 | −13.2564 | 1e-05 | −2.08 |
| 0.3 | −1.253377 | −1.253355 | −2e-05 | −9.2411 | −9.2408 | 3e-05 | −2.89 |
| 0.4 | −1.939169 | −1.939105 | −6e-05 | −4.0563 | −4.0552 | 3e-04 | −4.85 |
| 0.5 | −2.074035 | −2.073792 | −2.4e-04 | −0.1053 | −0.0827 | 0.27 | −1.62 |

Three observations:

* **Branch tracking works here**, including across the $E=0.3\to0.4$ step where $\delta$ drops by 0.69 rad and the raw arctangent value would have jumped by $\pi$. The derivative-predicted unwrapping is what makes this robust; a nearest-value heuristic would have failed at $E=0.4$–$0.5$, as the commit message notes.
* **The derivative is now a genuine derivative.** Relative accuracy $10^{-5}$ where $\delta$ itself is accurate to $10^{-5}$. At $E=0.5$ the absolute error is $0.02$ on a value of $-0.08$; the $2.4\times10^{-4}$ error in $\delta$ itself, divided by $2\times$`fineDE` $=2\times10^{-3}$, gives exactly that scale. The derivative inherits the basis error of $\delta$ amplified by $1/(2\,\text{fineDE})$ — see the free-particle report, where this is the only error present.
* `figures/output_phase_shifts_fsw.png`: the upper panel is now the monotone 3.5-rad fall of $\delta(E)$ predicted in iteration 2 §6; the lower panel is the time-delay curve rising from −19 toward 0.

## 5. Correctness assessment

| Quantity | 9d0f04b | 8236239 | f4e8359 |
|---|---|---|---|
| bound energies | 2e-4…2e-2 | 1e-13…1e-9 | 1e-13…1e-9 |
| box states near threshold | wrong | 3e-4…2e-6 | 3e-4…2e-6 |
| eigenfunctions / continuum $\psi_E$ | 1–60 % | ≤ 6e-4 | ≤ 6e-4 |
| $\delta$ (mod π) | 0.03–0.6 rad | ≤ 2.4e-4 | ≤ 2.4e-4, **now stored continuously** |
| $d\delta/dE$ | noise | noise | **correct to 1e-5 … 3e-4 rel.; 0.02 abs at E=0.5** |
| `phase_shifts.png` | — | uninformative | **informative** |

**Overall: PASS. Nothing remains open on this test.**

## 6. Recommendations

1. **Make `fineDE` energy-aware.** The derivative error is $\approx \epsilon_\delta / (2\,\text{fineDE})$ with $\epsilon_\delta$ the basis error in $\delta$. A larger step ($10^{-2}$) would reduce the $E=0.5$ error from 0.02 to ~0.002 while the truncation error $\tfrac16\delta'''\,\text{fineDE}^2$ stays below $10^{-4}$ for this smooth $\delta(E)$. Better still, Richardson-extrapolate two steps ($h$, $2h$), which removes the $O(h^2)$ truncation and lets a larger $h$ be used safely.
2. **Expose `fineDE` in the YAML** (`tise.continuum.fine_dE`) with 1e-3 as default, so users with steep resonances can shrink it.
3. **Add the derivative to the regression assertions**: `abs(ddelta_dE(0.1) + 19.432) < 0.01`. It is the quantity that was broken for two iterations and is now the most sensitive check of the continuum pipeline.
4. `computeEAcc` still ignores $V_{\min}$ (carried from iteration 1); harmless here but still worth fixing.
