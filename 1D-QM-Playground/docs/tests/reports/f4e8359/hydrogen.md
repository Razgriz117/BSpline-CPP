# Test report: `tests/hydrogen.yaml` — iteration 3

**Verified against:** the hydrogen radial equation for $\ell = 1$: $E_n = -1/2n^2$ ($n\ge2$), Laguerre radial functions, and — now that the solver matches to Coulomb functions — the exact result that a pure Coulomb potential has **zero** phase shift relative to $F_\ell, G_\ell$. Reference Coulomb functions from `mpmath` at 25–30 digits.
**Branch / commit:** `TISE-Generalization` @ `f4e8359`. **Previous:** `8236239`, `9d0f04b`.
**Date:** 2026-09-11.
**Verdict:** **PASS, with a quantified accuracy floor.** Bound states unchanged (2p–5p to $\le10^{-6}$, 6p+ box-confined and correctly flagged). The new Coulomb-tail matching (commit `f4e8359`) produces $\delta = 0.0012$–$0.0023$ rad against an exact value of 0 — a 10–50× improvement on the flat-asymptote values (0.02–0.07 rad off) and the continuum wavefunctions now match $\sqrt{2/\pi k}\,F_1(\eta,kr)$ to $2\times10^{-5}$–$2\times10^{-3}$. The residual $2\times10^{-3}$ rad is **entirely** the hand-rolled Coulomb-function evaluator, not the B-spline solution: feeding the exact $F_1$ through the solver's own $(F,G)$ pair reproduces the reported $\delta$ to $10^{-5}$. Its Wronskian is $-0.9982$ instead of $-1$, and it converges as $O(h^2)$, not the $O(h^4)$ Numerov should give.

Companion files: `verify_known_solutions.py`, `coulomb_wave_check.cpp`, `coulomb_wave_stepsize.cpp`, `figures/`.

---

## 1. Input

```yaml
potential:  "{'domain': '(0, inf)', 'function': '-1/x + 1/x^2'}"     # was (0, 100]
tise.continuum:  ... l: 1                                             # new field
```

Two changes since 8236239. Declaring the potential on $(0,\infty)$ is what makes `classifyAsymptote` run on the right edge (it fits $V\sim -0.999996/x$, classifies **Coulomb**), and `tise.continuum.l: 1` tells `matchAsymptotic` which $F_\ell, G_\ell$ to use. Everything else (51 nodes, order 12, $R=100$, five energies 0.1–0.5) is unchanged.

## 2. What the solver actually poses

Bound problem: unchanged — $u(0)=u(100)=0$, $\ell=1$ centrifugal term baked into $V$. Continuum: $\psi_E$ is built as before, but at $R=100$ it is now matched to

$$\psi_E(R) = A_E\big[\cos\delta\,F_1(\eta,kR) + \sin\delta\,G_1(\eta,kR)\big], \qquad \eta = C/k = -1/k,$$

with $F_1, G_1$ from `evaluateCoulombFunctions` (WKB-corrected asymptotic start at $50\times$ the requested $\rho$, Numerov integration inward with step 0.1), and $A_E$ fixed by the Wronskian so that $\psi_E$ is energy-normalised. For a *pure* Coulomb potential the exact $\psi_E$ **is** $F_1$, so the exact answer for this test is $\delta(E) \equiv 0$ and $d\delta/dE \equiv 0$ — a much sharper target than iteration 1's "self-consistent to $10^{-7}$ but $0.02$–$0.07$ rad from the Coulomb phase".

## 3. The known solution

Bound states and radial functions: `../9d0f04b/hydrogen.md` §3.1. Continuum: the regular solution is $F_1(\eta, kr)$ with $\eta=-1/k$; the Coulomb phase $\sigma_1 = \arg\Gamma(2+i\eta)$ is absorbed into the definition of $F, G$, so the short-range phase shift of hydrogen relative to Coulomb functions is exactly zero at every energy.

## 4. What the solver produced

### 4.1 Bound states — unchanged

2p, 3p, 4p to $10^{-14}$–$10^{-11}$; 5p to $1.1\times10^{-6}$; 6p–8p box-confined and flagged by `checkWellContainment` ($\psi'(100) = 0.0025, 0.0104, 0.0188$). Radial functions to $\le 2\times10^{-6}$ for 2p–4p.

### 4.2 Continuum — the new result

| E | $\delta$ (solver) | exact | $d\delta/dE$ (solver) | exact | rel. $L^2$ err vs $\sqrt{2/\pi k}F_1$ | at 8236239 (flat) |
|---|---|---|---|---|---|---|
| 0.1 | 0.002088 | 0 | +0.336 | 0 | 2.1e-05 | δ off by 0.075; L2 2.4 % |
| 0.2 | 0.002255 | 0 | −0.085 | 0 | 1.1e-03 | 0.047; 0.36 % |
| 0.3 | 0.001186 | 0 | −0.212 | 0 | 1.7e-03 | 0.022; 0.80 % |
| 0.4 | 0.002007 | 0 | +0.102 | 0 | 3.4e-04 | 0.023; 0.38 % |
| 0.5 | 0.001982 | 0 | −0.093 | 0 | 1.5e-03 | 0.019; 0.33 % |

The wavefunction agreement at $E=0.1$ improved from 2.4 % to $2\times10^{-5}$ because the normalisation $A_E$ is now taken from the Coulomb Wronskian rather than from $\sqrt{2/\pi k}$ against a flat wave — the round-1 report attributed that 2.4 % to "slow amplitude drift of a Coulomb wave", which is exactly what the new matching accounts for.

### 4.3 Where the 0.002 rad comes from

`coulomb_wave_check.cpp` calls the solver's `evaluateCoulombFunctions(l=1, η, kR)` at the five energies and compares with `mpmath.coulombf/coulombg`:

| E | $F$ error | $G$ error | $F'$ error | $G'$ error | $W = FG'-F'G$ | $\delta$ obtained by pushing the **exact** $F_1$ through the solver's $(F,G)$ |
|---|---|---|---|---|---|---|
| 0.1 | −6.8e-05 | 2.0e-03 | 2.1e-03 | 1.9e-03 | −0.99817 | 0.00206 |
| 0.2 | −1.1e-03 | −8.3e-04 | −2.3e-03 | 1.1e-04 | −0.99825 | 0.00223 |
| 0.3 | −1.3e-03 | 1.1e-04 | −1.7e-03 | 1.5e-03 | −0.99828 | 0.00116 |
| 0.4 | −5.7e-04 | −1.2e-03 | −1.9e-03 | −9.5e-04 | −0.99829 | 0.00199 |
| 0.5 | 1.2e-03 | 5.3e-04 | 2.1e-03 | −6.1e-04 | −0.99830 | 0.00197 |

The last column reproduces the solver's reported $\delta$ (0.00209, 0.00225, 0.00119, 0.00201, 0.00198) to $10^{-5}$: **the B-spline continuum state is correct to $\sim10^{-5}$ rad; all of the residual is the Coulomb-function evaluator.** ($\sigma_1$ from `coulombPhaseShift` is exact to $10^{-13}$.) The Wronskian, which must be exactly $\pm1$, is off by $1.8\times10^{-3}$ — the same size as the phase error, and a diagnostic the code could compute for free.

`coulomb_wave_stepsize.cpp` varies the evaluator's `stepSize` and `farMultiplier` at $E=0.1$:

| far | h | $F$ | $G$ | $W$ | time |
|---|---|---|---|---|---|
| 50 | 0.1 | 0.9763645 | 0.0240391 | −0.998169 | 0.5 ms |
| 50 | 0.05 | 0.9763933 | 0.0235914 | −0.999542 | 0.9 ms |
| 50 | 0.025 | 0.9763962 | 0.0235621 | −0.999886 | 1.9 ms |
| 50 | 0.0125 | 0.9763965 | 0.0235601 | −0.999971 | 3.7 ms |
| 100 | 0.0125 | 0.9764147 | 0.0227966 | −0.999971 | 7.8 ms |
| **exact** | | **0.9764322** | **0.0220322** | **−1** | |

Two independent error sources are visible. The Wronskian error falls by $4\times$ per halving of $h$ — $O(h^2)$, so the derivative is being formed with a second-order formula and the fourth-order accuracy of Numerov is being thrown away at the end. And at fixed $h\to0$ the values still depend on `farMultiplier` (0.976397 at 50, 0.976415 at 100, 0.976432 exact): the WKB-corrected asymptotic *starting values* carry an $O(1/\rho_{\rm far})$ phase error that the inward integration faithfully preserves. $G$, being small here, is the more sensitive of the two (3 % off at the shipped settings).

## 5. Correctness assessment

| Quantity | Status |
|---|---|
| bound energies / radial functions | correct (unchanged) |
| box-confinement diagnostic | correct (6p, 7p, 8p flagged) |
| Coulomb-tail classification, $C=-0.999996$ | correct |
| continuum $\psi_E$ and normalisation | correct to $2\times10^{-5}$–$2\times10^{-3}$ |
| $\delta$ (should be 0) | $1$–$2\times10^{-3}$ rad — limited by the $F/G$ evaluator, not the solver |
| $d\delta/dE$ (should be 0) | $\pm0.1$–$0.3$: the $2\times10^{-3}$ $\delta$ error divided by $2\,\text{fineDE}=2\times10^{-3}$ |
| `phase_shifts.png` | not produced (`visualization.phase_shifts` not set in this YAML) |

**Overall: PASS.** The physics is right and the accuracy floor is understood and localised.

## 6. Recommendations

1. **Compute and report the Wronskian** $FG' - F'G$ in `evaluateCoulombFunctions` and warn if $|W \mp 1| > 10^{-6}$. It costs nothing and would have exposed both issues below.
2. **Form the derivative to Numerov order.** With $u_{n\pm1}$ and the ODE's $f(\rho)=1-2\eta/\rho-\ell(\ell+1)/\rho^2$, the central difference $D = (u_{n+1}-u_{n-1})/2h = u' + \tfrac{h^2}{6}u''' + O(h^4)$ can be corrected using $u''' = -(fu)' = -f'u - fu'$, giving
$$u'_n = \frac{D + \tfrac{h^2}{6}\,f'_n u_n}{1 - \tfrac{h^2}{6}\,f_n} + O(h^4).$$ This turns the $O(h^2)$ Wronskian error into $O(h^4)$ at no extra integration cost.
3. **Replace the WKB start by the asymptotic series** (Abramowitz & Stegun 14.5.1–14.5.10, or NIST DLMF 33.11) evaluated to convergence at $\rho_{\rm far}$; its error is controlled and exponentially small for $\rho_{\rm far}\gg\eta^2$, eliminating the `farMultiplier` dependence. Alternatively, evaluate $F,G$ directly by Steed's continued-fraction algorithm (Barnett's COULFG), which is ~100 lines and gives $10^{-12}$ everywhere; or link GSL's `gsl_sf_coulomb_wave_FG_e`. Any of these would bring $\delta$ for hydrogen to $\lesssim10^{-6}$ and the derivative to $\lesssim10^{-3}$.
4. **Enable `visualization.phase_shifts: true`** in this YAML now that the plot is meaningful, and assert `abs(delta) < 5e-3` (tightening to `1e-5` once item 3 lands) in the e2e test.
5. **Name the test `hydrogen_l1`** and add an $\ell=0$ companion — carried from iteration 1; the new `tise.continuum.l` field makes the mismatch risk concrete (an `l` inconsistent with the centrifugal term in `potential` would silently give a nonzero "phase shift").
6. **Sign convention for plotted states.** `eigenstate_1.png` now plots the raw 2p function, which comes out negative (eigenvector sign is arbitrary). Normalising each state so that its first antinode is positive would make plots comparable across runs.
