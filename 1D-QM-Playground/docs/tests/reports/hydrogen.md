# Test report: `tests/hydrogen.yaml`

**Verified against:** the hydrogen atom radial equation for **ℓ = 1** (p states): Bohr energies $E_n = -1/2n^2$, $n \ge 2$; Laguerre radial functions; regular Coulomb functions $F_1(\eta, kr)$ for the continuum.
**Branch / commit:** `TISE-Generalization` @ `9d0f04b`.
**Date:** 2026-09-04.
**Verdict:** **PASS.** 2p–5p energies and radial functions are reproduced to $10^{-6}$ or better (2p, 3p to $10^{-13}$); higher states show textbook confinement by the $r = 100$ wall. Continuum wavefunctions are the exact Coulomb waves to $< 1$ %. The reported phase shifts are self-consistent to $10^{-7}$ but are "flat-asymptote" values at $R = 100$, not Coulomb phase shifts — a documented limitation, quantified here at 0.02–0.07 rad. The `dδ/dE` column is unusable.

Companion files: `known_solutions_report.html`, `verify_known_solutions.py`, `figures/`.

---

## 1. Input

```yaml
run:        output_dir: "./data/hydrogen"
bspline:    n_nodes: 51,  order: 12,  domain: [0.0, 100.0]
potential:  "{'domain': '(0, 100]', 'function': '-1/x + 1/x^2'}"
tise:       n_pts_eigenstate: 301, error_threshold: 1.0e-10
            continuum: enabled: true, E_threshold: 0.0, E_max: 0.5, n_energies: 5, n_pts: 500
visualization: eigenstates: true
```

Node spacing $h = 2$, $Z = 1$, box radius $R = 100$ bohr.

## 2. What the solver actually poses — and which hydrogen problem this is

Separating the three-dimensional Coulomb problem in spherical coordinates with $\psi(\mathbf r) = \dfrac{u(r)}{r}\,Y_{\ell m}(\theta,\varphi)$ leaves a one-dimensional equation for the reduced radial function:

$$-\tfrac12 u'' + \Big(-\frac{1}{r} + \frac{\ell(\ell+1)}{2r^2}\Big)u = E\,u, \qquad u(0) = 0.$$

The $1/r^2$ term is the centrifugal barrier, the price of carrying angular momentum. The YAML's potential $-1/x + 1/x^2$ has centrifugal coefficient $\ell(\ell+1)/2 = 1$, i.e. **ℓ = 1**. This test is therefore the **p-wave** radial problem, not the s-wave one; there is no $-0.5$ (1s) in the spectrum and there should not be. The solver's dropped edge splines supply precisely the $u(0) = 0$ condition the radial equation needs, plus a wall at $r = 100$. `eigenstate_N.png` plots $|u|^2$, the radial probability density.

## 3. The known solution, step by step

### 3.1 Bound states

**Step 1 — asymptotics.** For $E < 0$ and large $r$, $u \sim e^{-\kappa r}$ with $\kappa = \sqrt{-2E}$. Near $r \to 0$ the centrifugal term dominates and the regular solution is $u \sim r^{\ell+1} = r^2$.

**Step 2 — factor and series.** Write $u = r^{\ell+1}e^{-\kappa r}\,L(r)$. The series for $L$ terminates only when $1/\kappa = n$ is an integer with $n \ge \ell + 1$; the polynomial is an associated Laguerre polynomial. Hence

$$E_n = -\frac{1}{2n^2}, \qquad n = \ell+1, \ell+2, \dots = 2, 3, 4, \dots \text{ for } \ell = 1:$$

$$E = -0.125,\ -0.0\overline{5},\ -0.03125,\ -0.02,\ -0.013\overline{8},\ -0.0102041,\ -0.0078125,\ \dots$$

**Step 3 — radial functions.**

$$u_{n1}(r) = r\,R_{n1}(r), \qquad R_{n\ell}(r) = \sqrt{\Big(\frac{2}{n}\Big)^3\frac{(n-\ell-1)!}{2n\,(n+\ell)!}}\; e^{-r/n}\Big(\frac{2r}{n}\Big)^{\ell} L^{2\ell+1}_{n-\ell-1}\!\Big(\frac{2r}{n}\Big),$$

normalised so that $\int_0^\infty u^2\,dr = 1$. $u_{n1}$ has $n - 2$ radial nodes. Landmarks visible in the plots: $|u_{2p}|^2$ peaks at $r = 4$; $|u_{3p}|^2$ has one node and its outer peak near $r = 12$; $\langle r\rangle_{n\ell} = \tfrac12[3n^2 - \ell(\ell+1)]$ gives 5, 12.5, 23, 36.5, 53, 72.5 for 2p–7p.

**Step 4 — box confinement.** A state with $\langle r\rangle \approx 1.5n^2$ and an $e^{-r/n}$ tail needs a box several times $n^2$. With $R = 100$: 2p–5p fit comfortably; 6p ($\langle r\rangle = 53$, tail $e^{-r/6}$ still $\sim 10^{-3}$ of its peak at $r = 100$) is slightly squeezed; 7p, 8p and above are strongly squeezed and their energies rise toward, and then above, zero. Confinement always *raises* the energy (the wall removes the low-kinetic-energy tail).

### 3.2 Continuum

For $E > 0$, $k = \sqrt{2E}$, the regular solution of the Coulomb radial equation is the Coulomb wave function $F_\ell(\eta, kr)$ with Sommerfeld parameter $\eta = -Z/k = -1/k$ (negative for attraction). Its asymptotic form is

$$F_\ell(\eta,\rho) \to \sin\!\big(\rho - \ell\tfrac{\pi}{2} - \eta\ln 2\rho + \sigma_\ell\big), \qquad \sigma_\ell = \arg\Gamma(\ell + 1 + i\eta),$$

with $\rho = kr$. The logarithmic term means a Coulomb wave never settles to a fixed phase relative to $\sin kr$: the phase keeps drifting, slowly, forever. "The phase shift" of a pure Coulomb potential is therefore only defined *relative to Coulomb functions* (where it is zero), not relative to free waves. Energy normalisation multiplies $F_\ell$ by $\sqrt{2/\pi k}$.

## 4. What the solver produced

### 4.1 Bound energies

| file index | orbital | $E$ (solver) | $-1/2n^2$ | rel. error |
|---|---|---|---|---|
| 1 | 2p | −0.1249999999999988 | −0.125 | 9.3e-15 |
| 2 | 3p | −0.0555555555555592 | −0.0555555555555556 | 6.5e-14 |
| 3 | 4p | −0.0312499999997379 | −0.03125 | 8.4e-12 |
| 4 | 5p | −0.0199999786301715 | −0.02 | 1.1e-06 |
| 5 | 6p | −0.0138718175739423 | −0.0138888888888889 | 1.2e-03 |
| 6 | 7p | −0.0096531596064853 | −0.0102040816326531 | 5.4e-02 |
| 7 | 8p | −0.0048446185441568 | −0.0078125 | 3.8e-01 |
| 8 | (box) | +0.0013231169449417 | — | — |

Seven negative eigenvalues, then box states. The error jumps from $10^{-12}$ (4p) to $10^{-6}$ (5p) to $10^{-3}$ (6p): this is the onset of box confinement from §3.1 step 4, not basis error — the Coulomb functions are smooth and the $h = 2$ grid resolves them easily (2p's wavelength scale is ~10 bohr). Every confined energy lies above the exact one, as it must.

### 4.2 Radial functions

Compared to $u_{n1}(r)$ on the 301-point output grid (sign fixed):

| orbital | $\int u^2\,dr$ | rel. $L^2$ difference |
|---|---|---|
| 2p | 1.0000002 | 1.1e-08 |
| 3p | 1.0000001 | 1.2e-09 |
| 4p | 1.0000000 | 1.8e-06 |
| 5p | 1.0000000 | 7.4e-04 |
| 6p | 1.0000000 | 3.5e-02 |
| 7p | 1.0000000 | 0.30 |

The output images (`figures/output_hydrogen.png`) show $|u_{2p}|^2$ peaking at $r \approx 4$ with the correct $r^4$ rise from the origin, one more node per state, outer peaks at ~12, 24, 39, 58, 77 bohr, and — for 6p and 7p — the density forced to zero at $r = 100$ before its natural tail has decayed. The overlays (`figures/hydrogen_states.png`) are indistinguishable for 2p–4p and show the wall bending the 6p tail down early. (The $2\times10^{-7}$ deviation of the 2p norm from 1 is trapezoid-rule error on the steep $r^2$ rise near the origin, sampled at 0.33-bohr spacing; the $S$-normalisation itself is exact.)

### 4.3 Continuum and phase shifts

Two independent comparisons were made, so that "is the wavefunction right?" is separated from "is the flat-asymptote phase-shift formula appropriate?":

1. The solver's $\psi_E(r)$ was laid over $\sqrt{2/\pi k}\,F_1(-1/k, kr)$ computed with `mpmath.coulombf` at 30 digits.
2. The solver's own matching formula, $\delta = \arctan(k\psi/\psi')|_{R} - kR$, was applied to the *exact* $F_1$ and its derivative at $R = 100$.

| E | $\delta \bmod \pi$ (solver) | same formula on exact $F_1$ | difference | true Coulomb phase $-\tfrac{\pi}{2} - \eta\ln 2kR + \sigma_1$ (mod π) | $d\delta/dE$ (solver) | rel. $L^2$ err vs $F_1$ |
|---|---|---|---|---|---|---|
| 0.1 | 0.807586 | 0.807586 | 4.8e-09 | 0.7325 | +211.4 | 2.4 % |
| 0.2 | −1.036321 | −1.036322 | 2.2e-07 | −1.0830 | +1.32 | 0.36 % |
| 0.3 | 1.150641 | 1.150642 | −4.2e-07 | 1.1283 | −6.64 | 0.80 % |
| 0.4 | 0.553418 | 0.553425 | −7.3e-06 | 0.5305 | −2.83 | 0.38 % |
| 0.5 | 0.120704 | 0.120710 | −5.9e-06 | 0.1022 | −3.37 | 0.33 % |

Interpretation:

* **The numerical continuum wave is the Coulomb function.** The solver's phase shift agrees with the exact Coulomb function pushed through the same formula to $10^{-8}$–$10^{-6}$, and the wavefunctions agree to $< 1$ % for $E \ge 0.2$ (`figures/hydrogen_continuum.png` shows the characteristic shortening of the wavelength toward the nucleus reproduced exactly). The 2.4 % at $E = 0.1$ is mostly the slow amplitude drift of a Coulomb wave, which the flat normalisation $\sqrt{2/\pi k}$ at $R = 100$ does not track.
* **The reported $\delta$ is a box-dependent quantity.** It differs from the true Coulomb asymptotic phase at $R = 100$ by 0.02–0.07 rad, because at $r = 100$ the potential is still $-0.01$ — a tenth of $E = 0.1$ — so the wave is not yet a free sinusoid there. This is the "flat-asymptote matching assumes a regular right boundary" caveat in `TISE/README.md`, now made quantitative. Bachau et al. (2001) match to Coulomb functions $F_\ell, G_\ell$ for exactly this reason.
* **$d\delta/dE$ is meaningless here.** The value +211 at $E = 0.1$ arises from a central difference of $\sin 2\delta$ over $\Delta E = 0.1$ divided by $2\cos 2\delta$ with $2\delta \approx 1.615$, $\cos 2\delta \approx -0.044$. Both the coarse difference and the division are unsound.

The solver also warned that $E = 0.1, 0.2, 0.3, 0.5$ each lie within $< 0.008$ of a box eigenvalue. As in the free-particle test this is unavoidable with a $\Delta E = 0.1$ grid in a 100-bohr box, and it did not harm the continuum states.

## 5. Correctness assessment

| Quantity | Status |
|---|---|
| $E_{2p}, E_{3p}, E_{4p}$ | correct to $10^{-14}$–$10^{-11}$ |
| $E_{5p}$ | correct to $10^{-6}$ (onset of box effect) |
| $E_{6p}$ and above | box-confined; qualitatively correct (raised), not comparable to free-atom values |
| $u_{2p}$–$u_{4p}$ | correct to $< 2\times10^{-6}$ |
| identification of the problem | ℓ = 1, not ℓ = 0 — the YAML is correct for p states; a test named "hydrogen" that omits 1s may surprise a reader |
| continuum $\psi_E$ | correct Coulomb waves to $< 1$ % ($E \ge 0.2$) |
| $\delta \bmod \pi$ | self-consistent to $10^{-7}$; differs from Coulomb phase by 0.02–0.07 rad by construction |
| $d\delta/dE$ | **unusable** |

**Overall: PASS** for everything the solver claims to compute; the phase shifts must be read with the flat-asymptote caveat.

## 6. Recommendations

1. **Say what the test is.** Rename or document the YAML as `hydrogen_l1` / "hydrogen p-wave", and add an ℓ = 0 companion (`function: '-1/x'`, expected $-0.5, -0.125, -0.0556, \dots$) so the ground state of hydrogen is actually certified. The two together also exercise the $r^{\ell+1}$ behaviour at the origin for two different ℓ.
2. **Coulomb matching for Coulomb tails.** For potentials with a $-Z/r$ tail, `matchAsymptotic` should match to $F_\ell(\eta,kr)$ and $G_\ell(\eta,kr)$ (and their derivatives) at $R$ rather than to $\sin, \cos$. Then $\delta$ for pure hydrogen is zero by construction and, for hydrogen-like model potentials, becomes the physically meaningful non-Coulomb phase shift. `classifyAsymptote` already detects the tail type; this is the natural place to branch. The `PHY5606_F25_ContinuumEigenstates` notes and Bachau et al. §4 give the formulas.
3. **Enlarge the box or use a non-uniform grid** for bound-state tests above 5p: $R = 400$ with the same $h$ certifies through ~10p; the strategic-node-placement path in `tise.cpp` (dense near the origin, sparse outside) would do so without quadrupling the basis.
4. **Fix $d\delta/dE$** as in the finite-well report: fine local energy grid with branch unwrapping, or analytic differentiation of the matching condition; remove the $\sin 2\delta/\cos 2\delta$ construction.
5. **Plot/ship only meaningful states.** Of the 59 eigenstate PNGs, seven are bound and five of those are free-atom-accurate. Gate on $E < 0$ (or on a confinement diagnostic: $|u(R - \text{few bohr})|$ above a threshold flags a squeezed state).
6. **Automated regression.** Assert $E_{2p,3p,4p}$ within $10^{-10}$ of $-1/2n^2$ and $E_{5p}$ within $10^{-5}$; assert $\int u_{2p}^2 = 1$ within $10^{-6}$ and the position of the 2p density peak at $r = 4 \pm h$; assert $\delta \bmod \pi$ within $10^{-5}$ of the flat-formula-on-$F_1$ value. The energy, norm and phase-shift checks are implemented in `verify_known_solutions.py` (which needs `mpmath` for the Coulomb functions).
