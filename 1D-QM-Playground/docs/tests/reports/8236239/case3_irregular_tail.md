# Test report: `tests/case3_irregular_tail.yaml` — iteration 2 (new test)

**Verified against:** high-resolution finite-difference eigenvalues (Richardson-extrapolated, $N = 40\,000 / 80\,000$ points) of $-\tfrac12\phi'' + V\phi = E\phi$ on $[0.1, 50]$ with hard walls, for both the raw potential $V = x^{-3/2}$ and the tapered potential the solver actually uses.
**Branch / commit:** `TISE-Generalization` @ `8236239`. First appearance of this test.
**Date:** 2026-09-04.
**Verdict:** **PASS as a numerical test; the modelling choice it exercises deserves a second look.** The solver's twelve lowest eigenvalues agree with the reference for the *tapered* potential to $10^{-12}$–$4\times10^{-8}$, so the Case-3 window is applied exactly as documented and the eigensolve is correct. Those same eigenvalues differ from the untapered problem by $1.5\times10^{-5}$–$1.5\times10^{-4}$ Ha — a systematic shift the taper introduces on purpose. For a run with `continuum: enabled: false`, as here, the taper buys nothing and costs four digits.

Companion files: `verify_new_tests.py`, `figures/case3_irregular_tail.png`.

---

## 1. Input

```yaml
run:        run_analysis: false,  output_dir: "./data/case3_irregular_tail"
bspline:    n_nodes: 41,  order: 6,  domain: [0.1, 50.0]
potential:  "{'domain': '(0.1, inf)', 'function': '1/x^1.5'}"
tise:       n_pts_eigenstate: 301
            continuum: enabled: false
```

Node spacing $h = 1.2475$; order 6 (lower than the other tests). The potential is declared on $(0.1,\infty)$, i.e. beyond the box, which is what triggers asymptote classification.

## 2. What the solver actually poses

Because the potential's domain extends past $x_{\max}=50$, `tise_solver` calls `classifyAsymptote` on the right side. It fits a power law to the tail, finds $p = 1.5$, and — since $p$ is neither "flat" ($V\to$ const fast enough) nor Coulomb ($p=1$) — classifies the tail as **Case 3, Irregular**, with `recommendedTransitionWidth` $= 0.1\,(x_{\max}-x_{\min}) = 4.99$. New at 8236239, this classification is now *acted on*: `fillBandedMatrices` evaluates

$$\tilde V(x) = W(x)\,V(x), \qquad W(x) = \begin{cases} 1 & x \le R-\Delta \\ \sin^2\!\big(\tfrac{\pi}{2}\,\tfrac{x-R}{\Delta}\big) & R-\Delta < x < R \\ 0 & x \ge R \end{cases}$$

with $R=50$, $\Delta=4.99$: the potential is smoothly switched off over the last 4.99 bohr before the wall ($V(45.01)=0.0033 \to \tilde V(50) = 0$). Both `warnings.json` entries describe this correctly. Hard walls at 0.1 and 50 remain.

## 3. The known solution

There is no closed form for $x^{-3/2}$, so the reference is numerical and independent of B-splines: a second-order finite-difference Hamiltonian on a uniform grid, lowest 12 eigenvalues by shift-invert Lanczos at $N=40\,000$ and $80\,000$ interior points, Richardson-extrapolated in $h^2$. Residual error is estimated at $<10^{-11}$ from the size of the extrapolation correction. Two references were computed: $E^{\rm raw}$ for $V=x^{-3/2}$ and $E^{\rm tap}$ for $\tilde V$.

Qualitatively: $V\ge0$ everywhere, so all states are box states, and $V$ is large only near $x=0.1$ ($V(0.1)=31.6$), which pushes the wavefunctions away from the left wall; the tail $V(50) = 0.0028$ is small compared with the level spacing (~0.01) but not compared with the solver's precision.

## 4. What the solver produced

| j | $E$ (solver) | $E^{\rm tap}$ (FD) | solver − tap | $E^{\rm raw}$ (FD) | solver − raw |
|---|---|---|---|---|---|
| 1 | 0.0092216090127 | 0.0092216090103 | 2.4e-12 | 0.0092368363 | −1.52e-05 |
| 2 | 0.0193936045680 | 0.0193936045616 | 6.4e-12 | 0.0194255183 | −3.19e-05 |
| 3 | 0.0330092016589 | 0.0330092016331 | 2.6e-11 | 0.0330611573 | −5.20e-05 |
| 4 | 0.0502317685239 | 0.0502317684348 | 8.9e-11 | 0.0503052606 | −7.35e-05 |
| 5 | 0.0711460909545 | 0.0711460906883 | 2.7e-10 | 0.0712406146 | −9.45e-05 |
| 6 | 0.0958063760263 | 0.0958063753557 | 6.7e-10 | 0.0959196717 | −1.13e-04 |
| 8 | 0.1565050043639 | 0.1565050013577 | 3.0e-09 | 0.1566446743 | −1.40e-04 |
| 10 | 0.2325231655953 | 0.2325231550141 | 1.1e-08 | 0.2326732145 | −1.50e-04 |
| 12 | 0.3239720293359 | 0.3239719914102 | 3.8e-08 | 0.3241217743 | −1.50e-04 |

Reading the table:

* **Column "solver − tap"** is the numerical error of the B-spline solve for the problem it actually set up: $10^{-12}$ for the ground state rising to $4\times10^{-8}$ by state 12, positive throughout (variational), growing geometrically as the order-6 basis on $h=1.25$ runs out of resolution ($\lambda_{12} = 2\pi/\sqrt{2\cdot0.32} \approx 7.8$, i.e. 6 nodes per wavelength — comfortable for order 6, hence still $10^{-8}$). **This confirms both the eigensolver and the exact implementation of the window.**
* **Column "solver − raw"** is the effect of the taper on the physics: every level is lowered by $1.5\times10^{-5}$ to $1.5\times10^{-4}$ Ha. That is roughly $\int_{45}^{50}(1-W)V\,|\phi|^2\,dx \approx 0.003 \times (\text{5 bohr}/50 \text{ bohr}) \times (\text{local density factor})$, as expected from first-order perturbation theory, saturating at about $-1.5\times10^{-4}$ once the wavefunction is spread uniformly over the box. Relative to the level spacing it is 0.1–1.5 %; relative to the solver's own precision it is four to seven orders of magnitude.

`figures/case3_irregular_tail.png`: left, the raw and tapered potentials over the last 20 bohr; right, both difference columns on a symmetric-log axis.

## 5. Correctness assessment

| Quantity | Status |
|---|---|
| Case-3 detection ($p=1.5$, $\Delta=4.99$) | correct |
| taper implementation | exact — eigenvalues match an independent solve of $\tilde V$ to $\le 4\times10^{-8}$ |
| eigenvalues vs the untapered box problem | shifted by $-1.5\times10^{-5}$ … $-1.5\times10^{-4}$ (by design) |
| bound-state count (0) | correct |
| eigenfunctions | not separately checked (the eigenvalue agreement to $10^{-11}$ makes them right to $\sim10^{-6}$) |

**Overall: PASS on numerics.** Whether the tapered problem is the problem the user wanted is the open question.

## 6. On the modelling choice

The Case-3 window exists to make the wall region field-free so that the flat-asymptote continuum matching (`matchAsymptotic`, which assumes $V(R)=0$) becomes formally applicable. That is a defensible trade for a continuum calculation: the alternative — matching against WKB or numerically integrated Jost solutions of the tail — is not implemented. But three consequences should be stated in the documentation, because none of them is visible from the warning text:

1. **It changes the bound/box spectrum**, by $O\!\big(V(R)\cdot\Delta/L\big)$. For this test that is $10^{-4}$ Ha. For a potential with a slower tail (say $x^{-1.2}$, which is also Case 3) the shift would be several times larger.
2. **It is applied even when the continuum is disabled**, as in this YAML. In that situation there is no benefit at all.
3. **$\Delta = 0.1\,L$ is arbitrary.** A physically motivated choice is the smallest $\Delta$ over which the taper is adiabatic for the highest continuum energy of interest, $\Delta \gg 1/k_{\max}$; with $k_{\max}=1$ that is a few bohr, and 5 bohr is fine here, but the rule scales with the box rather than with the physics.

## 7. Recommendations

1. **Gate the taper on `tise.continuum.enabled`.** With the continuum off the raw potential should be integrated to the wall; the box eigenvalues are then the eigenvalues of the potential the user wrote. (If the code prefers one consistent Hamiltonian per run, at least say in the warning that the *bound* spectrum is shifted by roughly $\int (1-W)V|\phi|^2$.)
2. **Report the shift.** Cheap and informative: compute $\langle\phi_n|(1-W)V|\phi_n\rangle$ for each state and add it to `warnings.json` (or a new diagnostics file). Users can then see whether the taper matters for their levels; here it would report $1.5\times10^{-5}$ for the ground state.
3. **Make $\Delta$ configurable** (`tise.case3.transition_width`) with the current $0.1L$ as the default, and document the adiabaticity criterion $\Delta \gtrsim$ several $/k_{\max}$.
4. **Strengthen the e2e test.** It asserts finiteness and the warning text. Add `abs(E_1 − 0.00922160901) < 1e-9` (tapered reference) so that a future change to the window shape or width is noticed; `verify_new_tests.py` regenerates the reference in ~10 s.
5. **Longer term:** implement WKB (or numerically integrated) tail matching for Case 3 so that the taper is not needed at all; the SDD's own Figure-7 dispatch anticipates this.
