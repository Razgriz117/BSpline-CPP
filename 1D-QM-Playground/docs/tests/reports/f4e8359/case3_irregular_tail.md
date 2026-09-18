# Test report: `tests/case3_irregular_tail.yaml` — iteration 3

**Verified against:** Richardson-extrapolated finite-difference eigenvalues ($N = 40\,000/80\,000$) of $-\tfrac12\phi'' + x^{-3/2}\phi = E\phi$ on $[0.1, 50]$ with hard walls — the **raw** potential, since the taper is no longer applied.
**Branch / commit:** `TISE-Generalization` @ `f4e8359`. **Previous:** `8236239`.
**Date:** 2026-09-11.
**Verdict:** **PASS.** The Case-3 taper is now gated on `tise.continuum.enabled` (iteration-2 recommendation 1). With the continuum off, the raw potential is integrated to the wall and the twelve lowest eigenvalues agree with the reference for the untapered problem to $1.4\times10^{-13}$–$3.8\times10^{-8}$. At 8236239 the same numbers were shifted by $1.5\times10^{-5}$–$1.5\times10^{-4}$ by an unrequested taper.

Companion files: `verify_new_tests.py`, `figures/case3_irregular_tail.png`.

---

## 1. Input

Unchanged: 41 nodes, order 6, $[0.1, 50]$, $V = x^{-1.5}$ declared on $(0.1,\infty)$, `continuum: enabled: false`, `run_analysis: false`.

## 2. What changed in the solver (commit `0bf17be`)

`classifyAsymptote` still runs (the potential extends beyond the box) and still classifies the tail as Case 3 / Irregular with $\Delta = 4.99$ — the first `warnings.json` entry is the classifier's own message and still says "the potential will be smoothly tapered". But `tise_solver_main.cpp` now only passes `case3RightR/Delta` into `fillBandedMatrices` when the continuum is enabled, and adds a second entry: "continuum is disabled, so the raw (untapered) potential is integrated to the wall — the bound-state spectrum is unaffected."

## 3. The known solution

No closed form; reference as in `../8236239/case3_irregular_tail.md` §3, now compared against $E^{\rm raw}$ only.

## 4. What the solver produced

| j | $E$ (f4e8359) | $E^{\rm raw}$ (FD, Richardson) | solver − raw | at 8236239: solver − raw |
|---|---|---|---|---|
| 1 | 0.0092368363230781 | 0.0092368363232204 | −1.4e-13 | −1.52e-05 |
| 2 | 0.0194255183178811 | 0.0194255183164458 | +1.4e-12 | −3.19e-05 |
| 3 | 0.0330611573027460 | 0.0330611572842069 | +1.9e-11 | −5.20e-05 |
| 4 | 0.0503052607250591 | 0.0503052606445518 | +8.1e-11 | −7.35e-05 |
| 6 | 0.0959196723854836 | 0.0959196717208110 | +6.6e-10 | −1.13e-04 |
| 8 | 0.1566446772612347 | 0.1566446742505112 | +3.0e-09 | −1.40e-04 |
| 10 | 0.2326732251125881 | 0.2326732145189044 | +1.1e-08 | −1.50e-04 |
| 12 | 0.3241218122798638 | 0.3241217743488975 | +3.8e-08 | −1.50e-04 |

The residual column is the same geometric growth ($10^{-13}\to4\times10^{-8}$, positive, variational) seen at 8236239 against the *tapered* reference: the order-6 basis on $h=1.25$ running out of resolution by state 12 (6 nodes per wavelength). The FD reference's own uncertainty is $\sim10^{-11}$, so the first three rows are at the reference's limit, not the solver's. `figures/case3_irregular_tail.png` (right) shows solver − raw at the $10^{-8}$ level versus the $10^{-4}$ offset of the previous iteration.

## 5. Correctness assessment

| Quantity | 8236239 | f4e8359 |
|---|---|---|
| Case-3 detection | correct | correct |
| taper applied with continuum off | yes (unrequested) | **no** |
| eigenvalues vs untapered problem | −1.5e-5 … −1.5e-4 | **1e-13 … 4e-8** |
| warnings | said "will be tapered" only | says both; first message now misleading in isolation |

**Overall: PASS.**

## 6. Recommendations

1. **Make the classifier's warning conditional or reword it.** The first `warnings.json` entry still promises tapering that the second entry then retracts. Either suppress the classifier's remediation sentence when the continuum is off, or have the classifier report detection only ("tail is irregular, p = 1.5") and let the orchestrator state what it does about it.
2. **Add a continuum-enabled twin of this test** (`case3_irregular_tail_continuum.yaml`) so the taper path — now only exercised when the continuum is on — stays covered by a numeric assertion. The tapered reference values from iteration 2 ($E_1^{\rm tap} = 0.0092216090103$) serve directly.
3. Carried from iteration 2: report $\langle\phi_n|(1-W)V|\phi_n\rangle$ per state when the taper *is* applied, make $\Delta$ configurable, and eventually replace the taper by WKB tail matching.
