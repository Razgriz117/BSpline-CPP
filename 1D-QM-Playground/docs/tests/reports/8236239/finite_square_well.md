# Test report: `tests/finite_square_well.yaml` — iteration 2

**Verified against:** the finite square well with a hard wall at the origin (odd-parity states of the symmetric well / $\ell=0$ spherical well), transcendental matching conditions.
**Branch / commit:** `TISE-Generalization` @ `8236239` ("TISE release-readiness pass — unify grid/dropset, wire diagnostics, add CI").
**Previous iteration:** `9d0f04b` — see `../9d0f04b/finite_square_well.md`.
**Date:** 2026-09-04.
**Verdict:** **PASS.** The round-1 defect is fixed. With the strategic knot at the potential step now reaching `tise_solver`, the four bound energies agree with the exact roots of $K\cot Ka=-\kappa$ to $10^{-13}$–$10^{-9}$ (previously $10^{-4}$–$10^{-2}$), the near-threshold box states are on the correct side of the empty-box ladder and accurate to $10^{-5}$ (previously off by 10×), and the continuum phase shifts match the closed form to $2\times10^{-4}$ rad at every energy (previously up to 0.6 rad). Two round-1 issues remain: the $d\delta/dE$ column is still unusable, and the new `phase_shifts.png` plots the unwrapped $\delta \approx -kR$, which hides the physics.

Companion files: `verify_known_solutions.py`, `figures/`.

---

## 1. Input

```yaml
run:        output_dir: "./data/finite_square_well"
bspline:    n_nodes: 51,  order: 12,  domain: [0.0, 100.0]
potential:  "{'domain': '[0, 10)',  'function': '-1.0'}"
            "{'domain': '[10, 100]', 'function': '0'}"
tise:       n_pts_eigenstate: 301, error_threshold: 1.0e-10
            continuum: enabled: true, E_threshold: 0.0, E_max: 0.5, n_energies: 5, n_pts: 500
visualization: eigenstates: true, phase_shifts: true      # phase_shifts is new at 8236239
```

Unchanged physics: $V_0=1$, $a=10$, box 100, base node spacing $h=2$.

## 2. What changed in the solver between 9d0f04b and 8236239

`tise_solver_main.cpp` no longer builds a plain uniform grid. It calls the new shared `buildStrategicGridAndDropSet` (ADR-0009), which runs `detectPotentialStructure` on the piecewise potential, classifies the join at $x=10$ as a **Step**, and inserts $\text{order}-3 = 9$ extra copies of the knot at $x=10$. Total multiplicity 10 at that knot means the basis is only $C^{12-10-1}=C^1$ there — exactly the smoothness the true eigenfunction has (continuous $\phi$, $\phi'$; jump in $\phi''$). The basis grows from 61 to 70 B-splines, so `eigenvalues.dat` now has 68 rows instead of 59. This is precisely recommendation 2 of the round-1 report.

Also new and visible in this test's output: `warnings.json` reports "4 of 68 computed states are below E=0.0" (`classifyBoundStates`), `analysis.py` writes `phase_shifts.png`, and `computeEAcc` now uses `minInterNodeGap` (still 2.0 here).

## 3. The known solution

The derivation is unchanged from round 1 (`../9d0f04b/finite_square_well.md` §3) and is summarised here for completeness.

**Bound states.** With $\phi(0)=0$: $\phi=A\sin Kx$ inside ($K=\sqrt{2(E+V_0)}$), $Be^{-\kappa x}$ outside ($\kappa=\sqrt{-2E}$); matching gives

$$K\cot(Ka)=-\kappa,$$

one root per interval $Ka\in((m-\tfrac12)\pi, m\pi)$. Since $\sqrt{2V_0}\,a = 14.1421 > 4.5\pi = 14.1372$ there are five roots on an infinite domain: $-0.956996254,\ -0.828522650,\ -0.616571751,\ -0.326731918,\ -0.000023092$. The fifth has decay length $1/\kappa = 147$ bohr; with the wall at 100 the matching condition becomes $K\cot Ka = -\kappa\coth\kappa(100-a)$, which has **no** root in the fifth bracket. **Four bound states** is the exact answer for this input.

**Box states above threshold.** $100k+\delta(E)=j\pi$ with $\delta$ below gives $5.5253\times10^{-5},\ 1.16434\times10^{-3},\ 3.37370\times10^{-3},\ 6.70885\times10^{-3},\dots$ — below the empty-box ladder, the first being the ghost of the fifth bound state.

**Continuum.** $\delta(E)=\arctan\!\big(\tfrac{k}{K}\tan Ka\big)-ka \pmod\pi$; $d\delta/dE$ by differentiating this.

## 4. What the solver produced

### 4.1 Bound energies

| n | $E_n$ (8236239) | exact | rel. error | rel. error at 9d0f04b |
|---|---|---|---|---|
| 1 | −0.956996254368899 | −0.956996254368992 | 9.7e-14 | 2.5e-04 |
| 2 | −0.828522650434678 | −0.828522650434757 | 9.5e-14 | 1.2e-03 |
| 3 | −0.616571750670465 | −0.616571750671813 | 2.2e-12 | 4.4e-03 |
| 4 | −0.326731917714706 | −0.326731918129711 | 1.3e-09 | 2.0e-02 |

The improvement is nine to ten orders of magnitude on the same 51-node grid — the signature of restoring spectral convergence by giving the basis the right smoothness class. The remaining $10^{-9}$ on $E_4$ is ordinary resolution (its $\kappa = 0.81$ tail and $K=1.16$ interior on $h=2$).

### 4.2 Box states above threshold

| j | $E$ (8236239) | exact ($100k+\delta=j\pi$) | rel. error | at 9d0f04b |
|---|---|---|---|---|
| 1 | 5.52701e-05 | 5.52530e-05 | 3.1e-04 | 5.57e-04 (10× too high, wrong side of ladder) |
| 2 | 1.164353e-03 | 1.164340e-03 | 1.1e-05 | 2.22e-03 |
| 3 | 3.373712e-03 | 3.373700e-03 | 3.6e-06 | 4.96e-03 |
| 4 | 6.708862e-03 | 6.708850e-03 | 1.7e-06 | 8.71e-03 |

The near-threshold spectrum — the most delicate part of this problem — is now right, including the ghost state at $5.5\times10^{-5}$.

### 4.3 Eigenfunctions

| n | rel. $L^2$ diff (8236239) | at 9d0f04b |
|---|---|---|
| 1 | 1.5e-06 | 0.9 % |
| 2 | 5.5e-06 | 2.0 % |
| 3 | 1.1e-05 | 3.5 % |
| 4 | 1.6e-05 | 6.2 % |

`figures/fsw_states.png`: solver and analytic curves are indistinguishable through the well edge at $x=10$, where round 1 showed a visible lag. The fifth panel now shows the solver's $E=+5.5\times10^{-5}$ box state next to the would-be marginal bound state; they are the same object squeezed by the wall. (The trapezoid norm on the 301-point output grid is $1-3\times10^{-5}$ rather than $1-10^{-6}$: the output sampler now sees a genuine $\phi''$ kink at 10, which the 0.33-bohr trapezoid rule integrates less accurately. The $S$-norm is exact; this is a diagnostic-of-the-diagnostic, not an error.)

### 4.4 Continuum and phase shifts

| E | $\delta \bmod\pi$ (8236239) | exact | diff | diff at 9d0f04b | $d\delta/dE$ (solver) | $d\delta/dE$ (exact) | rel. $L^2$ err vs analytic $\psi_E$ |
|---|---|---|---|---|---|---|---|
| 0.1 | 1.46421 | 1.46421 | 2e-06 | −0.110 | +2.41 | −19.43 | 5.1e-05 |
| 0.2 | −0.13073 | −0.13072 | −1e-05 | −0.031 | −2.08 | −13.26 | 9.0e-05 |
| 0.3 | −1.25338 | −1.25335 | −3e-05 | −0.050 | −2.89 | −9.24 | 1.8e-04 |
| 0.4 | 1.20242 | 1.20249 | −7e-05 | −0.213 | −4.85 | −4.06 | 3.7e-04 |
| 0.5 | 1.06756 | 1.06780 | −2.4e-04 | −0.624 | −1.62 | −0.08 | 6.0e-04 |

Two things to record honestly:

* **The phase shifts are now correct to $\le 2.4\times10^{-4}$ rad at all five energies on the shipped grid.** Round 1 attributed the 0.6 rad error at $E=0.5$ partly to node density (1.8 nodes per wavelength inside the well). This result shows the smoothness defect at $x=10$ was the dominant cause; with it removed, $h=2$ is adequate to $10^{-4}$ even at $E=0.5$. The round-1 recommendation to raise `n_nodes` is therefore no longer necessary for this test.
* **$d\delta/dE$ is still wrong at every energy**, and still for the same reason: a central difference of $\sin2\delta$ over $\Delta E = 0.1$, divided by $2\cos2\delta$. The release-readiness plan's Part C investigated the *sign* at the last grid point and concluded (correctly) that the formula has no sign bug — but the formula is only valid when $\delta$ changes by $\ll 1$ rad between grid points, and here it changes by 1–2 rad. The numbers in this column (+2.4 where the exact slope is −19.4) are not approximations of anything.

`figures/fsw_continuum.png`: at $E=0.5$ the round-1 mismatch is gone; solver and analytic curves coincide inside and outside the well.

### 4.5 The new `phase_shifts.png`

`figures/output_phase_shifts_fp_fsw.png` (right panel) shows what `analysis.py` now draws: the upper panel is the *raw* $\delta$ from `phase_shifts.dat`, i.e. $\arctan(\cdot)-kR$, which runs from −46 to −100. It is a straight line of slope $-R/k$ with the physics ($\pm 1$ rad of structure) invisible on it. The lower panel plots the unusable derivative. As it stands the figure conveys nothing about the well.

## 5. Correctness assessment

| Quantity | 9d0f04b | 8236239 |
|---|---|---|
| number of bound states (4) | correct | correct |
| bound energies | 2e-4 … 2e-2 rel. | **1e-13 … 1e-9** |
| box states near threshold | qualitatively wrong | **correct to 3e-4 … 2e-6** |
| bound eigenfunctions | 1–6 % | **1e-6 … 2e-5** |
| continuum $\psi_E$ | 3–60 % | **5e-5 … 6e-4** |
| $\delta \bmod \pi$ | 0.03–0.6 rad | **≤ 2.4e-4 rad** |
| $d\delta/dE$ | unusable | unusable |
| `phase_shifts.png` | — | plots raw $-kR$-dominated $\delta$; uninformative |

**Overall: PASS on every physical quantity the test is meant to certify; the derivative column and the new plot are the open items.**

## 6. Recommendations

1. **Close the $d\delta/dE$ item properly.** The Part-C test (`MatchAsymptoticDDeltaDESignTest`) proves the formula on a "narrow pole-avoiding grid", i.e. a grid fine enough that $\Delta\delta \ll 1$. Production runs use `n_energies: 5` over 0.5 Ha. Either (a) compute $d\delta/dE$ from an internal fine grid ($\Delta E \sim 10^{-3}$) around each requested energy after unwrapping $\delta$; (b) differentiate the matching condition analytically ($\partial_E\psi(R)$ and $\partial_E\psi'(R)$ are available from the $(E-E_n)^{-1}$ expansion that builds $\psi_E$); or (c) omit the column and the lower panel of the plot. The $\sin2\delta/\cos2\delta$ construction should go in any case — it divides by zero at $\delta=\pi/4$ and is the reason hydrogen's value at $E=0.1$ is +211.
2. **Store and plot a physically meaningful $\delta$.** Unwrap the branch across the energy grid so that $\delta(E)$ is continuous and starts in $(-\pi/2,\pi/2]$, and drop the $-kR$ term into a separate column if it is wanted. Then `phase_shifts.png` will show what this well actually does: a monotone decrease of about 3.5 rad, from 1.46 at $E=0.1$ through 0 near $E=0.19$ to −2.07 at $E=0.5$ (the exact continuous values are 1.464, −0.131, −1.253, −1.939, −2.074).
3. **Keep the strategic-knot path under regression test with a tolerance that would have caught round 1.** The current e2e test for this YAML asserts "4 negative eigenvalues" and warnings content; add `abs(E_1 + 0.956996254) < 1e-9` and `abs(delta_mod_pi(E=0.5) − 1.06780) < 1e-3`. Both numbers are in `verify_known_solutions.py`.
4. **E_acc rule.** Unchanged from round 1: `computeEAcc` should subtract $V_{\min}$ and warn at a fraction of the Nyquist value; it stays silent here even though the interior local energy at $E=0.5$ is 1.5.
