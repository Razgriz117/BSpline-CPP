# Test report: `tests/finite_square_well.yaml`

**Verified against:** the finite square well with a hard wall at the origin (equivalently, the odd-parity states of the symmetric finite well, or the $\ell = 0$ spherical well), solved by the standard transcendental matching conditions.
**Branch / commit:** `TISE-Generalization` @ `9d0f04b`.
**Date:** 2026-09-04.
**Verdict:** **PASS on physics, MARGINAL on numerics.** The solver finds the right number of bound states with the right structure — including the subtle behaviour of the fifth, marginally bound state — and its continuum machinery converges to the exact phase shifts. On the shipped 51-node grid, however, bound energies are only accurate to $2\times10^{-4}$–$6\times10^{-3}$, the lowest positive-energy box state is misplaced by a factor of ten (and on the wrong side of the empty-box ladder), and the continuum at $E = 0.5$ is visibly under-resolved. Both are grid/basis limitations that disappear on refinement, not solver errors. The `dδ/dE` column is not usable.

Companion files: `known_solutions_report.html`, `verify_known_solutions.py`, `figures/`.

---

## 1. Input

```yaml
run:        output_dir: "./data/finite_square_well"
bspline:    n_nodes: 51,  order: 12,  domain: [0.0, 100.0]
potential:  "{'domain': '[0, 10)',  'function': '-1.0'}"
            "{'domain': '[10, 100]', 'function': '0'}"
tise:       n_pts_eigenstate: 301, error_threshold: 1.0e-10
            continuum: enabled: true, E_threshold: 0.0, E_max: 0.5, n_energies: 5, n_pts: 500
visualization: eigenstates: true
```

Well depth $V_0 = 1$, width $a = 10$, box length 100, node spacing $h = 2$. The potential step at $x = 10$ falls exactly on a node ($10 = 5h$), so Gauss–Legendre quadrature never straddles the discontinuity — good — but see §5 for why the step still hurts.

## 2. What the solver actually poses

The dropped edge B-splines impose $\phi(0) = \phi(100) = 0$. The wall at $x = 0$ changes which textbook problem this is: it is **not** the symmetric well on $[-a, a]$, it is the right half of it with $\phi(0) = 0$. Only the odd states of the symmetric well vanish at the origin, so this configuration has exactly the odd-parity spectrum of a symmetric well of half-width 10 and depth 1. (Equivalently, it is the $\ell = 0$ radial equation for a spherical square well, which is why the plots look like a nucleus.) The wall at $x = 100$ additionally converts the $E > 0$ continuum into discrete box states and slightly perturbs any bound state whose tail reaches it.

## 3. The known solution, step by step

### 3.1 Bound states, $-V_0 < E < 0$

**Step 1 — inside the well** ($x < a$), the particle has positive kinetic energy $E + V_0$ and oscillates:

$$\phi_{\rm in}(x) = A\sin Kx, \qquad K = \sqrt{2(E + V_0)}.$$

The sine rather than the cosine because $\phi(0) = 0$.

**Step 2 — outside** ($x > a$), $E < 0$ and the wave must decay:

$$\phi_{\rm out}(x) = B e^{-\kappa x}, \qquad \kappa = \sqrt{-2E}.$$

Picture the wave sloshing inside the well and leaking a decaying tail through the wall.

**Step 3 — glue at $x = a$.** Continuity of $\phi$ and $\phi'$; dividing the two conditions eliminates $A, B$:

$$K\cot(Ka) = -\kappa.$$

Both sides depend on $E$ alone; a bound state is a crossing. $\cot(Ka)$ sweeps from $+\infty$ to $-\infty$ on each interval $Ka \in \big((m-\tfrac12)\pi,\ m\pi\big)$, while $-\kappa$ is a smooth negative curve, so there is exactly one crossing in each such interval that $Ka$ can reach before $E$ hits 0.

**Step 4 — count the states.** At $E \to 0^-$, $Ka \to \sqrt{2V_0}\,a = \sqrt{200} = 14.1421$. Compare $4.5\pi = 14.1372$. The fifth threshold is exceeded by only 0.005, so on an infinite domain this well has **five** bound states, the fifth extraordinarily shallow.

**Step 5 — solve numerically** (Brent's method in each bracket):

| m | $E_m$ (infinite domain) | $\kappa$ | decay length $1/\kappa$ |
|---|---|---|---|
| 1 | −0.956996254 | 1.383 | 0.72 |
| 2 | −0.828522650 | 1.287 | 0.78 |
| 3 | −0.616571751 | 1.110 | 0.90 |
| 4 | −0.326731918 | 0.808 | 1.24 |
| 5 | −0.000023092 | 0.0068 | **147** |

**Step 6 — account for the wall at $x = 100$.** A tail with decay length 147 bohr cannot fit in the 90-bohr corridor between the well and the wall. With a hard wall the outside solution is $\sinh\kappa(100 - x)$ instead of $e^{-\kappa x}$ and the matching condition becomes

$$K\cot(Ka) = -\kappa\coth\!\big(\kappa(100-a)\big).$$

For $m = 1$–4 this changes nothing at the $10^{-9}$ level ($\coth \to 1$). For $m = 5$ the equation has **no root** — checked on a 200 000-point scan of the bracket, the left-hand side minus right-hand side stays positive throughout. The wall pushes the fifth state up into the positive-energy box spectrum. **So the correct answer for this input is four bound states**, followed by box states.

**Step 7 — where the box states sit.** For $E > 0$ the outside wave is $\sin(kx + \delta(E))$ with $\delta$ from §3.2, and the wall requires $\sin(100k + \delta) = 0$. Solving $100k + \delta(E) = j\pi$ gives the exact box levels

$$E = 5.53\times10^{-5},\ 1.164\times10^{-3},\ 3.374\times10^{-3},\ 6.709\times10^{-3},\ \dots$$

all *below* the empty-box ladder $j^2\pi^2/(2\cdot100^2) = 4.93\times10^{-4},\ 1.97\times10^{-3},\ 4.44\times10^{-3},\ 7.90\times10^{-3}$, as they must be for an attractive well — and the lowest one, at $5.5\times10^{-5}$, is the ghost of the marginally bound fifth state.

### 3.2 Continuum, $E > 0$

Both regions oscillate: $\sin Kx$ inside and $\sin(kx + \delta)$ outside, $k = \sqrt{2E}$. Matching the logarithmic derivative at $x = a$, $K\cot Ka = k\cot(ka + \delta)$, gives the standard result

$$\delta(E) = \arctan\!\Big(\frac{k}{K}\tan Ka\Big) - ka \pmod\pi.$$

Two features to see in any continuum plot: the inside wave has a shorter wavelength (higher local momentum) and a smaller amplitude (probability current $\propto K|A|^2$ must equal $k|B|^2$), and the outside wave is displaced relative to a free wave — that displacement is $\delta$. The exact $d\delta/dE$ was obtained by differentiating this formula numerically with $h = 10^{-5}$.

## 4. What the solver produced

### 4.1 Bound energies — and their convergence with node count

| n_nodes (h) | $E_1$ | $E_2$ | $E_3$ | $E_4$ |
|---|---|---|---|---|
| 51 (2.0) — *shipped* | −0.9567575 (+2.4e-4) | −0.8274903 (+1.0e-3) | −0.6138839 (+2.7e-3) | −0.3203209 (+6.4e-3) |
| 101 (1.0) | −0.9569401 (+5.6e-5) | −0.8282904 (+2.3e-4) | −0.6160201 (+5.5e-4) | −0.3256849 (+1.0e-3) |
| 201 (0.5) | −0.9569881 (+8.1e-6) | −0.8284900 (+3.3e-5) | −0.6164983 (+7.3e-5) | −0.3266035 (+1.3e-4) |
| 401 (0.25) | −0.9569952 (+1.0e-6) | −0.8285185 (+4.1e-6) | −0.6165625 (+9.2e-6) | −0.3267161 (+1.6e-5) |
| **analytic** | **−0.9569963** | **−0.8285227** | **−0.6165718** | **−0.3267319** |

(Parentheses: error relative to the analytic value. The 101/201/401 rows were produced by editing only `n_nodes` in the YAML.)

Observations:

* Exactly four negative eigenvalues on the shipped grid, precisely as §3.1 step 6 predicts. The count and the fifth-state behaviour are **correct**.
* The first box states, however, are badly placed on the shipped grid. Solver (51 nodes): $E_5 = +5.57\times10^{-4}$, $E_6 = 2.22\times10^{-3}$, $E_7 = 4.96\times10^{-3}$, $E_8 = 8.71\times10^{-3}$. Exact (§3.1 step 7): $5.53\times10^{-5}$, $1.164\times10^{-3}$, $3.374\times10^{-3}$, $6.709\times10^{-3}$. The shipped values sit *above* the empty-box ladder, the exact ones *below* it — a qualitative error in the near-threshold spectrum, with the lowest box state off by a factor of ten. Refinement repairs it: 101 nodes gives $2.5\times10^{-4}$, 201 gives $8.3\times10^{-5}$, 401 gives $5.86\times10^{-5}$, $1.167\times10^{-3}$, $3.376\times10^{-3}$, $6.711\times10^{-3}$. The mechanism is the same $\phi''$-kink defect as for the bound states (§5), but its consequence is larger here because the ghost of the fifth state is exquisitely sensitive to the matching at $x = 10$.
* Every solver energy lies *above* the exact one. The method is variational (Rayleigh–Ritz in a B-spline subspace), so it can never undershoot; this is a consistency check that passes.
* Each halving of $h$ reduces the error by a factor rising from ~4 toward ~8: algebraic convergence like $h^3$, in stark contrast to the $10^{-13}$ the same grid achieves for the free particle. See §5 for the cause.

### 4.2 Eigenfunctions

Compared to the normalised piecewise $\sin Kx\,/\,e^{-\kappa(x-a)}$ solution built from the analytic $E_m$:

| n | $\int\phi_n^2\,dx$ | rel. $L^2$ diff (51 nodes) |
|---|---|---|
| 1 | 0.9999989 | 0.9 % |
| 2 | 1.0000051 | 2.0 % |
| 3 | 0.9999860 | 3.5 % |
| 4 | 1.0000387 | 6.2 % |
| 5 (would-be) | — | 49 % (not a bound state in the box; comparison is meaningless) |

Node counts (0, 1, 2, 3 interior zeros) and the shape — sinusoid inside, exponential shoulder just outside $x = 10$ — are right (`figures/output_finite_square_well.png`, `figures/fsw_states.png`). The visible discrepancy is a small lag near $x = 10$ where the spline cannot bend sharply enough. The 4-in-10⁻⁵ deviations of the norm from 1 are a symptom of the same thing: the eigenvector is $S$-normalised exactly, but the trapezoid integral of the output samples is not, because the output grid (301 points, spacing 0.33) samples a function with a poorly represented kink.

### 4.3 Continuum and phase shifts

| E | $\delta \bmod \pi$ (solver, 51) | $\delta \bmod \pi$ (solver, 401) | $\delta \bmod \pi$ (exact) | $d\delta/dE$ (solver) | $d\delta/dE$ (exact) | rel. $L^2$ err (51) |
|---|---|---|---|---|---|---|
| 0.1 | 1.3543 | 1.4641 | **1.4642** | +4.06 | **−19.43** | 11.5 % |
| 0.2 | −0.1617 | −0.1307 | **−0.1307** | −2.45 | **−13.26** | 3.4 % |
| 0.3 | −1.3037 | −1.2534 | **−1.2534** | −3.59 | **−9.24** | 5.5 % |
| 0.4 | 0.9896 | 1.2024 | **1.2025** | −8.09 | **−4.06** | 21.8 % |
| 0.5 | 0.4442 | 1.0675 | **1.0678** | −1.12 | **−0.08** | 59.5 % |

* On the shipped grid the phase shifts are right to 0.03–0.1 rad for $E \le 0.3$ and wrong by 0.2 and 0.6 rad at $E = 0.4, 0.5$. At 401 nodes all five agree with the exact values to $\le 5\times10^{-4}$: **the continuum construction and the asymptotic matching are correct**, the shipped resolution is not adequate.
* Why $E = 0.5$ is so bad: inside the well the local energy is $E + V_0 = 1.5$, $K = 1.73$, wavelength 3.6 bohr, and $h = 2$ provides **1.8 nodes per wavelength** — below even the Nyquist limit. The solver's own accuracy check, $E_{\rm acc} = \pi^2/2h^2 = 1.23$, did not fire because it compares $E_{\max} = 0.5$ to the grid without subtracting $V_{\min} = -1$.
* `figures/fsw_continuum.png` shows both cases: at $E = 0.1$ the outside phase and the reduced inside amplitude are reproduced; at $E = 0.5$ the inside oscillation is under-sampled and the outside phase is off.
* The $d\delta/dE$ column is unrelated to the exact derivative at every energy. This is independent of resolution: `matchAsymptotic` computes a central difference of $\sin 2\delta$ over $\Delta E = 0.1$, then divides by $2\cos 2\delta$. Between neighbouring grid points $\delta$ changes by 1–2 rad, so the difference quotient samples a function it cannot see. The exact derivative is large and negative at low energy (−19 at $E = 0.1$) and rises smoothly toward zero — the familiar slope of a square-well $\delta(E)$ curve.

## 5. Why this test converges slowly: the step in $V$

The exact eigenfunction is continuous with continuous first derivative at $x = 10$, but its **second** derivative jumps by $2V_0\phi(10)$ (read directly off the TISE: $\phi'' = 2(V - E)\phi$). An order-12 B-spline basis on *simple* knots is $C^{10}$-continuous across every node, so it cannot represent that jump; it can only approximate it by bending over the adjacent intervals. The energy error from a function-space defect of this kind scales algebraically with $h$ — the observed $h^3$ — and no amount of polynomial order helps.

The standard cure, described in the B-spline notes and by Bachau et al., is to **repeat the knot at $x = 10$** with multiplicity $k - 2 = 10$ (or simply $k - 1$), which reduces the basis to $C^1$ (or $C^0$) exactly there and nowhere else. The branch already contains this machinery — commit `387f384` "wire strategic node placement and singular-join B-spline removal into `solveTISE`" — but `tise_solver_main.cpp` always calls `buildUniformRadialGrid` and never uses it.

## 6. Correctness assessment

| Quantity | Status |
|---|---|
| Number of bound states in the box (4) and fifth-state behaviour | correct |
| Bound energies (51 nodes) | correct to 2.4e-4 … 6.4e-3; converge as $h^3$ to the exact values |
| Bound eigenfunctions | correct shape/nodes; 1–6 % $L^2$ error at 51 nodes, concentrated at $x = 10$ |
| Box states above threshold (51 nodes) | **qualitatively wrong** — lowest one 5.6e-4 vs exact 5.5e-5, all above instead of below the empty-box ladder; converges on refinement |
| Continuum $\psi_E$, $E \le 0.3$ | correct to 3–11 % at 51 nodes; converges |
| Continuum $\psi_E$, $E = 0.4, 0.5$ | under-resolved at 51 nodes (1.8 nodes/λ inside the well) |
| Phase shifts $\delta \bmod \pi$ | correct once resolved (401 nodes: ≤ 5e-4 rad) |
| $d\delta/dE$ | **wrong** at every energy; method defect, not resolution |

**Overall: the physics is right and the code is right; the test configuration is too coarse to show it, and one derived quantity is computed by an unsound method.**

## 7. Recommendations

1. **Refine the shipped grid.** `n_nodes: 201` (h = 0.5) brings the bound energies to $10^{-5}$–$10^{-4}$ and the phase shifts to $10^{-3}$ rad at ~25 s run time; `401` gives $10^{-6}$ / $5\times10^{-4}$. Alternatively keep 51 nodes and lower `E_max` to 0.3 so the test only claims what it can deliver.
2. **Expose knot multiplicity at potential steps in `tise_solver`.** Add a config field (e.g. `bspline.breakpoints: [10.0]` or `bspline.knot_multiplicity_at_pieces: true`) that routes the piece boundaries of the `potential` list into the strategic-node-placement path that `solveTISE` already has. This is the change that would turn the $h^3$ convergence back into spectral convergence and let 51 nodes deliver $10^{-10}$.
3. **Fix $E_{\rm acc}$.** Compute it from the *local* kinetic energy: warn when $E_{\max} - V_{\min} > E_{\rm acc}$ (and consider a safety factor: the oscillator and free-particle data indicate ~3.3 nodes per wavelength for $10^{-6}$, i.e. warn at about $E_{\rm acc}/2.7$).
4. **Replace the $d\delta/dE$ estimator.** Options in increasing order of effort: (a) evaluate $\delta$ on a fine local energy grid ($\Delta E \sim 10^{-3}$) around each requested energy and difference that, after unwrapping the branch; (b) differentiate the matching formula analytically — $\delta = \arctan(k\psi/\psi') - kR$ requires $\partial_E\psi(R)$ and $\partial_E\psi'(R)$, which are available in closed form from the same $(E - E_n)^{-1}$ expansion that builds $\psi_E$; (c) drop the column until (a) or (b) is in place. The present $\sin 2\delta / \cos 2\delta$ formulation should go regardless: it divides by $\cos 2\delta$, which vanishes at $\delta = \pi/4$.
5. **Unwrap $\delta(E)$** across the energy grid so the stored phase shift is continuous rather than an arbitrary multiple of $\pi$ below zero; this is a prerequisite for 4(a).
6. **Output-grid normalisation check.** Consider having `analysis.py` report $\int\phi_n^2$ on the output grid; a deviation from 1 of more than ~$10^{-6}$ is a cheap resolution alarm, as it was here.
7. **Automated regression.** Assert four negative eigenvalues; assert $E_1$–$E_4$ within a tolerance appropriate to the chosen grid (e.g. $10^{-4}$ at 201 nodes) of the analytic roots of $K\cot Ka = -\kappa$; assert $\delta \bmod \pi$ within $10^{-2}$ of the closed form. All three are implemented in `verify_known_solutions.py` (with the 51-node tolerances it reports rather than asserts).
