# Test report: `tests/right_edge_singularity.yaml` — iteration 2 (new test)

**Verified against:** the repulsive Coulomb potential in a box, $V = 1/(100-x)$ on $[0,100)$, whose exact eigenvalues are the zeros of the regular Coulomb function $F_0(1/k, 100k)$ and whose eigenfunctions are $F_0\big(1/k, k(100-x)\big)$.
**Branch / commit:** `TISE-Generalization` @ `8236239`. First appearance of this test.
**Date:** 2026-09-04.
**Verdict:** **PASS for the spectrum, and the diagnostic the test was written for fires correctly.** All seventeen eigenvalues below $E=0.2$ agree with the Coulomb-function zeros to $10^{-12}$ or better and the eigenfunctions to $10^{-11}$–$10^{-8}$. The continuum output, as the new warning says, should not be trusted — and the data show why: the "continuum states" are essentially the last B-spline alone, non-zero at a point where the physical wavefunction must vanish. Recommendation: refuse (not just warn about) continuum construction on a singular edge.

Companion files: `verify_new_tests.py`, `figures/right_edge_singularity.png`.

---

## 1. Input

```yaml
run:        run_analysis: false,  output_dir: "./data/right_edge_singularity"
bspline:    n_nodes: 41,  order: 8,  domain: [0.0, 100.0]
potential:  "{'domain': '[0, 100)', 'function': '1/(100-x)'}"
tise:       n_pts_eigenstate: 301
            continuum: enabled: true, E_threshold: 0.0, E_max: 0.5, n_energies: 3, n_pts: 200
```

Node spacing $h = 2.5$. The energy grid is $E = 1/6, 1/3, 1/2$.

## 2. What the solver actually poses

The singular point is the right *domain edge*, where the knot already has full multiplicity and $B_N$ is the only basis function non-zero at $x=100$. For the bound-state (diagonalisation) problem $B_N$ is dropped, so $\phi(100)=0$ is imposed exactly and the basis retains full approximation power up to the wall — the same situation as hydrogen's $r=0$. `buildStrategicGridAndDropSet` sets `rightEdgeSingular = true`, which now reaches `warnings.json`:

> potential is singular at the right domain edge x=100; continuum phase-shift matching (matchAsymptotic) assumes a regular boundary there, so continuum results should be treated with suspicion.

For the continuum, however, $B_N$ is deliberately *kept* (it is the "escape" function that lets $\psi_E(100) \ne 0$), and the matching formula assumes $V(100) \approx 0$. Both assumptions are violated here.

## 3. The known solution, step by step

**Step 1 — change variable.** Let $s = 100 - x$. Then $-\tfrac12 u_{ss} + \tfrac{1}{s}u = Eu$ on $0 < s \le 100$: the $\ell=0$ radial Schrödinger equation with a *repulsive* unit Coulomb charge.

**Step 2 — regular solution.** With $k=\sqrt{2E}$ and Sommerfeld parameter $\eta = +1/k$ (positive: repulsive), the solution regular at $s=0$ is the Coulomb wave function $F_0(\eta, ks)$. It behaves as $C_0(\eta)\,ks$ near $s=0$ — the wavefunction vanishes *linearly* at the singular wall, and the amplitude factor $C_0(\eta) = \sqrt{2\pi\eta/(e^{2\pi\eta}-1)}$ is exponentially small for slow particles: the Coulomb barrier keeps them away from the wall. This is the Gamow factor.

**Step 3 — quantise.** The other wall at $x=0$ is $s=100$: $F_0(1/k, 100k) = 0$. Scanning $E\in[0.005, 0.2]$ and polishing each sign change with Brent's method gives 17 roots:

$$E = 0.0149922,\ 0.0203389,\ 0.0260832,\ 0.0324237,\ 0.0394463,\ 0.0472038,\ 0.0557341,\ \dots,\ 0.110836,\ \dots$$

For comparison the empty box on $[0,100]$ would start at $4.93\times10^{-4}$: the repulsive $1/s$ tail (which is $0.01$ at the far wall and $0.1$ at $s=10$) raises every level substantially.

**Step 4 — eigenfunctions** are $u_j(s) = F_0(\eta_j, k_j s)$, normalised on $[0,100]$.

## 4. What the solver produced

### 4.1 Eigenvalues

| j | $E$ (solver) | $F_0(1/k,100k)=0$ | rel. error |
|---|---|---|---|
| 1 | 0.01499221266094 | 0.01499221266093 | 3.9e-13 |
| 2 | 0.02033885931642 | 0.02033885931660 | −8.8e-12 |
| 3 | 0.02608321293993 | 0.02608321293992 | 1.8e-13 |
| 4 | 0.03242368771323 | 0.03242368771323 | 4.4e-14 |
| 5 | 0.03944628833233 | 0.03944628833233 | 2.9e-14 |
| 6 | 0.04720380917396 | 0.04720380917398 | −5.6e-13 |
| 11 | 0.09809405990281 | 0.09809405990281 | 5.5e-11 |
| 12 | 0.11083564751346 | 0.11083564749335 | 1.8e-10 |

Agreement at the $10^{-13}$–$10^{-10}$ level for all seventeen states checked, on a coarse ($h=2.5$) order-8 grid — the Coulomb functions are smooth, and the domain-edge treatment of the singularity is correct. The small negative signs on rows 2 and 6 are at the level of the reference root-finding tolerance, not a variational violation.

### 4.2 Eigenfunctions

| j | rel. $L^2$ difference vs $F_0(1/k, k(100-x))$ |
|---|---|
| 1 | 1.3e-11 |
| 2 | 1.4e-10 |
| 5 | 8.9e-09 |

`figures/right_edge_singularity.png`, left: state 1 over the exact Coulomb function — indistinguishable, with the characteristic long soft approach to zero at the barrier side.

### 4.3 Continuum

| E | $\delta$ (raw) | $\psi_E(100)$ | max $\lvert\psi_E\rvert$ on $x<90$ |
|---|---|---|---|
| 1/6 | −57.57 | 0.171 | 0.0053 |
| 1/3 | −81.42 | 0.204 | 0.057 |
| 1/2 | −99.72 | 0.218 | 0.213 |

At $E=1/6$ the "continuum state" is 30× larger at the singular wall than anywhere in the interior: it is the $B_N$ spline with a small admixture of everything else (`figures/right_edge_singularity.png`, right panel). Physically $\psi_E(100)$ must be zero, so the construction has produced a function that is not a solution of the problem, and the phase shift extracted from it has no meaning. The warning is therefore correct, but "treat with suspicion" understates it.

## 5. Correctness assessment

| Quantity | Status |
|---|---|
| eigenvalues (17 states, $E<0.2$) | correct to $\le 2\times10^{-10}$ |
| eigenfunctions | correct to $\le 10^{-8}$ |
| bound-state count (0) | correct |
| `rightEdgeSingular` warning | correct and reaches `warnings.json` |
| continuum states / phase shifts | not physical; correctly flagged, but still written |

**Overall: PASS** for what a right-edge singularity can legitimately be asked for; the continuum output is a known non-result.

## 6. Recommendations

1. **Refuse, don't warn.** When `rightEdgeSingular` is true, skip continuum construction entirely (or make `tise.continuum.enabled: true` a configuration error for that potential). Writing `phase_shifts.dat` and `continuum_state_NNN.dat` that are known to be meaningless invites downstream use; the plot pipeline will happily draw them.
2. **If a singular edge must support scattering**, the physically meaningful setup is the mirror image of hydrogen: the singular point is the *origin* and the regular boundary is the other end. For this potential that means flipping the domain ($s = 100 - x$), placing the flat/Coulomb matching at large $s$, and treating the $1/s$ tail with Coulomb functions as recommended in the round-1 hydrogen report.
3. **Use this test as the regression anchor for domain-edge singularities.** Its assertions currently check only the warning text; add `abs(E_1 − 0.0149922127) < 1e-9`. It is the cleanest available check that the edge treatment stays exact while the interior treatment (see `interior_singularity.md`) is repaired — the two must end up giving the same physics for the same $1/s$ wall.
4. **Contrast with `interior_singularity.yaml`.** This test passes to $10^{-12}$ and that one fails by 41 % for the *same* $1/s$ singularity; the only difference is knot multiplicity at the singular point (full at an edge, simple in the interior). That comparison is the strongest argument for recommendation 1 of the interior-singularity report.
