# Test report: `tests/interior_singularity.yaml` — iteration 2 (new test)

**Verified against:** the exact solution of a domain split by a non-integrable singularity: a particle in a box on $[0,20]$ plus a repulsive Coulomb region on $[20,40]$, whose spectra are the union of $n^2\pi^2/800$ and the zeros of the regular Coulomb function $F_0(1/k, 20k)$.
**Branch / commit:** `TISE-Generalization` @ `8236239`. First appearance of this test.
**Date:** 2026-09-04.
**Verdict:** **FAIL.** The A4b "singular-join B-spline removal" that this test was written to exercise produces eigenvalues that are wrong by **41 %** for the states on the field-free side ($E_1 = 0.01742$ vs exact $0.01234$) and by 0.1–1 % on the Coulomb side. The mechanism is identified below: dropping every B-spline whose support touches $x=20$ on a simple-knot grid forces $\phi \equiv 0$ on the whole interval $[19, 21]$ and forces $\phi$ to vanish like $(19-x)^7$ beyond it. The correct treatment — a knot of multiplicity $\text{order}-1$ at the singular point and removal of the *single* B-spline that is non-zero there — reproduces all twelve exact eigenvalues to $10^{-11}$ with the same node count; this is demonstrated numerically in §5. The test's own assertions (finite eigenvalues; basis smaller than 45) pass, which is why the defect was not caught.

Companion files: `verify_new_tests.py`, `figures/interior_singularity.png`.

---

## 1. Input

```yaml
run:        run_analysis: false,  output_dir: "./data/interior_singularity"
bspline:    n_nodes: 41,  order: 8,  domain: [0.0, 40.0]
potential:  "{'domain': '[0, 20)',  'function': '0'}"
            "{'domain': '(20, 40]', 'function': '1/(x-20)'}"
tise:       n_pts_eigenstate: 301
            continuum: enabled: true, E_threshold: 0.0, E_max: 0.5, n_energies: 5, n_pts: 500
```

Node spacing $h=1$; $x=20$ is a grid node. No plots are produced (`run_analysis: false`), so the figures here come from the `.dat` files.

## 2. What the solver actually poses

`detectPotentialStructure` classifies the join at $x=20$ as **Singular**. Per `strategicKnotsFromJoins` a Singular join gets *no* extra knot; its remediation is `bSplinesTouchingX(20)`: every B-spline whose support $[t_{i-8}, t_i]$ contains 20 — nine of them for order 8 on simple knots — is dropped from the basis. With $B_1$ and $B_N$ also removed, 47 − 9 − 2 = 36 states are computed (the number the e2e test asserts is below 45).

The consequence, which the test does not check, is what the surviving basis can represent near $x=20$. On a simple-knot grid with spacing 1 the last kept B-spline to the left has support $[11,19]$ and the first to the right has support $[21,29]$. Therefore:

* every basis function — and hence every eigenfunction — is **identically zero on $[19, 21]$**;
* approaching 19 from the left the basis vanishes like $(19-x)^{7}$ (a degree-7 spline meeting zero with $C^6$ continuity), so the eigenfunctions cannot have finite slope at the effective wall.

The solver is therefore solving a box of length 19 with a "soft" wall that suppresses the wavefunction over the last few bohr, not the problem the YAML describes.

## 3. The known solution, step by step

**Step 1 — what a non-integrable singularity does.** $V = 1/(x-20)$ for $x>20$ diverges like $1/s$, $s = x-20$. $\int_0 ds/s$ diverges, so $\langle\phi|V|\phi\rangle$ is finite only if $\phi(20)=0$ — the wavefunction must vanish at the singular point (linearly, as the $s^1$ behaviour of the regular solution). A one-sided $1/s$ singularity is thus an **impenetrable wall**: the two sides decouple completely, and the spectrum is the union of two independent problems.

**Step 2 — the left side, $[0,20]$, $V=0$, $\phi(0)=\phi(20)=0$.** A particle in a box of length 20:

$$E_n^{\rm box} = \frac{n^2\pi^2}{2\cdot20^2} = 0.0123370\,n^2 : \quad 0.012337,\ 0.049348,\ 0.111033,\ 0.197392,\ 0.308425,\ 0.444132,\ \dots$$

**Step 3 — the right side, $[20,40]$, $V=1/s$, $u(0)=u(20)=0$.** This is the $\ell=0$ radial equation with a *repulsive* Coulomb potential ($Z=-1$), whose regular solution is the Coulomb function $F_0(\eta, ks)$ with Sommerfeld parameter $\eta = +1/k$. $F_0$ is an entire function of $ks$ starting as $C_0(\eta)\,ks\,(1+\eta ks+\dots)$ — it vanishes linearly at $s=0$, no logarithm (that lives in the irregular $G_0$). The hard wall at $s=20$ quantises: $F_0(1/k, 20k)=0$. Roots (mpmath, 25 digits, bracketed on a scan and polished with Brent):

$$E^{\rm Coul} = 0.100541,\ 0.168654,\ 0.254096,\ 0.359297,\ 0.485622,\ \dots$$

These lie *above* the corresponding empty-box values on $[20,40]$ because the Coulomb barrier costs energy near $s=0$.

**Step 4 — the exact spectrum** is the sorted union, with each state living entirely on one side:

| # | E | side |
|---|---|---|
| 1 | 0.0123370 | box |
| 2 | 0.0493480 | box |
| 3 | 0.1005405 | Coulomb |
| 4 | 0.1110330 | box |
| 5 | 0.1686541 | Coulomb |
| 6 | 0.1973921 | box |
| 7 | 0.2540958 | Coulomb |
| 8 | 0.3084251 | box |
| 9 | 0.3592967 | Coulomb |
| 10 | 0.4441322 | box |
| 11 | 0.4856224 | Coulomb |
| 12 | 0.6045133 | box |

There are no bound states ($V\ge0$ everywhere), consistent with the solver's "0 of 36 below E=0".

## 4. What the solver produced

| # | $E$ (solver) | exact | side (solver eigenvector) | rel. error |
|---|---|---|---|---|
| 1 | 0.017421 | 0.012337 | left | **+41 %** |
| 2 | 0.069726 | 0.049348 | left | **+41 %** |
| 3 | 0.100627 | 0.100541 | right | +0.09 % |
| 4 | 0.157052 | 0.111033 | left | **+41 %** |
| 5 | 0.170372 | 0.168654 | right | +1.0 % |
| 6 | 0.263653 | 0.197392 | right* | — |
| 7 | 0.279629 | 0.254096 | left* | — |
| 8 | 0.387645 | 0.308425 | right* | — |
| 9 | 0.437800 | 0.359297 | left* | — |

(*From state 6 on the solver's ordering no longer matches the exact ordering because the left-side states are shifted up by 41 % and interleave differently; the side column is where the solver's eigenvector actually lives.)

The left-side states are all high by the same factor $1.412 \approx (20/19)^2 \times 1.27$: a box of length 19 would give $+10.8\,\%$; the additional factor comes from the $(19-x)^7$ suppression, which acts like a further soft wall. The right-side states are much less affected because $F_0$ already vanishes at $s=0$ and its first lobe peaks 8 bohr away, so the dead zone $[20,21]$ costs little.

Every eigenvector is exactly zero on $[19.07, 20.93]$ in the `eigenstate_NNN.dat` output (the output grid has no point at exactly 19 or 21). `figures/interior_singularity.png`, left panel: the solver's ground state (blue) is squeezed into $[0,19]$ and bows away from the exact $\sqrt{2/20}\sin(\pi x/20)$ (grey); the shaded band is the dead zone. Middle panel: the Coulomb-side ground state is nearly right because it is naturally small there. Right panel: the whole spectrum, solver vs exact.

The continuum output at $E=0.1$–$0.5$ is not analysed further: with a dead zone at $x=20$ the "continuum state" is a superposition of two decoupled boxes and the flat-asymptote matching at $x=40$ — where $V = 0.05$, not flat — is not meaningful either.

## 5. What the correct treatment gives (demonstration)

The right way to impose $\phi(x_s)=0$ at an interior point in a B-spline basis is the same one the code already uses at the domain edges: make $x_s$ a knot of multiplicity $k-1$ (here 7), so that the basis is merely $C^0$ there and exactly **one** B-spline is non-zero at $x_s$, then drop that one. Every other basis function is unaffected; the two halves become independent bases that each approach $x_s$ with full polynomial freedom (linear, quadratic, …), which is what the true solutions need.

`verify_new_tests.py` implements this with `scipy.interpolate.BSpline` on the same uniform $h=1$ grid, order 8, Gauss–Legendre quadrature per interval, generalized symmetric eigensolve (50 basis functions kept):

| # | exact | mult-7 knot at 20 + drop 1 | rel. error |
|---|---|---|---|
| 1 | 0.0123370055 | 0.0123370055 | 1.7e-15 |
| 2 | 0.0493480220 | 0.0493480220 | 1.1e-14 |
| 3 | 0.1005405216 | 0.1005405216 | 2e-14 |
| 5 | 0.1686540851 | 0.1686540851 | 1.9e-13 |
| 8 | 0.3084251375 | 0.3084251375 | 1.4e-12 |
| 12 | 0.6045132696 | 0.6045132698 | 3.2e-10 |

For comparison the *naive* treatment (simple knots, drop nothing, integrate $1/(x-20)$ with Gauss–Legendre, i.e. approximately the pre-8236239 behaviour) gives 0.012258 for state 1 — 0.6 % low. So the pre-fix result was mildly wrong and the "fix" is severely wrong.

## 6. Correctness assessment

| Quantity | Status |
|---|---|
| bound-state count (0) | correct |
| basis reduction (36 < 45) | as designed — but the design is the problem |
| left-side (box) eigenvalues | **wrong by 41 %** |
| right-side (Coulomb) eigenvalues | 0.1–1 % high |
| eigenfunctions | identically zero on $[19,21]$; wrong shape near 19 and 21 |
| continuum / phase shifts | not meaningful for this configuration |

**Overall: FAIL — the A4b remediation as wired at 8236239 is not a correct discretisation of an interior singular join.**

## 7. Recommendations

1. **Replace `bSplinesTouchingX`-removal by knot-multiplicity + single-drop.** In `buildStrategicGridAndDropSet`, for each Singular join at $x_s$: add $x_s$ to the strategic knots with `extraMultiplicity = order − 2` (total multiplicity $k-1$ if $x_s$ is already a grid point), rebuild the grid, then drop only the B-spline index $i$ with $B_i(x_s) \ne 0$ (there is exactly one). This is what the domain edges already get for free from the clamped end knots, and it is why hydrogen's $r=0$ singularity has never caused trouble. Expected effect: all twelve eigenvalues in §4 exact to $\sim10^{-11}$ (as in §5) with `nEnBound = 47 + 6 − 3 = 50`.
2. **If the singularity is one-sided (as here) consider multiplicity $k$ instead of $k-1$**, which makes the two halves completely independent bases (no shared function at all) and lets the two sides even use different node densities. Either works for the eigenvalues; multiplicity $k$ is cleaner for a truly impenetrable point.
3. **Make the e2e test discriminate.** `test_bspline_removal_actually_shrinks_the_basis_below_classic` guards the plumbing, not the physics. Add: `abs(E_1 − π²/800) < 1e-9` and `abs(E_3 − 0.1005405216) < 1e-8` (the first Coulomb-side root). With the current code these fail by 41 % and 0.09 %.
4. **Reconsider the continuum for split domains.** When a Singular interior join is detected, the continuum construction should either be refused (the right box is not connected to the left one) or be built per sub-domain. At minimum, warn.
5. **Add a right-edge test of the same physics** with a *smooth* interior repulsion (e.g. `1/((x-20)^2+1)`) so that the Step/Kink/Singular detection thresholds are also exercised on a case that must *not* trigger removal.
