# Test report: `tests/free_particle.yaml` — iteration 2

**Verified against:** the particle in a box (infinite square well), exact closed form.
**Branch / commit:** `TISE-Generalization` @ `9d0f04b` ("test: e2e known solution config files").
**Date:** 2026-09-04.
**Verdict:** **PASS (unchanged).** Bound energies, eigenfunctions, normalisation, continuum states, and phase shifts agree with the exact solution to machine precision for every state the basis can resolve. The one new artefact, `phase_shifts.png`, plots the unwrapped $-kR$-dominated $\delta$ and is uninformative as drawn.

Companion files: `verify_known_solutions.py` (the script that produced every reference number and figure cited below); figures in `figures/`.

---

## 1. Input

```yaml
run:        output_dir: "./data/free_particle"
bspline:    n_nodes: 51,  order: 12,  domain: [0.0, 100.0]
potential:  "{'domain': '[0, 100]', 'function': '0'}"
tise:       n_pts_eigenstate: 301, error_threshold: 1.0e-10
            continuum: enabled: true, E_threshold: 0.0, E_max: 0.5, n_energies: 5, n_pts: 500
visualization: eigenstates: true
```

Node spacing $h = 100/50 = 2$ bohr. Atomic units throughout ($\hbar = m = 1$).

## 2. What the solver actually poses

Three conventions of `tise_solver` decide what "the known solution" is for this input:

1. **Hard walls.** The first and last B-spline are dropped (drop-set $\{1, N\}$), so every eigenstate satisfies $\phi(0) = \phi(100) = 0$. A "free particle" with hard walls is the particle in a box of length $L = 100$.
2. **Discrete positive spectrum.** Because the box is finite, *all* eigenvalues are discrete box states. There are $n_{\rm En} = N_{\rm Bsplines} - 2 = 59$ of them; the `eigenstate_N.png` files run through all 59.
3. **Plot conventions.** `eigenstate_N.png` plots the density $|\phi_n(x)|^2$; `continuum_N.png` plots $\psi_E(x)$ itself. Continuum states are energy-normalised, $\psi_E \to \sqrt{2/\pi k}\,\sin(kx + \delta)$ with $k = \sqrt{2E}$, and the phase shift written to `phase_shifts.dat` is

$$\delta = \arctan\!\left(\frac{k\,\psi(R)}{\psi'(R)}\right) - kR, \qquad R = 100,$$

with the principal-branch arctangent, so $\delta$ is only meaningful modulo $\pi$ and is reported as a large negative number.

## 3. The known solution, step by step

**Step 1 — reduce the equation.** With $V = 0$ the TISE $-\tfrac12\phi'' = E\phi$ becomes $\phi'' = -k^2\phi$ with $k^2 = 2E$. General solution: $\phi = A\sin kx + B\cos kx$.

**Step 2 — apply the walls.** $\phi(0) = 0$ forces $B = 0$. $\phi(L) = 0$ forces $\sin kL = 0$, so $k_n = n\pi/L$, $n = 1, 2, 3, \dots$ Intuitively, the wave is pinned at both ends like a guitar string: only an integer number of half-wavelengths fits.

**Step 3 — energies and normalised states.**

$$E_n = \frac{k_n^2}{2} = \frac{n^2\pi^2}{2L^2} = n^2 \times 4.9348022\times10^{-4}, \qquad \phi_n(x) = \sqrt{\frac{2}{L}}\,\sin\frac{n\pi x}{L}.$$

The density $|\phi_n|^2$ has $n$ bumps of equal height $2/L = 0.02$ and $n-1$ interior zeros.

**Step 4 — continuum states.** At any $E > 0$ the regular solution is $\psi_E = \sqrt{2/\pi k}\,\sin kx$. Nothing scatters it, so the phase shift is $\delta(E) \equiv 0 \pmod\pi$ and $d\delta/dE = 0$ exactly. The amplitude at the five requested energies is $\sqrt{2/\pi k}$ = 1.19312, 1.00329, 0.90657, 0.84366, 0.79788.

## 4. What the solver produced

### 4.1 Bound (box) energies

| n | $E_n$ (solver) | $n^2\pi^2/2L^2$ | relative error |
|---|---|---|---|
| 1 | 4.934802200547e-04 | 4.934802200545e-04 | 4.2e-13 |
| 2 | 1.973920880218e-03 | 1.973920880218e-03 | 1.6e-13 |
| 3 | 4.441321980490e-03 | 4.441321980490e-03 | 5.8e-14 |
| 4 | 7.895683520871e-03 | 7.895683520871e-03 | 3.8e-14 |
| 5 | 1.233700550136e-02 | 1.233700550136e-02 | 2.1e-14 |
| 10 | 4.934802200545e-02 | 4.934802200545e-02 | 5.3e-15 |
| 15 | 1.110330495123e-01 | 1.110330495123e-01 | 4.5e-15 |
| 20 | 1.973920880218e-01 | 1.973920880218e-01 | 5.4e-14 |
| 30 | 4.441322007224e-01 | 4.441321980490e-01 | 6.0e-09 |
| 40 | 7.896170267119e-01 | 7.895683520871e-01 | 6.2e-05 |
| 50 | 1.3327 | 1.2337 | 8.0e-02 |
| 59 | 32.65 | 1.7178 | 18 |

Interpretation: states 1–25 are exact to double-precision round-off. From $n \approx 30$ the error rises steeply and the top few of the 59 states are meaningless ($E_{59} = 32.7$ instead of 1.72). This is the normal signature of a finite basis: 59 order-12 B-splines on 50 intervals can only represent about 30 half-wavelengths accurately. The high box states are not physical content of this test and should not be read as a defect; they should also not be plotted or used downstream (see recommendations).

### 4.2 Eigenfunctions

Solver eigenfunctions were compared to $\sqrt{2/L}\sin(n\pi x/L)$ on the 301-point output grid (overall sign fixed, since eigenvectors are sign-arbitrary).

| n | $\int\phi_n^2\,dx$ | relative $L^2$ difference |
|---|---|---|
| 1 | 1.0000000000000 | 9.1e-14 |
| 2 | 1.0000000000000 | 9.5e-14 |
| 3 | 1.0000000000000 | 8.1e-14 |
| 5 | 1.0000000000000 | 1.2e-13 |
| 10 | 1.0000000000000 | 3.5e-12 |

The `eigenstate_N.png` images (`figures/output_free_particle.png`) show densities with $n$ equal bumps of height 0.0200 and the right node count; the overlay `figures/free_particle_states.png` shows solver and analytic curves indistinguishable.

### 4.3 Continuum states and phase shifts

| E | $\delta$ (raw) | $\delta \bmod \pi$ | $d\delta/dE$ | max $\vert\psi\vert$ | $\sqrt{2/\pi k}$ | rel. $L^2$ error vs $\sin kx$ |
|---|---|---|---|---|---|---|
| 0.1 | −43.9823 = −14π | +4.8e-10 | −6.7e-07 | 1.19312 | 1.19312 | 6.8e-10 |
| 0.2 | −62.8319 = −20π | −6.6e-08 | +8.6e-06 | 1.00328 | 1.00329 | 1.6e-07 |
| 0.3 | −78.5398 = −25π | +1.7e-06 | −1.2e-04 | 0.90657 | 0.90657 | 1.7e-06 |
| 0.4 | −87.9646 = −28π | −2.4e-05 | −6.6e-05 | 0.84366 | 0.84366 | 1.5e-05 |
| 0.5 | −100.5310 = −32π | −1.1e-05 | +1.3e-04 | 0.79797 | 0.79788 | 1.0e-04 |

All five raw phase shifts are integer multiples of $\pi$, i.e. exactly zero modulo $\pi$ as required. The amplitudes match the energy normalisation to five digits. The wavefunction error grows smoothly from $10^{-9}$ to $10^{-4}$ as $E$ increases: at $E = 0.5$, $k = 1$, the wavelength is $2\pi \approx 6.3$ bohr and $h = 2$ gives only 3.1 nodes per wavelength. See `figures/free_particle_continuum.png`.

The solver also emitted "pole proximity" warnings for $E = 0.1, 0.2, 0.3, 0.5$ (each lies within $\lesssim 0.008$ of a box eigenvalue). These are expected: with $L = 100$ the box levels near $E \sim 0.3$ are only $\Delta E \approx k\pi/L \approx 0.024$ apart, so *any* energy on a coarse grid is close to one. The warnings did not affect the results here — the continuum state built at a near-pole energy is still an exact solution at that energy.

## 5. Correctness assessment

| Quantity | Status |
|---|---|
| Box energies $n \le 25$ | correct to ~$10^{-13}$ |
| Box energies $n \gtrsim 30$ | basis-limited, increasingly wrong; not physical content |
| Eigenfunction shape and normalisation | correct to ~$10^{-13}$ |
| Continuum wavefunctions, amplitude | correct to $10^{-9}$–$10^{-4}$ (resolution-limited at high $E$) |
| Phase shifts $\delta \bmod \pi$ | correct (all ≡ 0) |
| $d\delta/dE$ | correct (≈ 0) — but only because the exact answer is a constant; see hydrogen and finite-well reports |

**Overall: PASS.** This test exercises the kinetic-energy matrix, the overlap matrix, the generalised eigensolver, the eigenstate reconstruction, the continuum construction and the asymptotic matching, and every one of them behaves exactly as the closed-form solution demands.

## 6. What changed at 8236239 for this test

* **Output data:** `eigenvalues.dat`, all `eigenstate_NNN.dat`, `continuum_state_NNN.dat` and `phase_shifts.dat` are byte-identical to 9d0f04b. `buildStrategicGridAndDropSet` finds no interior join in `[0,100] → 0` and returns the same uniform grid and `{1}` drop-set.
* **New diagnostics:** `warnings.json` now carries "0 of 59 computed states are below E=0.0" (`classifyBoundStates`), which is correct. No well-containment warnings, correct (no bound states to check).
* **New plot:** `visualization.phase_shifts: true` was added to this YAML and `analysis.py` now writes `phase_shifts.png` (`figures/output_phase_shifts_fp_fsw.png`, left). The upper panel shows the raw $\delta$ from −44 to −101, i.e. the $-kR$ line; the fact that it is exactly $-14\pi, -20\pi, \dots$ (the physics) is invisible. The lower panel shows $d\delta/dE \sim 10^{-4}$, which for the free particle is correctly "zero" but only because the exact answer is constant.
* **Round-1 recommendations status:** (1) wrap $\delta$ — not done, now visible in the plot; (2) cap plotted states — not done (59 PNGs); (3) tighten $E_{\rm acc}$ — not done (`minInterNodeGap` is now used, which is the right input, but the formula and the missing $V_{\min}$ are unchanged); (4) denser continuum grid — not done; (5) automated check — partially: the YAML is now loaded by an e2e test and PNG existence is asserted, but no energy is compared to $n^2\pi^2/2L^2$.

## 7. Recommendations

1. **Wrap the reported phase shift.** `phase_shifts.dat` stores $\arctan(\cdot) - kR$, which is $-14\pi$ etc. here. Storing $\delta$ reduced to $(-\pi/2, \pi/2]$ (or better, tracking the branch across the energy grid so $\delta(E)$ is continuous) would make the file readable and make $d\delta/dE$ well-defined.
2. **Cap what `analysis.py` plots.** All 59 eigenstates are plotted although only ~25 are accurate. Either plot only states below a user-set energy, or below the basis accuracy limit $E_{\rm acc}$ (see recommendation 3), or apply `error_threshold` the way `main.cpp`'s hydrogen demo does.
3. **Tighten $E_{\rm acc}$.** `computeEAcc` returns $\pi^2/2h^2 = 1.23$ for $h = 2$, which is the Nyquist limit of 2 nodes per wavelength. The data above (and in the oscillator report) show that ~3.3 nodes per wavelength are needed for $10^{-6}$ accuracy, i.e. a practical ceiling about $E_{\rm acc}/2.7$. Consider reporting both, or warning at a fraction of the Nyquist value.
4. **Denser continuum grid.** Five energies over $[0, 0.5]$ are fine for a smoke test but too few for anything derived from $\delta(E)$; this matters for the other three tests more than here.
5. **Add an automated check.** The comparison in `verify_known_solutions.py` (energies vs $n^2\pi^2/2L^2$, $\delta \bmod \pi$ vs 0, amplitude vs $\sqrt{2/\pi k}$) is cheap and deterministic; it would make a good end-to-end pytest for this YAML.
