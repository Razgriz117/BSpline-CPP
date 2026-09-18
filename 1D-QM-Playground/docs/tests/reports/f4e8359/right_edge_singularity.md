# Test report: `tests/right_edge_singularity.yaml` — iteration 3

**Verified against:** repulsive Coulomb in a box, $V = 1/(100-x)$: eigenvalues are the zeros of $F_0(1/k, 100k)$, eigenfunctions $F_0\big(1/k, k(100-x)\big)$.
**Branch / commit:** `TISE-Generalization` @ `f4e8359`. **Previous:** `8236239`.
**Date:** 2026-09-11.
**Verdict:** **PASS, and the open item is closed.** Eigenvalues and eigenfunctions are byte-identical to 8236239 (17 states to $\le 2\times10^{-10}$, functions to $10^{-8}$). The solver now **refuses** continuum construction on a singular right edge instead of warning and writing unphysical states: no `phase_shifts.dat` or `continuum_state_NNN.dat` is produced, and `warnings.json` says why. Iteration-2 recommendation 1, implemented.

Companion files: `verify_new_tests.py`, `figures/right_edge_singularity.png`.

---

## 1. Input

Unchanged: 41 nodes, order 8, $[0,100]$, `continuum: enabled: true` with three energies, `run_analysis: false`.

## 2. What changed in the solver (commit `0bf17be`)

When `sgr.rightEdgeSingular` is true the continuum block is skipped entirely. The warning now reads: "potential is singular at the right domain edge x=100; continuum phase-shift matching (matchAsymptotic) assumes a regular boundary there, which does not hold here — skipping continuum construction entirely (no phase_shifts.dat/continuum_state_NNN.dat written)". The e2e test asserts $E_1 = 0.0149922127$.

## 3. The known solution

`../8236239/right_edge_singularity.md` §3. Zeros of $F_0(1/k,100k)$: $0.0149922,\ 0.0203389,\ 0.0260832,\ 0.0324237,\ 0.0394463,\ \dots$

## 4. What the solver produced

| j | $E$ (solver) | exact | rel. error |
|---|---|---|---|
| 1 | 0.01499221266094 | 0.01499221266093 | 3.9e-13 |
| 2 | 0.02033885931642 | 0.02033885931660 | −8.8e-12 |
| 3 | 0.02608321293993 | 0.02608321293992 | 1.8e-13 |
| 12 | 0.11083564751346 | 0.11083564749335 | 1.8e-10 |

Eigenfunction $L^2$ differences vs $F_0$: $1.3\times10^{-11}$ (state 1), $1.4\times10^{-10}$ (2), $8.9\times10^{-9}$ (5). `figures/right_edge_singularity.png` now shows states 1 and 5 over the exact Coulomb functions (the former right panel, the unphysical continuum state, no longer exists to plot).

Output directory contents: `eigenvalues.dat`, `eigenvectors.dat`, 45 `eigenstate_NNN.dat`, `hamiltonian.dat`, `overlap.dat`, `warnings.json` — no continuum files. `analysis.py` (if `run_analysis` were on) tolerates their absence.

## 5. Correctness assessment

| Quantity | Status |
|---|---|
| eigenvalues, eigenfunctions | correct (unchanged) |
| continuum on a singular edge | correctly refused |
| warning text | accurate |

**Overall: PASS. Nothing open.**

## 6. Recommendations

1. Carried from iteration 2, now lower priority: if scattering off a singular wall is ever wanted, mirror the domain so the singularity is the origin and the regular boundary is the outer edge, then use the new Coulomb-tail matching with $C=+1$ (repulsive) — `matchAsymptotic`'s Coulomb branch already handles the sign through $\eta = C/k$.
2. Keep this YAML as the guard that domain-edge and interior singular treatments agree; with `interior_singularity` now exact, the two tests together pin the $1/s$ physics at both placements.
