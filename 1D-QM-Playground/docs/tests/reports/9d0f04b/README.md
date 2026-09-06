# Known-solution test reports — iteration 1 (commit `9d0f04b`)

Verification of the `tests/*.yaml` end-to-end configurations against exact solutions, at commit `9d0f04b` (`TISE-Generalization`), 2026-09-04.

| YAML | Report | Verdict |
|---|---|---|
| `tests/free_particle.yaml` | [free_particle.md](free_particle.md) | PASS — exact to 1e-13 |
| `tests/finite_square_well.yaml` | [finite_square_well.md](finite_square_well.md) | PASS on physics, marginal numerics on shipped grid (h³ convergence at the potential step; near-threshold box states misplaced; dδ/dE unusable) |
| `tests/harmonic_oscillator.yaml` | [harmonic_oscillator.md](harmonic_oscillator.md) | PASS — exact to 1e-6 for n ≤ 7, resolution-limited above |
| `tests/hydrogen.yaml` | [hydrogen.md](hydrogen.md) | PASS — ℓ = 1 spectrum and Coulomb waves; δ is flat-matched at R = 100; dδ/dE unusable |

Also here:

* `known_solutions_report.html` — the four tests in one self-contained document (figures embedded, MathJax loaded from cdnjs; open in a browser).
* `verify_known_solutions.py` — computes every reference value and overlay figure. Run from `1D-QM-Playground/` after `python3 controller.py --config tests/<name>.yaml` for all four tests; set `D` at the top to the `data/` directory and `OUT` to where figures should go. Needs numpy, scipy, matplotlib, mpmath.
* `figures/` — overlay plots (solver vs analytic) and contact sheets of the solver's own `eigenstate_N.png` / `continuum_N.png` output.

Cross-cutting findings (details in each report):

1. `computeEAcc = π²/2h²` is the Nyquist limit (2 nodes per wavelength) and ignores `V_min`; the data indicate ~3.3 nodes per local wavelength for 1e-6 accuracy, i.e. `h ≲ 2π / (3.3 √(2 (E − V_min)))`.
2. The `dδ/dE` column in `phase_shifts.dat` (central difference of sin 2δ over ΔE = 0.1, divided by cos 2δ) is noise for every potential except the free particle.
3. Repeated-knot / strategic node placement exists in `tise.cpp` but `tise_solver` always builds a uniform grid, so potentials with steps converge only algebraically.
4. `analysis.py` plots every eigenstate, most of which are basis-limited box states; gate on an energy ceiling.
