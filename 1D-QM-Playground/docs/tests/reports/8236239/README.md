# Known-solution test reports — iteration 2 (commit `8236239`)

Verification of every `tests/*.yaml` configuration at commit `8236239` (`TISE-Generalization`, "TISE release-readiness pass — unify grid/dropset, wire diagnostics, add CI"), 2026-09-04. Iteration 1 (commit `9d0f04b`) is in `../9d0f04b/`.

| YAML | Report | Verdict at 8236239 | Change vs 9d0f04b |
|---|---|---|---|
| `free_particle.yaml` | [free_particle.md](free_particle.md) | PASS — exact to 1e-13 | output byte-identical; new `phase_shifts.png` uninformative |
| `finite_square_well.yaml` | [finite_square_well.md](finite_square_well.md) | **PASS** — bound energies 1e-13…1e-9, δ to 2e-4 rad | **fixed**: strategic knot at the step now reaches `tise_solver` (was 1e-4…1e-2 / 0.6 rad; box states misplaced) |
| `harmonic_oscillator.yaml` | [harmonic_oscillator.md](harmonic_oscillator.md) | PASS — exact to 1e-6 for n ≤ 7 | output byte-identical |
| `hydrogen.yaml` | [hydrogen.md](hydrogen.md) | PASS — ℓ=1 spectrum, Coulomb waves | output byte-identical; new wall-collision diagnostic flags exactly 6p/7p/8p (correct) |
| `interior_singularity.yaml` (new) | [interior_singularity.md](interior_singularity.md) | **FAIL** — field-free-side eigenvalues 41 % high | A4b singular-join B-spline removal zeroes ψ on [19, 21]; correct treatment demonstrated (knot multiplicity k−1 + drop one) |
| `right_edge_singularity.yaml` (new) | [right_edge_singularity.md](right_edge_singularity.md) | PASS — 17 eigenvalues to 1e-12 vs zeros of F₀ | continuum output unphysical (ψ(100) ≠ 0), correctly warned; should be refused |
| `case3_irregular_tail.yaml` (new) | [case3_irregular_tail.md](case3_irregular_tail.md) | PASS on numerics — matches FD solve of the tapered V to 1e-8 | taper shifts the box spectrum by 1e-5…1e-4 even with continuum disabled |

Also here:

* `verify_known_solutions.py` — reference values and overlays for the four original tests (unchanged from iteration 1 apart from output paths).
* `verify_new_tests.py` — references for the three new tests: exact split-domain spectrum (box ∪ Coulomb zeros), Coulomb-function zeros for the right-edge case, Richardson-extrapolated finite-difference eigenvalues for the tapered and raw Case-3 potentials, and a scipy B-spline demonstration of the correct interior-singularity treatment. Run both from `1D-QM-Playground/` after running all seven YAMLs through `controller.py`; set `D` at the top to the `data/` directory. Needs numpy, scipy, matplotlib, mpmath.
* `figures/` — overlay plots for all seven tests plus contact sheets of the solver's own PNG output where `run_analysis` was on.

Cross-cutting findings at 8236239:

1. **Fixed since iteration 1:** potential steps (`finite_square_well`) now get a C¹ knot; hydrogen's box-confined states are now flagged by `checkWellContainment`; `computeEAcc` now takes `minInterNodeGap`.
2. **New defect:** interior Singular joins are remediated by dropping every B-spline touching the point on a simple-knot grid, which kills the basis on a 2h-wide band and forces (k−1)-order vanishing beyond it. Replace with knot multiplicity k−1 (or k) at the point plus removal of the single non-zero B-spline — the same treatment the domain edges already receive, which is why `right_edge_singularity` passes to 1e-12 and `interior_singularity` fails by 41 %.
3. **Still open from iteration 1:** `dδ/dE` (finite difference of sin 2δ over ΔE = 0.1, divided by cos 2δ) is noise for every non-trivial potential — the Part-C investigation established the formula's sign, not its applicability at the shipped grid spacing; `phase_shifts.dat`/`.png` store and draw the raw −kR-dominated δ; `computeEAcc` ignores V_min and warns only at the Nyquist limit; `analysis.py` plots every basis-limited state; no e2e test compares an energy to an analytic value.
4. **Modelling choice to revisit:** the Case-3 taper is applied whether or not the continuum is requested, and its width 0.1 L is not tied to any physical scale.
