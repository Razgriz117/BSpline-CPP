# Known-solution test reports — iteration 3 (commit `f4e8359`)

Verification of every `tests/*.yaml` at commit `f4e8359` (`TISE-Generalization`, "implement Coulomb-tail continuum matching…"), 2026-09-11. Covers the twelve commits since `8236239`, chiefly `0bf17be` ("fix: correct interior singular-join B-splines, dDeltaDE, and Case-3/right-edge continuum handling"), `5cd28e4` (raw-wavefunction plots) and `f4e8359` (Coulomb-tail matching). Earlier iterations: `../8236239/`, `../9d0f04b/`.

| YAML | Report | Verdict at f4e8359 | Change vs 8236239 |
|---|---|---|---|
| `free_particle.yaml` | [free_particle.md](free_particle.md) | PASS | data identical; δ now stored ≈ 0 (unwrapped); new dδ/dE ≤ 1e-2 vs exact 0 (basis error × 1/2·fineDE) |
| `finite_square_well.yaml` | [finite_square_well.md](finite_square_well.md) | **PASS, nothing open** | **dδ/dE fixed** (−19.432 vs exact −19.432); δ unwrapped matches continuous trajectory to 3e-5; `phase_shifts.png` now informative |
| `harmonic_oscillator.yaml` | [harmonic_oscillator.md](harmonic_oscillator.md) | PASS | data identical; plots opt into \|φ\|² via new `bound_states_squared` |
| `hydrogen.yaml` | [hydrogen.md](hydrogen.md) | PASS | **Coulomb-tail matching**: δ = 0.001–0.002 vs exact 0 (was 0.02–0.07 off); ψ_E to 2e-5; residual traced 100 % to the hand-rolled F/G evaluator (Wronskian −0.9982, O(h²) derivative, WKB start) |
| `interior_singularity.yaml` | [interior_singularity.md](interior_singularity.md) | **PASS** (was FAIL 41 %) | knot multiplicity k−1 + single drop, exactly as recommended; all 12 eigenvalues to ≤ 1.4e-12; continuum still written for a split domain |
| `right_edge_singularity.yaml` | [right_edge_singularity.md](right_edge_singularity.md) | **PASS, nothing open** | continuum now refused (no files written) instead of warned |
| `case3_irregular_tail.yaml` | [case3_irregular_tail.md](case3_irregular_tail.md) | **PASS** | taper gated on `continuum.enabled`; eigenvalues match untapered reference to ≤ 4e-8 (was shifted 1e-5…1e-4) |

Also here:

* `verify_known_solutions.py`, `verify_new_tests.py` — reference values and overlays (as in iteration 2, output paths updated; the right-edge script no longer expects continuum files).
* `coulomb_wave_check.cpp`, `coulomb_wave_stepsize.cpp` — small programs linking `libtise_lib.a` that compare `evaluateCoulombFunctions` with `mpmath` and scan its step size / start radius (build line in the hydrogen report §4.3; needs the `TISE/build` static libraries).
* `figures/` — overlays for all seven tests, plus `output_phase_shifts_fsw.png` (the new, meaningful phase-shift plot) and `output_plots_f4e8359.png` (raw-wavefunction plots).

Cross-cutting state at f4e8359:

1. **Every iteration-1 and iteration-2 defect that was a solver defect is fixed:** step joins (8236239), interior singular joins, dδ/dE, δ unwrapping, Case-3 gating, right-edge refusal (0bf17be), Coulomb-tail matching (f4e8359). All seven tests pass on physics.
2. **New accuracy floor:** the Coulomb-function evaluator limits hydrogen's δ to ~2e-3 rad. Fixes in order of effort: report the Wronskian; use an O(h⁴) derivative consistent with Numerov; replace the WKB start by the DLMF 33.11 asymptotic series or Steed's continued fractions (or link GSL). Expected gain: δ → ≲1e-6.
3. **Still open (minor):** dδ/dE amplifies basis error by 1/(2·fineDE) — Richardson or a larger step; continuum output for split domains should be refused; the Case-3 classifier's "will be tapered" message contradicts the orchestrator's "not tapered" when the continuum is off; `computeEAcc` ignores V_min; `analysis.py` plots every basis-limited state; free-particle and oscillator e2e tests still have no analytic assertions; plotted eigenstates have arbitrary sign.
