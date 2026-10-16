# Known-solution test reports

One directory per iteration, named by the short hash of the solver commit the runs were made against. Each contains one report per `tests/*.yaml` (same base name), a README with the verdict table, the verification scripts that produced every reference number, and `figures/`.

| Iteration | Commit | Tests covered | Headline |
|---|---|---|---|
| 1 | [`9d0f04b`](9d0f04b/README.md) | free_particle, finite_square_well, harmonic_oscillator, hydrogen | 3 PASS; finite well physically right but numerically under-resolved (step in V on simple knots) |
| 2 | [`8236239`](8236239/README.md) | the four above + interior_singularity, right_edge_singularity, case3_irregular_tail | finite well fixed (1e-13); **interior_singularity FAIL** (41 %, A4b removal defect); right-edge and Case-3 pass on numerics with modelling caveats |
| 3 | [`f4e8359`](f4e8359/README.md) | all seven | **all PASS**: interior singularity, dδ/dE, δ unwrapping, Case-3 gating, right-edge refusal fixed; Coulomb-tail matching brings hydrogen δ to 2e-3 rad, floor traced to the F/G evaluator |


Review verifications (same method applied to a code review rather than to a test YAML):

| Directory | Subject | Headline |
|---|---|---|
| [`pr-18-433e11a`](pr-18-433e11a/README.md) | `docs/planning/PR-18-phy5606-projects-review.md` (PR #18 @ `433e11a`) | all 6 correctness findings confirmed e2e; near-duplicate deltas corrupt silently by up to 10¹² with Cholesky still succeeding — the review's `LLT` fix would not catch it, delta validation + a cond(S) guard would |

Iteration 1 also has a self-contained HTML version (`9d0f04b/known_solutions_report.html`).
