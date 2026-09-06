# TISE release-readiness: gap closure + e2e test plan

> Companion doc to `engineer-a-plan-*.md`/`engineer-a-plan-A4-wiring-design.md`, same convention: plan first, "Implementation status (post-execution)" appended once the work actually lands. Produced from a 3-agent audit (requirements/ADRs/planning docs; C++ implementation; Python pipeline + test coverage) run against the repo as of 2026-08-28, comparing it to `docs/SDD.md`'s requirements (REQ-F-010..050, REQ-NF-010), the 8 existing ADRs, and the planning docs, ahead of shipping the TISE pipeline end-to-end before TDSE work begins.

## Context

Two scope decisions were made explicitly with the user before finalizing this plan: (1) the interior-potential-discontinuity gap (G1 below) gets a **full fix** — unifying `solveTISE`'s and `tise_solver`'s grid+dropset construction — not just a warning; (2) **CI gets set up** as part of this pass, since none currently exists.

**Methodology note, worth recording**: the requirements/ADR audit agent was documentation-only, and several of its claims (that the `E_acc` warning doesn't reach `warnings.json`, that `tise.continuum.E_threshold`/`E_max` aren't wired to config) are **stale** — they describe a state from `docs/planning/engineer-a-plan-A5.md` that predates this session's own wiring work (`tise_solver_main.cpp` already reads these fields and already routes warnings through `addWarning`→`warnings.json`, confirmed directly in code both by the C++-focused audit agent and by this session's own prior direct verification). Only the code-focused audit's findings and this session's own directly-verified facts are treated as ground truth in the gap inventory below — this is itself evidence that `docs/planning/`'s per-task status sections need a periodic re-verification pass, not just a one-time write.

## Gap inventory

| # | Gap | Source | Disposition |
|---|---|---|---|
| G1 | `tise_solver_main.cpp` always uses a plain uniform grid + classic `{1}` drop-set — never `detectPotentialStructure`/strategic node placement/A4b singular-B-spline removal, unlike `solveTISE`. An interior Step/StitchedKink/Singular potential join produces **silently wrong eigenvalues with `EXIT_SUCCESS`**, zero warning. Confirmed: zero test anywhere (C++ or Python) exercises this path. | ADR-0008; C++ audit §1/§3/§6 | **Fix — Part A** |
| G2 | Neither `README.md` nor `TISE/README.md` documents the real pipeline (`config.yaml`→`controller.py`→`tise_solver`→`analysis.py`) — both describe only the legacy `H-BoundStates` JSON-argv demo; `tise_solver` isn't named in either file. | Python audit §3 | **Fix — Part E** |
| G3 | `visualization.phase_shifts` (default `true`, documented TISE-only/non-TDSE-gated) has **zero plotting implementation** in `analysis.py` — silently produces no plot. | Python audit §1 | **Fix — Part D** |
| G4 | The 4 known-solution configs (`tests/{free_particle,finite_square_well,harmonic_oscillator,hydrogen}.yaml`) are **never loaded by any test** — orphaned. | Python audit §5 item 2 | **Fix — Part G** |
| G5 | No integration test asserts real PNG/plot output exists. | Python audit §5 item 1 | **Fix — Part G** |
| G6 | No CI configuration anywhere in the repo. | Python audit §6 | **Fix — Part F** |
| G7 | `matchAsymptotic`'s `dDeltaDE` finite-difference has a self-flagged `// TODO: will this sign be wrong` at the last energy-grid point of every continuum run (`tise.cpp:640`). | C++ audit §4 | **Fix — Part C** |
| G8 | `evaluateFunction` resolves overlapping potential-piece domains via `std::map` iteration order, no validation. | C++ audit §4 | **Fix — Part C** |
| G9 | `classifyBoundStates`/`checkWellContainment` — fully implemented, fully unit-tested, zero production callers. | C++ audit §1 | **Fix — Part B** |
| G10 | `evaluateWindowedPotential`/`case3WindowFunction` — Case-3 detected+warned but never remediated. | C++ audit §1 | **Fix — Part B** |
| G11 | `minInterNodeGap` has zero production callers, even inside `solveTISE`. | C++ audit §1 | **Fix — Part B** (falls out of Part A) |
| G12 | `physics.mass`/`physics.hbar` documented, never consumed, no error if set. | C++ audit §2 | **Fix — Part C** (guard-rail) |
| G13 | `buildEnergyGrid` — critical path, zero direct unit test. | C++ audit §5 | **Fix — Part G** |
| G14 | No dedicated explicit real-subprocess test for `run.run_analysis: false`. | Python audit §5 item 3 | **Fix — Part G** |
| G15 | SDD §7.1 never lists `controller.py`'s own CLI contract; §5.1.3 pseudocode stale. | Python audit §2 | **Fix — Part E** |
| G16 | Stale comment (`tise.hpp:299`) says A4b is "out of scope" — implemented since. | C++ audit §4 | **Fix — Part E** |
| — | Coulomb-tail continuum matching (REQ-F-030), FEDVR (ADR-0001), WKB density (ADR-0002), multi-particle/3D (ADR-0003), CAP/ECS, TDSE-dependent `analysis`/`visualization` fields | requirements audit §1/§4 | **Explicitly deferred, no action** |

## Part A: Unify `solveTISE`/`tise_solver` grid + drop-set construction (closes G1, G11)

Extract `solveTISE`'s grid/dropset-construction block (`TISE/tise.cpp:1081-1146`, from `detectPotentialStructure` through computing `nEnBound`, before `fillBandedMatrices` is called) into a new shared function:

```cpp
// tise.hpp
struct StrategicGridResult {
    std::vector<Real> grid;
    bspline::BSpline bs;
    int nBSplines;
    std::vector<int> fillDropSet;   // for fillBandedMatrices; B_N excluded (kept for continuum)
    int nEnBound;
    bool rightEdgeSingular;
};
StrategicGridResult buildStrategicGridAndDropSet(int nNodes, int order, Real rMin, Real rMax,
                                                   const std::map<std::string, std::string> &potential);
```

- `solveTISE` (`tise.cpp:1072-1174`) calls it instead of inlining; everything after (continuum construction, its own file-writing) is unchanged.
- `tise_solver_main.cpp` (`:160-168`) calls the same function instead of `buildUniformRadialGrid`+hardcoded `{1}`.
- **Critical follow-through**: every downstream call in `tise_solver_main.cpp` currently relying on the nullopt→classic-`{1}` default must receive the real `fillDropSet`/`nBSplines` explicitly, or a non-classic drop-set silently misattributes coefficients (the exact bug class fixed earlier this session, at a different call site). Checklist: `fillBandedMatrices`, `buildContinuumState`, `matchAsymptotic`, `writeContinuumInfo`, the `eigenstateCoefficients`/`eigenstate_NNN.dat` loop, and `writeEigenvectors` (confirm/add a `dropSet` parameter if it doesn't already have one).
- Surface `rightEdgeSingular` through `tise_solver_main.cpp`'s existing `addWarning`/`warnings.json` path.
- `minInterNodeGap` (G11): switch `computeEAcc`'s node-spacing input from the flat uniform formula to `minInterNodeGap(sgr.grid)` in both orchestration paths.
- Add ADR-0009 documenting this as resolving ADR-0008's stated revisit trigger; mark ADR-0008 superseded (not deleted).
- New tests: real-subprocess Step-join (`finite_square_well.yaml`) and interior-Singular-join (new config) tests are the ones that actually prove `tise_solver_main.cpp` wires this correctly end-to-end (Part G items 1-2) — these are what would have caught a missed threading site.

## Part B: Wire the remaining orphaned diagnostics (closes G9, G10)

- `classifyBoundStates`: call in `tise_solver_main.cpp` post-diagonalization, threshold `0.0`. Surface count via a `warnings.json` "physics" informational entry. Do **not** change `eigenvalues.dat`'s file format.
- `checkWellContainment`: call per bound state; on `notWellContained`, one `warnings.json` "physics" entry naming the state index.
- Case-3 remediation: wire `evaluateWindowedPotential` into `fillBandedMatrices`'s potential evaluation when `classifyAsymptote` returns `Irregular`, using its `recommendedTransitionWidth`. Real numerical-behavior risk — needs its own test (irregular-tail config, finite eigenvalues + taper actually applied) before enabling unconditionally.

## Part C: Correctness fixes (closes G7, G8, G12)

- G7: work out the correct one-sided finite-difference sign at the last `dDeltaDE` grid point against the free-particle closed form (δ≡0 ⇒ dδ/dE≡0 everywhere); fix; regression test for continuity with the neighboring point.
- G8: validate potential pieces don't overlap (raise on detected overlap), replacing the silent map-order tiebreak.
- G12: `tise_solver_main.cpp` reads `physics.mass`/`physics.hbar` if present and raises a clear config error if either isn't `1.0` — guard-rail, not full generalization.

## Part D: Implement `visualization.phase_shifts` plotting (closes G3)

Mirror the `plot_eigenstates`/`visualization.eigenstates` pattern from earlier this session: `plot_phase_shifts(tise_data, tise_dir)` in `analysis.py`, gated by `cfg.get("visualization", {}).get("phase_shifts", False)` in `run()`. TDD as before; free-particle config (δ≈0 flat line) as the real-subprocess sanity check.

## Part E: Documentation (closes G2, G15, G16)

Rewrite `README.md`/`TISE/README.md` to document the real pipeline (build `tise_solver`, `config.yaml` schema, `controller.py --config ...`, output/plot locations) while keeping `H-BoundStates`'s own docs, clearly distinguished. Update SDD §7.1's interface table + §5.1.3's stale pseudocode note. Fix the stale `tise.hpp:299` comment.

## Part F: CI (closes G6)

`.github/workflows/ci.yml`: apt-install LAPACK/BLAS/yaml-cpp/muparser/nlohmann-json, `cmake`+`ctest` for C++, Python venv + `pytest tests/`. Add `matplotlib` to `requirements.txt`. Gate on push/PR to default branch.

## Part G: Integration/e2e test plan (closes G4, G5, G13, G14; regression backstop for Part A)

1. Load `tests/{free_particle,finite_square_well,harmonic_oscillator,hydrogen}.yaml` by name (not inline-mutated equivalents) in real-subprocess tests — closes G4, and `finite_square_well.yaml` becomes the Step-join/strategic-grid regression test for Part A.
2. New interior-Singular-join config + real-subprocess test (mirrors `SolveTISETest.InteriorSingularityTriggersBSplineRemoval`'s `1/(x-20)` construction) — exercises the A4b removal path Part A wires into `tise_solver_main.cpp`.
3. Full-chain plot-existence assertions (closes G5): after a real `controller.py` run with `run_analysis: true`, assert `continuum_*.png`/`eigenstate_*.png` exist with nonzero size.
4. Explicit `run.run_analysis: false` real-subprocess test (closes G14).
5. `warnings.json` content tests for each warning type reachable from the real binary (E_acc, pole-proximity, Part A/B's new ones), each via a config engineered to trigger it.
6. `buildEnergyGrid` direct unit test (closes G13).
7. Full existing suite (146 C++ / 166 Python tests as of this session) re-run to confirm zero regressions.

## Verification

1. `cmake --build TISE/build && ctest --test-dir TISE/build` — all pass.
2. `pytest tests/` — all pass; confirm the finite-square-well/interior-singularity tests actually discriminate (would fail against a deliberately-reintroduced classic-dropset regression).
3. Manually run all 4 known-solution configs via `controller.py`, spot-check `warnings.json` for the new diagnostic entries.
4. Push to a branch, confirm the new GitHub Actions workflow runs and passes.

## Implementation status (post-execution)

All of Parts A-G landed. Final regression state: all 4 C++ suites green (`ctest --test-dir TISE/build`: BSplineTests, UtilsTests, TISETests, TimeEvolutionTests — TISETests alone now carries 165 individual `TEST`/`TEST_F` cases across 49 test suites), and all 185 Python tests green (`pytest tests/`). TDD (RED-then-GREEN) was followed for every new capability; every regression risk surfaced was caught by an existing or new automated test before being called done. Net working-tree diff for this phase (tracked files only): 15 files changed, 1194 insertions, 109 deletions, plus 8 new untracked files (this plan doc, ADR-0009, `.github/workflows/ci.yml`, and 3 new known-solution/regression YAML configs, `docs/bug-fixes/`).

**Part A — unify `solveTISE`/`tise_solver` grid+dropset construction: fully done.**
`buildStrategicGridAndDropSet` (new struct `StrategicGridResult` in `tise.hpp`) extracted exactly as scoped; `solveTISE` and `tise_solver_main.cpp` both call it. Every downstream call site in `tise_solver_main.cpp` was rethreaded with the real drop-set/B-spline-count instead of the classic-`{1}` default: `fillBandedMatrices`, `buildContinuumState`, `matchAsymptotic`, `writeContinuumInfo`, the `eigenstate_NNN.dat` loop, and `writeEigenvectors`. `minInterNodeGap(sgr.grid)` now feeds `computeEAcc` in place of the flat uniform-spacing formula. ADR-0009 written; ADR-0008 marked Superseded (not deleted).
*Deviation from plan*: the plan's checklist treated "the drop-set" as one thing to thread through. Implementation surfaced a real distinction the plan didn't anticipate — `writeEigenvectors`/`eigenstateCoefficients` need the drop-set used for *zero-padding a full eigenvector* (`fillDropSet ∪ {nBSplines}`, since `er.vectors` truncates B_N separately from `fillDropSet`), which is a strict superset of the `fillDropSet` used everywhere else. This was caught by reasoning through the code before running anything (not by a failing test), and fixed with a locally-constructed `fullDropSet` in `tise_solver_main.cpp` scoped to just those two call sites.

**Part B — wire orphaned diagnostics: fully done.**
`classifyBoundStates` runs post-diagonalization and surfaces the bound-state count as a `warnings.json` "physics" entry; `eigenvalues.dat`'s format is unchanged as scoped. `checkWellContainment` runs per state and warns on `notWellContained`. Case-3 remediation (`evaluateWindowedPotential`) is wired into `fillBandedMatrices` via new optional `case3RightR`/`case3RightDelta` params, covered by `FillBandedMatricesCase3WindowTest` (2 cases) plus a real-subprocess `case3_irregular_tail.yaml` test.
*Deviation from plan*: `checkWellContainment` had to be restricted to states with `er.values[j-1] < 0.0` (genuinely bound), which the plan didn't specify. An unrestricted call against the real `config.yaml` produced 40+ warnings — almost every continuum-like pseudostate "collides with the wall," which is physically expected and not meaningful. Restricting to bound-only states cut this to exactly the physically meaningful cases (3, for the real config).

**Part C — correctness fixes: fully done, one finding inverted the plan's assumption.**
G8 (overlapping potential pieces): `validateNoOverlappingPotentialPieces` added and called at parse time in `tise_solver_main.cpp`, raising `std::runtime_error` on overlap; 6 new tests. G12 (mass/hbar guard-rail): `tise_solver_main.cpp` now raises a config error if either is set to anything but `1.0`; covered by `test_non_unity_mass_raises_and_leaves_no_partial_files`.
*Deviation from plan*: G7 (`dDeltaDE` sign at the last energy-grid point) — the plan assumed a sign bug existed and asked to "work out the correct sign... fix." Investigation found **no bug**: an initial verification attempt showed a same-sign ~10% magnitude discrepancy, traced (via a throwaway eigenvalue-printing helper) to the square well's densely-packed pole spectrum making δ(E) genuinely non-smooth in that window, not a formula error. `MatchAsymptoticDDeltaDESignTest` was added instead, using a narrow pole-avoiding grid to prove the existing formula is correct, and the stale `// TODO: will this sign be wrong` comment was replaced with an explanation referencing that test. Net effect matches the plan's intent (the ambiguity is resolved and tested) but the deliverable is a proof + comment fix rather than a code change.

**Part D — `visualization.phase_shifts` plotting: fully done** as scoped. `plot_phase_shifts` added to `analysis.py`, gated identically to `plot_eigenstates`; `TestPlotPhaseShifts` (2 unit tests) plus 2 `TestRun` gating tests plus real-subprocess coverage via `free_particle.yaml`/`finite_square_well.yaml` (both now set `visualization.phase_shifts: true`).

**Part E — documentation: fully done** as scoped. `README.md` and `TISE/README.md` rewritten to document the real `config.yaml`→`controller.py`→`tise_solver`→`analysis.py` pipeline (keeping `H-BoundStates`'s docs, clearly distinguished); `docs/SDD.md` §7.1 gained the `controller.py` interface row, §5.1.3 gained the illustrative-pseudocode note; the stale `tise.hpp:299` A4b comment was fixed.

**Part F — CI: fully done, with one caveat.** `.github/workflows/ci.yml` added: installs LAPACK/BLAS/Eigen3/nlohmann-json/muparser/yaml-cpp/gtest via apt, builds+`ctest`s the C++ suite, then sets up Python 3.12 and runs `pytest tests/`. `matplotlib` added to `requirements.txt`. Package names were cross-validated against this machine's actual installed packages (`dpkg -l`/`dpkg -L libgtest-dev`), not against a live GitHub Actions run.
*Not yet done*: the plan's own Verification item 4 ("push to a branch and confirm the workflow actually runs and passes") has **not** been executed — pushing is a shared/visible action requiring explicit user approval, which hasn't been sought yet. This is the one open item from the whole plan.

**Part G — integration/e2e test plan: fully done.** All 7 items landed: the 4 known-solution configs are now loaded by name (`TestKnownSolutionConfigFilesLoadAndRun`) rather than being orphaned; a new `interior_singularity.yaml` + real-subprocess test exercises the A4b removal path through the actual binary; `test_end_to_end_run_with_tise_and_analysis_enabled` now asserts `continuum_1.png`/`eigenstate_1.png` exist with nonzero size; `test_run_analysis_false_produces_no_analysis_output` covers the `run_analysis: false` path explicitly; `warnings.json` content is now checked per warning type via dedicated real-subprocess test classes (`TestFiniteSquareWellRealSubprocess`, `TestHydrogenConfigFileRealSubprocess`, `TestInteriorSingularityRealSubprocess`, `TestRightEdgeSingularWarningRealSubprocess`, `TestBoundStateDiagnosticsRealSubprocess`, `TestCase3RemediationRealSubprocess`); `BuildEnergyGridTest` adds 5 direct unit tests; the full suite re-run (item 7) is the 165/185-passing state reported above.

**Verification checklist**: items 1-3 done and confirmed (`ctest`/`pytest` both fully green; all 4 known-solution configs manually run via `controller.py` and `warnings.json` spot-checked for the new diagnostic entries earlier in this session). Item 4 (push + confirm CI run) is the one remaining action, deferred pending explicit user approval to push.

**Process note**: no git commits were made at any point during this phase — every change above is an uncommitted working-tree modification, left for the user to review/stage/commit on their own terms.
