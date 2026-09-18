# TISE-Generalization Branch Review

**Reviewed:** 2026-09-11
**Branch:** `TISE-Generalization` @ `b9f3200` vs. `main` @ merge-base `d01efad`
**Scope:** 52 commits, 146 files changed, ~18,649 insertions / 266 deletions. TDSE (time-dependent) work is intentionally out of scope for this review — it is not implemented on this branch, and is consistently documented as future work (`docs/adr/0006-tdse-propagator-choice.md`, `docs/planning/tdse-task-breakdown.md`).

## 1. Summary / Verdict

**Recommend merge, conditional on a small set of cheap documentation fixes (Section 5.1).** No correctness defects were found in the TISE solver itself as of HEAD. The branch takes the TISE module from a hydrogen-only demo to a general, config-driven 1D solver (strategic node placement, unified grid/drop-set construction, Coulomb-tail continuum matching, a warnings sidecar) with real end-to-end test coverage and three independently-verified rounds of physics validation against closed-form solutions. The process that produced this branch is visibly working as intended: two real regressions (a 41%-high eigenvalue bug, a `dδ/dE` sign/branch bug) were caught by the project's own verification-report discipline and fixed before this review, not found here.

The issues this review did find are entirely in the documentation layer, not the code: several docs (`SDD.md`, both `README.md`s, two C++ comments) still describe the *old*, pre-generalization behavior for capabilities this branch actually implements — most notably, three separate places still claim Coulomb-tail continuum matching "is not implemented," when it is (ADR-0013). One file (`analysis.py`) contradicts itself outright: its own module docstring says "No plotting... deliberately deferred per ADR-0005," directly above code that plots and saves three kinds of files. These are cheap to fix and worth fixing before merge because they're actively misleading to the next person who reads them — but none of them reflect a defect in the software's behavior.

Separately, this review also resolves the specific open question the team flagged going in: `docs/planning/engineer-a-plan-cleanup.md` item 8 ("verify node placement is at least as capable as what we discussed in original design meetings") has apparently just been investigated — two new ADRs (`0014`, `0015`) sitting **uncommitted** in the working tree formalize the answer (two real, narrow gaps, both reasonably deferred). See Section 5.1(c).

## 2. Scope of This Review

- **Compared:** `main` (last commit 2026-07-26) vs. `TISE-Generalization` (HEAD `b9f3200`, 2026-09-11), merge-base `d01efad`.
- **Verified directly** (not just agent-reported): ran the full C++ and Python test suites myself; read `SDD.md`, both `README.md`s, `analysis.py`, `engineer-a-plan-cleanup.md`, `docs/adr/README.md`, and the two uncommitted ADRs in full; spot-checked commit messages against the release-readiness plan's claims.
- **Out of scope:** TDSE propagation (not implemented, not expected to be — confirmed guarded, see Section 3.4).

## 3. Implementation & Test Coverage

### 3.1 What the branch adds

The TISE module goes from a hardcoded-hydrogen demo to a config-driven general solver:

- **`tise::solveTISE` / `SolveTISEResult`** — a real orchestration entry point returning `{eigen, bs, grid, nBSplines, dropSet}`, so both `H-BoundStates` and the production `tise_solver` binary consume the exact basis/grid/drop-set that was actually diagonalized instead of independently re-deriving (and risking diverging) copies of it.
- **Strategic node placement + generalized drop-set** — `buildStrategicGridAndDropSet` auto-detects potential structure at each piece-to-piece join (`Step`, `StitchedKink`, interior `Singular`) and either inserts degenerate knots or removes touching B-splines from an arbitrary index set, replacing the old hardcoded `{1}` drop-set. Unifies behavior between `H-BoundStates` and `tise_solver` (ADR-0009, supersedes ADR-0008 once the old split was found to *silently degrade eigenvalues*, not just duplicate code).
- **Classifier functions** (`classifyAsymptote`, `classifySequenceConvergence`, `isSingularApproaching`, `detectPotentialStructure`) — magic numbers extracted to named/optional parameters; a guard added against too-small `numSamples`/`windowSize`. `docs/planning/tise-classifier-configurability-exposure.md` catalogs ~20 such parameters and recommends which are worth exposing to `config.yaml` later (none are, yet — see Section 5.2).
- **Coulomb-tail continuum matching** (ADR-0013, supersedes deferred ADR-0010) — a hand-rolled Coulomb wave-function evaluator (Lanczos Γ, Numerov-integrated F_ℓ/G_ℓ, WKB start), dispatched from `matchAsymptotic` via a new `tise.continuum.l` config field. Reduces hydrogen phase-shift error from ~0.02–0.07 rad (flat-asymptote formula misapplied to a Coulomb tail) to ~0.001–0.002 rad. An earlier RK4-based version of this evaluator hit a 60-second timeout on the production energy grid; switching to Numerov fixed it.
- **R-bound / domain-validity warning** — warns (via `warnings.json`) when `|V(rMax)|` is non-negligible relative to the smallest requested continuum energy, i.e., when the box isn't big enough for the flat-asymptote assumption to hold. A practical, cheap mitigation for the case ADR-0011 (CAP/ECS) leaves otherwise unaddressed.

### 3.2 Test inventory

**Python** (`1D-QM-Playground/tests/`): `test_analysis_unit.py` (1172 lines), `test_analysis_integration.py` (1100 lines), `test_controller_unit.py` (782, mocked subprocess), `test_controller_integration.py` (373, real subprocess) — plus 7 YAML known-solution configs (`free_particle`, `finite_square_well`, `harmonic_oscillator`, `hydrogen`, `interior_singularity`, `right_edge_singularity`, `case3_irregular_tail`) and a session-scoped `conftest.py` fixture that builds the real `tise_solver` binary once per test run.

**C++** (`TISE/tests/test_tise.cpp`): 128 `TEST()` cases across suites including `ClassifyAsymptoteTest` (9), `SolveTISETest` (7), `ClassifySequenceConvergenceTest` (7), `BuildStrategicGridAndDropSetTest` (6), `WarnIfContinuumExceedsEAccTest` (4). Plus `test_bspline.cpp`, `test_utils.cpp`, `test_time_evolution.cpp` via GoogleTest/CMake (`BUILD_TESTING=ON`).

**I ran both suites myself**, the same sequence CI uses:

```
cmake -S TISE -B TISE/build -DBUILD_TESTING=ON && cmake --build TISE/build -j
ctest --test-dir TISE/build --output-on-failure
→ 100% tests passed, 4/4 (BSplineTests, UtilsTests, TISETests [113s], TimeEvolutionTests)

python3.10 -m pytest tests/ -q
→ 207 passed, 4 benign matplotlib UserWarnings (zero-width xlim on a degenerate plot)
```

`.github/workflows/ci.yml` (added in `8236239`) runs this identical sequence (LAPACK/BLAS/Eigen/nlohmann-json/muparser/yaml-cpp/gtest install → C++ build+ctest → `pip install` → `pytest tests/`), so the branch is very likely green in CI as well.

**End-to-end coverage is real, not mocked.** `test_controller_integration.py::TestRunEndToEnd::test_run_end_to_end_success` drives the complete `config.yaml → controller.py → tise_solver (real subprocess) → analysis.py` flow. `test_analysis_integration.py` has dedicated real-subprocess test classes per known-solution system (hydrogen, finite-square-well, interior-singularity, right-edge-singularity, Case-3, R-validity warning, CLI) plus round-trip tests with continuum on/off — all built on the real compiled binary, not a stub.

### 3.3 Physics validation evidence

Three successive rounds of verification reports (`docs/tests/reports/{9d0f04b,8236239,f4e8359}/`) compare numerical output against closed-form analytic references across all 7 known-solution systems — not smoke tests. Eigenvalues match harmonic-oscillator/infinite-and-finite-well/hydrogen (`E_n = -1/(2n²)`) closed forms to 1.4e-12–4e-8; hydrogen continuum phase shifts are cross-checked against `mpmath.coulombf`/`coulombg` (an independent reference implementation, not the project's own evaluator) to ≤6e-3 rad.

This process caught real bugs, which is the strongest evidence it's working: round 2 found `interior_singularity` **41% too high** on the field-free side (over-aggressive B-spline removal at an interior singular join) — fixed in `0bf17be`, then independently re-verified as fixed (to 1e-11–1e-13) in round 3 (`f4e8359`). The same fix pass also corrected a `dδ/dE` branch/sign bug and gated the Case-3 tail taper on `continuum.enabled` (it had previously shifted bound-state energies even when continuum output was off). All 7 known-solution configs now **PASS** as of the round-3 report.

### 3.4 TDSE guard

`controller.py` raises `ControllerError` if `run.run_tdse` is set (lines ~403-404); `tise_solver` never calls `runTimeEvolution`. TDSE is properly fenced off, not silently no-op'd. The one caveat is that `H-BoundStates` (the standalone demo binary, not the production `tise_solver` path) runs time evolution unconditionally after every solve — this is explicitly documented as a known, out-of-scope limitation in `TISE/README.md`'s "Known Limitations" section, not a gap this review needs to flag.

## 4. SDD & ADR Changes: What Changed and Why

### 4.1 `docs/SDD.md` diff (55 lines)

Mostly a post-hoc documentation catch-up, not upfront design — driven by commits `bbf80ec`, `7249dfc`, and `f4e8359`:

- §1.2/§1.5: adds ADR-0005 to the deferred-extension list; records that the referenced PHY5606/Bachau PDFs are now checked into `docs/references/` (this is what actually closed ADR-0012).
- §3.3 traceability table: adds an ADR-0005 row.
- §5.2.3 (TISE internal design): rewrites the "current implementation baseline" to reflect that `tise::` is now config-driven, not hardcoded hydrogenic; adds a paragraph describing `SolveTISEResult`, automatic strategic node placement, and the interior-vs-edge singularity distinction.
- §6.1 schema: adds `tise.continuum.l` and `visualization.bound_states_squared`; clarifies `error_threshold` is unconsumed by `tise_solver`.
- §6.3: adds `eigenstate_NNN.dat` and `warnings.json` rows; clarifies `hamiltonian.dat`/`overlap.dat` as LAPACK banded storage (previously left undecided).
- §7.1/§7.2.2/§7.2.3: adds `controller.py --config`; adds `eigenstate_NNN.dat` to Analysis inputs.
- §11.1/§12.B: records muparser as the adopted expression parser; points the CAP/Siegert/ECS discussion at the new ADR-0011.
- §8 (warnings taxonomy, previously miscited internally as "§12.D"): documents the `data/tise/warnings.json` sidecar contract.

### 4.2 ADR inventory

`main` has ADRs 0001–0004 (simple "defer-X" acceptances: FEDVR basis, WKB collocation, multi-particle/3D, visualization plot parameters). This branch adds:

| ADR | Status | What it decided |
|---|---|---|
| 0005 | Accepted (deferred) | Analysis produces no output artifact until REQ-F-060 is computable (Phase 7/8) |
| 0006 | Accepted (deferred) | Symmetric split-operator chosen for TDSE; Crank-Nicolson family deferred |
| 0007 | Accepted (deferred) | `tise_solver` writes all eigenvalues/vectors; bound-state filtering left to downstream consumers |
| 0008 | **Superseded by ADR-0009** | Originally accepted the `solveTISE`/`tise_solver` orchestration split as a known gap; superseded once that gap was found to cause silent eigenvalue degradation |
| 0009 | Accepted (implemented) | Unifies grid/drop-set construction (`buildStrategicGridAndDropSet`) between `solveTISE` and `tise_solver_main.cpp` |
| 0010 | **Superseded by ADR-0013** | Deferred Coulomb-tail continuum matching (flat-asymptote formula only) |
| 0011 | Accepted (deferred) | CAP/Siegert/ECS boundary conditions remain unimplemented; solver stays real-valued (`DSBGV`)-only |
| 0012 | Accepted (implemented) — verified against source | Confirms the A4 knot-multiplicity/continuity-order formula against `PHY5606_F25_Bsplines_v2.pdf`; no discrepancy |
| 0013 | Accepted (implemented) | Implements Coulomb-tail continuum matching; supersedes ADR-0010 |
| 0014* | Accepted (deferred) | *Uncommitted — see 4.4.* Defers delta-potential join detection in strategic node placement |
| 0015* | Accepted (deferred) | *Uncommitted — see 4.4.* Defers a user-supplied arbitrary node-placement formula |

`docs/adr/README.md`'s index table is accurate and correctly reflects the 0008→0009 and 0010→0013 supersession chain. One minor readability note: commit `f4e8359`'s title ("...close ADR-0012, add R-bound validity warning") doesn't name ADR-0010/0013 even though the commit's real headline change is the Coulomb-tail implementation — the commit body is correct (it does credit both), so this is a git-log-readability nit, not a documentation error.

### 4.3 Bug-fix doc cross-check

`docs/bug-fixes/0001-continuum-state-coefficient-basis.md` records the free-particle continuum coefficient-basis bug: introduced `9c318c1` (2026-08-12), fixed `9d0f04b` (2026-08-30), with regression tests in both `TISE/tests/test_tise.cpp` and `tests/test_analysis_integration.py`. Consistent with commit history; no inconsistency found. (Note: this is a *fixed* date of 2026-08-30, not 2026-08-28 — a minor discrepancy from an earlier informal note, not a doc defect.)

### 4.4 The uncommitted ADR-0014/0015 and cleanup item 8

`docs/planning/engineer-a-plan-cleanup.md` item 8 reads: *"Verify node placement implementation is at least as capable as what we discussed in original design meetings."* As of the last commit, this item is unstruck — genuinely open.

The working tree currently contains **uncommitted** files, all dated 2026-09-11 (today): `docs/adr/0014-defer-delta-potential-join-detection.md`, `docs/adr/0015-defer-user-supplied-node-placement-formula.md`, and a matching edit to `docs/adr/README.md`. These are a capability review against the original design docs (`docs/planning/architecture-06-20.md`'s 2026-07-03 stakeholder-feedback notes), explicitly prompted by cleanup item 8, and they answer it:

- **ADR-0014**: the original design's required-treatment table for strategic node placement has **four** join types (`architecture-06-20.md`, `tise-task-breakdown.md`'s A4 spec), including a delta-potential row (`δ(x−x₀)` → pile up degenerate knots, `extra = order-2`). The implementation's `JoinType` enum (`TISE/tise.hpp:521-527`) has only `Continuous`/`Step`/`StitchedKink`/`Singular` — no `Delta`. The math for this case was actually already derived and independently cross-checked (`engineer-a-plan-A4.md:381`, ADR-0012's verification), just never wired up — and `config.yaml`'s muparser-based potential DSL can't express a literal delta function anyway, so there's nothing to detect even if it were added. Formally deferred, with a revisit trigger tied to the config schema someday gaining distributional-potential support.
- **ADR-0015**: the 2026-07-03 stakeholder notes also said the program should "leave open the option for the user to supply an arbitrary sequence of points via a formula x(n)." This half of the request was never carried into REQ-F-050, `tise-task-breakdown.md`, or `config.yaml`'s `bspline` block (which exposes only `n_nodes`, `order`, `domain`) — and, unlike WKB-density (ADR-0002) or CAP/ECS (ADR-0011), it was never formally deferred either; it simply stopped appearing in any document after the initial notes. ADR-0015 formalizes the deferral now, with a revisit trigger.

Both gaps are narrow, low-impact, and reasonably deferred rather than blocking — the review's independent judgment agrees with what these ADRs conclude. **Recommendation: commit `0014`, `0015`, and the `docs/adr/README.md` update as part of finishing this branch, and strike item 8 in `engineer-a-plan-cleanup.md`** (see Section 6, item 3).

### 4.5 README changes

Both `README.md` (83 lines) and `TISE/README.md` (150 lines) were substantially rewritten to describe the real `config.yaml → controller.py → tise_solver → analysis.py` pipeline, replacing descriptions of the old ad-hoc hydrogen-only build. Genuinely useful additions: a "Quick start," a table of the 7 known-solution test configs, a clear split between the two binaries (`tise_solver` = production pipeline; `H-BoundStates` = standalone demo, potential now passed as JSON via argv), documentation of the potential DSL, and a coverage-gate section (`ENABLE_COVERAGE`, ≥80% on `tise.cpp`). Both are strong overall — the one defect is the stale Coulomb-tail claim covered in Section 5.1(a).

## 5. Gaps Identified

### 5.1 Should fix before merge (cheap, high-value doc-correctness fixes)

**(a) Stale "Coulomb-tail not implemented" claims, in three docs plus two code comments.** All were written before ADR-0013 and never updated:

- `docs/SDD.md:564` — *"The Coulomb-tail sub-branch... is required by REQ-F-030 but not yet worked out at this level of detail by any source document... This gap is formalized as [ADR-0010]..."* Note this is a narrative-prose staleness specifically — §6.1's schema table in the same file *does* correctly cite `tise.continuum.l` under ADR-0013, so this file is internally inconsistent, not just behind the code.
- `README.md:80` — *"Not yet implemented: ... the Coulomb-tail continuum-matching formula for potentials with a genuine `1/r`-type tail beyond the box..."*
- `TISE/README.md:273` — *"...the Coulomb-tail continuum-matching formula... is not implemented at all, and implementing it is not currently planned — see `docs/adr/0010-defer-coulomb-tail-continuum-matching.md`..."*
- `TISE/tise_solver_main.cpp:472` — comment: *"...which needs its own Coulomb-aware matching (see ADR-0010) rather than a bigger box."* Should cite ADR-0013, where that matching now lives.
- `TISE/tise_solver_main.cpp:268` — comment: *"{l, eta} for matchAsymptotic's Coulomb branch (Coulomb-tail continuum matching, ADR-0009 supersedes ADR-0010)"* — this cites the **wrong superseding ADR**; ADR-0009 supersedes ADR-0008 (grid/drop-set unification, an unrelated topic). Should read "ADR-0013 supersedes ADR-0010."

All five should be updated to reflect ADR-0013. This is the single clearest fix in this review — it's the branch's headline feature, described as absent in its own primary docs.

**(b) `ADR-0005`/`analysis.py`/`SDD.md §7.2.3` self-contradict on plotting.** `ADR-0005` and `SDD.md` (§3.3 row, §7.2.3) both assert Analysis produces no output artifact of any kind. `analysis.py`'s own module docstring (lines 15, 25, 29-31) repeats this: *"No plotting"* / *"No output artifact of any kind: no plot files, no derived-data files, no data/analysis/ directory... deliberately deferred per ADR-0005."* Yet the same file calls `plt.savefig(...)` three times (lines 464, 519, 544), writing `continuum_{idx}.png`, `eigenstate_{idx}.png`, and `phase_shifts.png` into `data/tise/` — implemented in commits `8236239`/`5cd28e4`. `docs/bug-fixes/0001-...md` independently references these `.png` files as existing artifacts used for visual verification, and the untracked `.dat` files in the working tree (Section 5.3) are further evidence of a real run producing this output.

This isn't just stale — it's a currently-false statement sitting in a module's own top-of-file docstring. Two ways to resolve it, either is acceptable, but one is needed:
1. Mark ADR-0005 superseded by a new ADR that documents the actual (now-implemented) plot-artifact behavior — following the same pattern ADR-0009/ADR-0013 already established for other superseded deferrals — and update `analysis.py`'s docstring + `SDD.md` §3.3/§7.2.3 to match.
2. Or, if the plotting in `analysis.py` was meant to stay experimental/non-canonical rather than a committed-to artifact format, say so explicitly in the docstring instead of flatly denying it exists.

**(c) Commit the uncommitted ADR-0014/ADR-0015 + `docs/adr/README.md` update**, and strike item 8 in `engineer-a-plan-cleanup.md` (see Section 4.4). This closes out a specifically-flagged open item with work that already exists and just needs to land.

**(d) Add the R-bound/domain-validity warning to `SDD.md`'s warning taxonomy** (§8, line ~1041). Currently lists only three physics warnings (Case-3 boundary discontinuity, `E_max > E_acc`, well-containment) — the `|V(rMax)|`-vs-`E_min` check added in `f4e8359` (Section 3.1) isn't mentioned. Minor, but a real omission.

### 5.2 Worth an explicit call, not a hard blocker

**Continuum/phase-shift output is still written for `interior_singularity`'s split-domain case, where it is physically meaningless.** Directly from the round-3 verification report (`docs/tests/reports/f4e8359/interior_singularity.md:6,41`): *"One item remains open: continuum states and phase shifts are still written for this split domain, where they have no meaning,"* with `phase_shifts.dat` containing numbers (δ = 0.92, −2.02, −4.39...) the report itself calls "flat-asymptote matching... applied to a state built from two decoupled boxes... describe[s] nothing physical." The same class of bug — silently producing bad output with no gate — was treated as worth fixing when it showed up for the right-edge-singular case (`0bf17be`'s `right_edge_singularity` fix: continuum construction is now *refused*, not just warned about). The interior-singularity analogue of that same fix was explicitly recommended ("iteration-2 recommendation 4: refuse or warn for split domains") and explicitly **not yet implemented**.

This is pre-existing behavior, not a regression introduced by this branch, and `interior_singularity.yaml` is a test config, not obviously a production-representative case — so it's reasonable to treat this as a fast-follow rather than a hard blocker. But because it's silent-bad-output on a config the test suite directly exercises, and the project's own review discipline has already flagged and fixed the identical failure mode once this branch, it deserves an explicit yes/no from reviewers rather than being silently left off the list. Minimum fix: emit a `warnings.json` caveat when continuum output is requested on a split (interior-singular) domain. Fuller fix: match the H5/right-edge-singular treatment and refuse continuum construction outright.

### 5.3 Safely deferred / no action needed

All of the following are quantified, low-impact, already tracked in their own planning docs with explicit priority, or already resolved and independently re-verified — no action needed for this merge:

- Coulomb-function evaluator's ~2e-3 rad accuracy floor (Wronskian −0.9982 vs. −1, O(h²) not O(h⁴) derivative, WKB start) — localized to the hand-rolled F/G evaluator, doesn't affect bound states; `docs/tests/reports/f4e8359/README.md` already documents the improvement path (Wronskian reporting → O(h⁴) derivative → DLMF 33.11 asymptotic series or Steed's continued fractions, expected gain to ≲1e-6).
- `dδ/dE`'s sensitivity to fine-step size (amplifies basis error by 1/(2·fineDE)) — same category, a numerics-polish item with a known fix (Richardson extrapolation or a larger step).
- Case-3 classifier's "will be tapered" warning text contradicting the orchestrator's "not tapered" message when continuum is disabled — a warning-wording bug, not a physics bug.
- `computeEAcc` ignoring `V_min` — conservative-direction inaccuracy in a warning heuristic, not silently-wrong output.
- `analysis.py` plotting every basis-limited state (not just physically-accurate ones) and plotted eigenstates having an arbitrary sign — cosmetic/UX.
- Missing analytic-value assertions on the free-particle and harmonic-oscillator e2e tests — both systems already PASS per the independent verification reports; this is test-hardening, not a known defect.
- `tise-classifier-configurability-exposure.md`'s two candidate config-exposure items (`edgeTolerance`, `coulombExponentTol`/`transitionWidthFraction`) — explicitly scoped as low/medium priority, "described only — not implemented," and nothing in the current `config.yaml` needs them yet.
- `engineer-a-plan-cleanup.md` item 3 ("review PR bot comments") — process-only, unverifiable from the repo alone, not a code-correctness question.
- CI "push and confirm it actually runs" verification — the workflow file exists and this review reproduced its exact command sequence locally with a clean pass; confirming a green Actions run is a one-click check for reviewers, not a code change.

### 5.4 Housekeeping

The working tree has three untracked files — `TISE/continuum_state_001.dat`, `continuum_state_002.dat`, `phase_shifts.dat` — confirmed as plain-text output from a manual run (not build artifacts, not already covered by `.gitignore`). The repo's `.gitignore` also has a pre-existing typo (`EigennState_*`, double "n") that doesn't match the actual `EigenState_XXX` naming pattern either. Recommend adding a `*.dat`-under-`TISE/`-style pattern (or similar) and fixing the typo while touching this file. Trivial, non-blocking.

## 6. Recommendations / Implementation Plan

Ordered checklist for the "should fix" items (Section 5.1), each independently small:

1. **Fix the 5 stale Coulomb-tail references** to cite ADR-0013 instead of (or in addition to, where historical context is useful) ADR-0010:
   - `docs/SDD.md:564` — rewrite the sentence to state the Coulomb-tail sub-branch *is* implemented (ADR-0013), removing the "not yet worked out... gap formalized as ADR-0010" framing.
   - `README.md:80` — remove the Coulomb-tail clause from the "Not yet implemented" list (or replace it with an accurate residual caveat, e.g. the accuracy-floor note from Section 5.3, if the team wants to keep a caveat there).
   - `TISE/README.md:273` — same: replace "is not implemented at all, and implementing it is not currently planned" with a pointer to ADR-0013 and, if desired, the accuracy-floor caveat.
   - `TISE/tise_solver_main.cpp:472` — change "(see ADR-0010)" to "(see ADR-0013)".
   - `TISE/tise_solver_main.cpp:268` — change "ADR-0009 supersedes ADR-0010" to "ADR-0013 supersedes ADR-0010".
2. **Resolve the ADR-0005/`analysis.py` plotting contradiction** — pick one of the two options in Section 5.1(b) and apply it consistently across `docs/adr/0005-*.md` (or a new superseding ADR), `analysis.py`'s module docstring (lines 15, 25, 29-33), and `docs/SDD.md` §3.3/§7.2.3.
3. **Commit `docs/adr/0014-defer-delta-potential-join-detection.md`, `docs/adr/0015-defer-user-supplied-node-placement-formula.md`, and the current uncommitted edit to `docs/adr/README.md`**; strike item 8 in `docs/planning/engineer-a-plan-cleanup.md`.
4. **Add the R-bound/domain-validity warning** to `docs/SDD.md`'s warning taxonomy (§8, near line 1041), as a fourth physics-warning bullet alongside the existing three.

Reviewers should separately decide on Section 5.2 (interior-singularity continuum output) — recommend at minimum a one-line `warnings.json` addition before merge if the team wants the cheap version, or an explicit "accepted, tracked as follow-up" note in `tise-known-solution-followup-plan.md` if not.

Section 5.4 (`.gitignore`) can be folded into whichever of the above commits touches nearby files, or done separately — no urgency either way.

## 7. Appendix: Verification Commands Run

```
$ cmake -S TISE -B TISE/build -DBUILD_TESTING=ON
$ cmake --build TISE/build -j
  → clean build, no errors/warnings of note

$ ctest --test-dir TISE/build --output-on-failure
  → 100% tests passed, 4/4 (BSplineTests, UtilsTests, TISETests [113s], TimeEvolutionTests)

$ python3.10 -m pytest tests/ -q
  → 207 passed, 4 benign matplotlib UserWarnings (zero-width xlim on a degenerate plot)
```

Both sequences match `.github/workflows/ci.yml` exactly.
