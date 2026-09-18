# TISE Engineer-A Tasks + Cleanup — Implementation Plan (Master)

> **Status:** Master reference plan covering all of Engineer A's tasks (A1–A5) plus shared infrastructure (Section 0) and the cleanup tasks (§4 of `tise-task-breakdown.md`). We are implementing one task at a time, pausing for review after each. See [`engineer-a-plan-A1.md`](engineer-a-plan-A1.md) for the actively-in-progress slice (currently task A1). Sections for A2–A5 and cleanup remain here as reference for when we pick them up in later rounds.

## Context

`docs/SDD.md` and `docs/planning/tise-task-breakdown.md` define the TISE solver's remaining numerics work (SDD §10.2 Phase 4), split across two parallel workstreams. We are picking up **Engineer A's** five tasks (A1–A5: boundary-condition asymptote classification, bound/continuum classification, well-containment diagnostics, strategic node placement, continuum-accuracy warning) plus the four **cleanup tasks** (§4) not assigned to either engineer. Engineer B's continuum-state-construction tasks (B1–B5) are explicitly out of scope and untouched by this plan.

REQ-NF-010 requires ≥80% test coverage maintained on every change. **No coverage tooling exists anywhere in this repo today** (confirmed directly: no `--coverage`/gcov/lcov flags anywhere, no CI config). Before any task can be verified against that bar, we need a minimal local coverage-measurement mechanism — that's Section 0 below, not scope creep, since without it "≥80%" is unverifiable rather than just unenforced.

**Corrections to the source framing, confirmed by direct inspection during planning:**
- `BSpline-CPP/` (the actual project root, one level above `1D-QM-Playground/`) **is already a git repository**, branch `TISE-Generalization`, clean working tree, tracks `origin` (`git@github.com:Razgriz117/BSpline-CPP.git`), latest commit is the task-breakdown doc itself. No `git init` needed. Per standing instructions, commits are only created when explicitly requested — commit timing/granularity across these sections is a separate decision.
- The actual dev environment (checked directly, not just a planning sandbox) is **missing `GTest`, `muparser`, `nlohmann-json`, and `gcovr`/`lcov`** — the same packages the project's own README already documents as dependencies, just not yet installed here. Nothing in `TISE/` can currently be built or tested until these are installed (Section 0.1). Installing them needs `sudo`/`pip`, so expect a permission prompt when that step runs.

**Working agreement:** this file has five independent sections, one per task (A1–A5), plus a shared Section 0 (infra + cleanup). Each section is implemented one at a time — code + tests + a coverage check — with a pause for review before starting the next.

---

## Section 0 — Environment, Coverage Tooling, Baseline, Cleanup (do first, shared by everything below)

### 0.1 Install missing dependencies, confirm a clean baseline
```bash
sudo apt-get install libgtest-dev libmuparser-dev nlohmann-json3-dev liblapack-dev libblas-dev libeigen3-dev
pip install --user gcovr   # or: sudo apt-get install gcovr lcov
```
Then, **before touching any code**, confirm the existing suite builds and passes as-is:
```bash
cd TISE && cmake -S . -B build -DBUILD_TESTING=ON && cmake --build build -j && ctest --test-dir build --output-on-failure
```

### 0.2 Coverage tooling (`TISE/CMakeLists.txt`)
Insert before line 9 (`option(BUILD_TESTING ...)`) — must land before `add_subdirectory(tests)` at line 23, since directory-scoped compile/link flags only propagate to targets/subdirectories added *after* they're set:
```cmake
option(ENABLE_COVERAGE "Build with --coverage instrumentation (gcov/gcovr)" OFF)
if(ENABLE_COVERAGE)
    add_compile_options(--coverage -O0 -g)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()
```
(Using `CMAKE_EXE_LINKER_FLAGS`/`CMAKE_SHARED_LINKER_FLAGS` rather than `add_link_options`, since the latter needs CMake ≥3.13 and this repo declares `cmake_minimum_required(VERSION 3.10)`.)

Measurement workflow (repeat after every section below lands):
```bash
cmake -S . -B build-coverage -DBUILD_TESTING=ON -DENABLE_COVERAGE=ON
cmake --build build-coverage -j
find build-coverage -name '*.gcda' -delete
ctest --test-dir build-coverage --output-on-failure
mkdir -p build-coverage/coverage
gcovr --root . --filter '.*tise\.cpp' --exclude '.*/tests/.*' --exclude '.*main\.cpp' \
      --print-summary --html --html-details -o build-coverage/coverage/index.html build-coverage
```
`build-coverage/` matches the existing `.gitignore`'s `build*` pattern, already excluded from commits. We gate ≥80% on **`tise.cpp`'s own reported percentage** — essentially all new A1–A5 code lands there — not on `BSpline.cpp`/`time_evolution.cpp`'s pre-existing coverage, which this task set doesn't touch and won't retroactively backfill.

### 0.3 Baseline backfill
`inInterval` and `evaluateFunction` (`tise.cpp:48-111`) currently have **zero direct unit tests** (only indirect coverage via `fillBandedMatrices` tests) and aren't declared in `tise.hpp` at all. Promote both to header declarations (pure addition, no behavior change) and add:
- `InIntervalTest.{ClosedBothEnds, HalfOpenVariants, InfinityBounds, ThrowsOnMalformedString}`
- `EvaluateFunctionTest.{PicksCorrectPiece, ThrowsWhenXUncovered}`

This gives an honest `tise.cpp` coverage baseline before any new A-series code lands, so later measurements reflect A1–A5's own contribution cleanly rather than conflating it with a pre-existing gap.

### 0.4 Cleanup tasks 1–4
No dependency on A1–A5 or on each other; can be done any time, e.g. right after 0.1–0.3.

**Cleanup 1 — `TISE/make_and_run.sh`.** Currently runs `H-BoundStates` with zero args; `main.cpp:60` dereferences `argv[1]` unconditionally (confirmed: no `argc` check anywhere in `main.cpp`) → crash. Fix the script:
```bash
cmake -S . -B build && cmake --build build -j6
cd build && ./H-BoundStates '[{"domain": "(0, 100]", "function": "-1/x + 1/x^2"}]' ; cd -
```
Also add a defensive `argc` check in `main.cpp` immediately before line 60 — a trivial, self-contained, zero-risk addition (unlike cleanup 4, it isn't blocked by any missing infrastructure), so this fixes it rather than merely flagging it:
```cpp
if (argc < 2)
{
    std::cerr << "Usage: " << argv[0] << " '<JSON array of {\"domain\":...,\"function\":...} pieces>'\n";
    return EXIT_FAILURE;
}
```
(Manual/CLI check only — `main()` isn't reachable from GTest, so this doesn't move the coverage number.)

**Cleanup 2 — `TISE/README.md`.** Replace the "Changing the Angular Momentum" section and the hydrogenic-only intro prose with a "Specifying the Potential" section documenting the JSON-argv `{"domain": "<interval>", "function": "<muparser expression in x>"}` DSL (`main.cpp:58-59`'s own comment is the authoritative example), and note `constexpr int L` now only feeds the analytic-comparison formula, not the potential's shape. Add a short "Coverage" subsection documenting the 0.2 workflow, and a "Known Limitations" bullet for cleanup 4.

**Cleanup 3 — reconcile the expression-parser record.** `docs/SDD.md` §11.1 and `docs/planning/resources.md` both still describe the FunctionParser-vs-NFParam choice as "not yet chosen definitively," but `TISE/CMakeLists.txt:16-18` already commits to muparser (a third option, not even shortlisted). Since all 4 existing ADRs record "defer X" decisions and this is the opposite (an already-adopted, already-in-use choice), a docs update fits the repo's own precedent better than a new ADR:
- SDD §11.1: replace the "not yet chosen definitively" sentence with a statement that muparser is adopted, wired via `pkg_check_modules(MUPARSER ...)`, used by `tise::evaluateFunction`; note FunctionParser/NFParam were superseded.
- `docs/planning/resources.md`: annotate the muparser line (add `https://beltoforion.de/en/muparser/`, marked adopted) in the file's existing terse link-list style.

**Cleanup 4 — flag only, no behavior change.** Add a code comment above the unconditional `tevol::runTimeEvolution(...)` call (`main.cpp:121-133`) noting it always runs regardless of intent, pointing at `docs/planning/tise-task-breakdown.md` §4 item 4, plus one "Known Limitations" bullet in the README (added alongside cleanup 2). Do not add gating logic — that requires out-of-scope config plumbing (Phase 1).

---

## Section A1 — Boundary-condition asymptote classifier (REQ-F-030)

**New types/functions (`tise.hpp`/`tise.cpp`):**
```cpp
enum class DomainSide { Left, Right };
struct SpatialDomain { Real xMin; Real xMax; };
enum class AsymptoteCase { HardWall, AnalyticAsymptote, Irregular };       // Case 1 / 2 / 3
enum class AsymptoteSubType { NotApplicable, Flat, Coulomb };              // Case-2 sub-branch

struct AsymptoteClassification
{
    AsymptoteCase    asymptoteCase;
    AsymptoteSubType subType;
    Real             fittedAsymptoticValue; // V_inf; NaN if HardWall
    Real             powerLawExponent;      // fitted p in V ~ V_inf + C/x^p; NaN if HardWall or Flat
    bool             warningEmitted;
};

struct ConvergenceFit { bool isDivergent; bool isFlat; Real fittedLimit; Real powerLawExponent; };

// Shared numeric core — also reused by A4's singularity detection (Section A4).
ConvergenceFit classifySequenceConvergence(const std::vector<Real> &V, Real ratio);

AsymptoteClassification classifyAsymptote(const std::map<std::string, std::string> &potential,
                                           const SpatialDomain &domain,
                                           DomainSide side,
                                           std::ostream &warnOut = std::cerr);
```
`warnOut` defaults to `std::cerr`, mirroring `writeEigenstate`'s existing `std::ostream&` convention; tests inject an `std::ostringstream` instead of swapping `std::cerr`'s buffer.

**Algorithm.** Sample `V` at 16 points growing geometrically (ratio 4) away from the boundary: `x_k = reference ± scale·4^k`. Feed to `classifySequenceConvergence`:
1. **Flat check:** if successive differences `dV[k]` are all `≲ 1e-10 + 1e-9·max|V|` → `isFlat=true`.
2. **Divergence check:** compare `|dV|` magnitude in the back half of the sample window vs. the middle; if it hasn't shrunk by ≥2× → `isDivergent=true` (catches `x²`, `x`).
3. **Power-law fit:** otherwise, fit exponent `p` from `ln(|dV[k]|/|dV[k+1]|)/ln(ratio)` over a tail window (median of ~5 estimates), extrapolate the limit via a two-point combination using the fitted `p`.

`classifyAsymptote` maps the fit to a case: divergent → `HardWall`; flat → `AnalyticAsymptote`/`Flat`; `|p−1| ≤ 0.15` → `AnalyticAsymptote`/`Coulomb`; else → `Irregular` (emits the Case-3 warning to `warnOut`, message mentions "discontinuity" and "approximate").

Only the **classification** of Case 2's Coulomb sub-branch is in scope here — the matching formula for it is explicitly undocumented anywhere in the source material (flagged in the task breakdown as a follow-up for whoever picks up boundary-condition work next). We still report it as a distinct sub-type so a future consumer has the information available.

**Tests to add (`test_tise.cpp`):**
- `ClassifyAsymptoteTest.Case1HardWallForQuadraticGrowth` — `{{"[0, inf)", "x*x"}}` → `HardWall`.
- `ClassifyAsymptoteTest.Case1HardWallForLinearGrowth` — `{{"[0, inf)", "x"}}` → `HardWall`.
- `ClassifyAsymptoteTest.Case2FlatForStepPotential` — `{{"[0,5)","0"},{"[5, inf)","10"}}` → `AnalyticAsymptote`/`Flat`, `fittedAsymptoticValue≈10`.
- `ClassifyAsymptoteTest.Case2CoulombForInverseR` — `{{"(0, inf)", "-1/x"}}` → `AnalyticAsymptote`/`Coulomb`, `powerLawExponent≈1`.
- `ClassifyAsymptoteTest.Case3IrregularForPowerLawOneAndHalf` — `{{"(0, inf)", "1/x^1.5"}}` → `Irregular`, `powerLawExponent≈1.5`, `warningEmitted==true`.
- `ClassifyAsymptoteTest.Case3WarningTextMentionsDiscontinuity` — assert the injected stream contains "discontinuity"/"approximate".
- `ClassifyAsymptoteTest.NoWarningForCase1AndCase2` — injected stream stays empty for both.
- `ClassifyAsymptoteTest.LeftSideSymmetric` — `{{"(-inf, inf)", "x*x"}}`, `Left` side → `HardWall` (confirms sign-flipped sampling direction).
- `ClassifySequenceConvergenceTest.{DetectsFlatSequence, DetectsDivergentSequence, FitsCoulombPowerLaw, FitsArbitraryPowerLaw}` — synthetic arrays, no muparser involved, exercising the numeric core directly.

**Verification:** build+run `TISETests`; all "Done when" example potentials (x², flat/step, 1/r^1.5) classify correctly with the Case-3 warning firing; `gcovr` shows the new functions covered.

---

## Section A2 — Ionization-threshold / bound-continuum classification (REQ-F-020)

```cpp
struct BoundStateClassification { std::vector<bool> isBound; int nBound; };
BoundStateClassification classifyBoundStates(const EigenResult &result, Real threshold);
```
Standalone pass over `result.values` (strict `<`), does not mutate `EigenResult`. Linear scan; does not assume pre-sortedness even though `values` is ascending by construction.

**Tests to add:**
- `ClassifyBoundStatesTest.{NBoundMatchesCount, AllAboveThresholdGivesZeroBound, AllBelowThresholdGivesAllBound}` — synthetic `EigenResult{values={-3,-2,-1}, dim=3}`.
- `ClassifyBoundStatesTest.HydrogenicThresholdZeroSeparatesCorrectly` — reuse `SolveEigenTest`'s exact fixture setup (`nNodes=31, order=8, L=0, rMin=0, rMax=60`, hydrogenic expression), `threshold=0.0`. Assert ground state is bound, `nBound < dim`, and the "prefix of trues" property (ascending values ⇒ once a `false` appears scanning forward, all later entries are also `false`).
- `ClassifyBoundStatesTest.FiniteSquareWellBoundedCount` — new fixture, `nNodes=41, order=8, rMin=0, rMax=30`, potential `{{"[0,2)","-50.0"},{"[2,30]","0"}}`, `threshold=0.0`. **During implementation, actually run this once to read off the exact bound-state count and assert it exactly** (rough analytic estimate is ~3 states for this well depth/width; tighten once observed) — we can execute the solver ourselves, so there's no need to settle for a loose range check.

**Verification:** `gcovr` shows `classifyBoundStates` covered across all three branches (all-bound / all-above-threshold / mixed).

---

## Section A3 — Well-containment diagnostic

```cpp
struct ContainmentCheck { Real psiPrimeAtBoundary; bool notWellContained; };
ContainmentCheck checkWellContainment(const bspline::BSpline &bs,
                                       const std::vector<Real> &coeffs,
                                       Real xBoundary,
                                       Real tol = 1e-3);
```
Thin wrapper: `bs.eval(xBoundary, coeffs.data(), coeffs.size(), /*derivativeOrder=*/1)` compared against `tol`. `coeffs` is whatever `eigenstateCoefficients` produces (already zero-padded at both ends).

**Tolerance rationale:** DSBGV-normalized eigenvectors have O(1) L2 norm; a well-contained state's tail (and its derivative) decays exponentially, landing many orders of magnitude below `1e-3` well before the wall, while a wall-colliding state's derivative is generically O(0.1–1) — a comfortable margin either way. Test design is comparative (small vs. large box, same potential) rather than pinned to one exact magnitude.

**Tests to add** (reuse the hydrogenic radial-potential expression pattern, `order=8, L=0, rMin=0.1`):
- `WellContainmentTest.SmallBoxFlagsGroundState` — `nNodes=31, rMax=2.0` → `notWellContained == true`.
- `WellContainmentTest.LargeBoxDoesNotFlagGroundState` — same potential, `rMax=60.0` (matches the existing converged `SolveEigenTest` box) → `notWellContained == false`.
- `WellContainmentTest.DerivativeMatchesDirectBsEvalCall` — cross-check against a direct `bs.eval(..., 1)` call, `EXPECT_NEAR` at `1e-15`.

**Verification:** small-box case flags, large-box case doesn't; `gcovr` shows both branches of the flag hit.

---

## Section A4 — Strategic node placement (REQ-F-050)

**Delta-potential scoping decision:** a true Dirac delta cannot be represented in this DSL — muparser evaluates ordinary real expressions pointwise; there's no way to encode "infinite at a point, zero width." This is documented as a known limitation rather than attempted. A user wanting delta-like behavior represents it as a narrow, tall rectangular spike given as its own domain piece — the generic step-detection logic below handles both edges of that spike automatically, with no delta-specific code required.

**Scope split — A4a required, A4b stretch:**
- **A4a (required):** detect step/stitched-derivative/singular joins from the `map<string,string>` piece boundaries, build a grid with degenerate knots, wire it into `BSpline::init`. Fully satisfies the task's literal "Done when" (discontinuity → degenerate knots → improved eigenvalue accuracy at matched node count).
- **A4b (stretch, attempt only once A4a is solid):** generalize the hardcoded "drop B_1/B_N" convention (currently baked into `fillBandedMatrices`/`eigenstateCoefficients`) to an arbitrary caller-supplied drop-set, needed for the table's 4th row (singular-potential B-spline removal) to be more than detection-only. Implement as new defaulted-parameter overloads so existing call sites/tests are untouched.

**Types/functions:**
```cpp
enum class JoinType { Continuous, Step, StitchedKink, Singular };
struct DetectedJoin { Real x; JoinType type; };
struct StrategicKnot { Real x; int extraMultiplicity; };

std::vector<DetectedJoin> detectPotentialStructure(const std::map<std::string, std::string> &potential);
std::vector<StrategicKnot> strategicKnotsFromJoins(const std::vector<DetectedJoin> &joins, int order);
std::vector<Real> buildStrategicRadialGrid(int nNodes, Real rMin, Real rMax,
                                            const std::vector<StrategicKnot> &knots);
```

**Detection algorithm per numerically-adjacent piece pair.** Sort pieces by parsed numeric lower bound first — `std::map`'s native iteration order is lexicographic on the domain *string* (e.g. `"[10,...)"` sorts before `"[5,...)"`), which would silently misorder joins if used directly:
1. Singularity check first: sample each side approaching the shared boundary with a *shrinking* offset sequence, reusing `classifySequenceConvergence` from A1 — divergence on either side ⇒ `Singular`.
2. Else, one-sided finite differences (`h=1e-4`, Richardson-extrapolated) for value and slope at the boundary: value jump beyond tolerance ⇒ `Step`; else slope jump beyond a looser tolerance ⇒ `StitchedKink`; else `Continuous`.
3. Also probe each piece's own domain endpoints (not just inter-piece joins), so a single-piece potential like `{"(0, inf)": "-1/x + 1/x^2"}` still gets its `x→0` singularity flagged even though it has no "join" at all.

**Multiplicity formula, derived from the physics (double-checked independently, not just taken on faith):** for a B-spline basis of order `k` (degree `k-1`, matching this codebase's convention — `coeff(m,...)` uses degree `m` from 0 to `order-1`), a knot of *total* multiplicity `m` gives continuity `C^{k-1-m}`. To allow ψ's derivative of order `n` to be the first discontinuous one, we need continuity `C^{n-1}`, i.e. `m = k - n`. Since an ordinary grid point already has base multiplicity 1, the *extra* multiplicity to add is `m - 1 = k - n - 1`:
- **Step** (V itself jumps ⇒ from `ψ'' = 2m(V−E)ψ`, a V value-jump forces a ψ'' jump, `n=2`): `extraMultiplicity = order - 3`.
- **StitchedKink** (V continuous but V′ jumps ⇒ differentiating once more, `ψ''' = 2m(V'ψ + (V−E)ψ')` picks up V′'s jump, `n=3`): `extraMultiplicity = order - 4`.
(As a consistency check, not implemented: a true delta potential forces a ψ′ jump, `n=1`, giving `extraMultiplicity = order - 2` — one less than full multiplicity, i.e. continuity `C^0`, matching the standard delta-potential matching condition. This cross-check is what confirms the formula is right, independent of the source PDF below.)

*Note: the original PHY5606 course-note PDF this table is transcribed from is not present anywhere in this repository (confirmed by direct search) — the derivation above is our own from the Schrödinger equation itself, not a line-for-line check against that source. It's internally consistent (verified against the delta-potential case above), but flagging the provenance rather than silently presenting it as independently confirmed.*

**Grid wiring:** `buildStrategicRadialGrid` builds the same uniform base as `buildUniformRadialGrid`, then inserts `extraMultiplicity` repeated copies of each strategic knot at the correct sorted position. Confirmed directly from `BSpline.cpp`: `init()`'s monotonicity check only rejects strict decreases (repeats pass), and span computations already guard against zero-width spans — **no changes needed to `BSpline.cpp` itself**; repeated knots are already handled correctly by the existing math, just never previously exercised by a non-uniform caller.

**Tests to add:**
- `DetectPotentialStructureTest.StepAtBothBarrierEdges` — `{{"[0,5)","0"},{"[5,6]","10"},{"(6,10]","0"}}` (box+barrier) → `Step` at `x=5` and `x=6`.
- `DetectPotentialStructureTest.StitchedKinkAtCornerJoin` — `{{"[0,5)","x"},{"[5,10]","10-x"}}` → `StitchedKink` at `x=5`.
- `DetectPotentialStructureTest.SingularAtCoulombOrigin` — `{{"(0, 100]", "-1/x + 1/x^2"}}` (main.cpp's own default potential) → `Singular` at `x=0`.
- `DetectPotentialStructureTest.ContinuousJoinNotFlagged` — `{{"[0,5)","x*x"},{"[5,10]","x*x"}}` → no/`Continuous` join.
- `StrategicKnotsFromJoinsTest.MultiplicityFormulaPerOrder` — check `order-3`/`order-4` across a few orders (4, 6, 8).
- `BuildStrategicRadialGridTest.{InsertsRepeatedKnotAtRequestedLocation, RemainsNonDecreasing}`.
- `StrategicNodePlacementAccuracyTest.ImprovesOverUniformGridForBoxBarrier` — matches "Done when" literally: build a fine uniform reference grid (e.g. 201 nodes) as ground truth for the box+barrier potential, then compare a coarse uniform grid vs. a coarse strategic grid (both e.g. 21 nodes) against that reference — the strategic grid's eigenvalue error should be smaller.

**Verification:** the box+barrier example produces degenerate knots at both edges; eigenvalue accuracy improves at matched node count vs. pure-uniform; `gcovr` shows all four `JoinType` branches exercised.

---

## Section A5 — `E_acc` continuum-accuracy warning (REQ-F-040, warning half)

```cpp
Real computeEAcc(Real nodeSpacing, Real mass); // pi^2 / (2 * mass * nodeSpacing^2)
bool warnIfContinuumExceedsEAcc(Real eMax, Real eAcc, std::ostream &warnOut = std::cerr);
```
Standalone; no dependency on Engineer B's (not-yet-existing) energy-grid loop — wiring this at that call site is a future integration step, not part of this task, per the task breakdown's own framing. Note for whoever wires that call site later: once A4 lands and grids can be non-uniform, the physically-correct `nodeSpacing` to pass is the *minimum* inter-node gap, not a naive average — that's a call-site concern, not something baked into `computeEAcc` itself.

**Tests to add:**
- `ComputeEAccTest.MatchesClosedForm` — `nodeSpacing=0.5, mass=1.0` → `EXPECT_NEAR(computeEAcc(...), 2*M_PI*M_PI, 1e-12)`.
- `ComputeEAccTest.ScalesInverselyWithMassAndSpacingSquared` — doubling spacing quarters `E_acc`; doubling mass halves it.
- `WarnIfContinuumExceedsEAccTest.{FiresWhenEMaxExceeds, SilentWhenEMaxBelow, BoundaryEqualDoesNotFire}` (strict `>`, so exact equality stays silent, per the task's own "fires exactly when E_max > E_acc" wording).

**Verification:** formula matches closed form; warning fires exactly when `E_max > E_acc`; `gcovr` shows both functions covered.

---

## File-by-file summary

**Modified:** `TISE/tise.hpp`, `TISE/tise.cpp` (all new declarations/implementations, grouped by task with comment banners), `TISE/tests/test_tise.cpp` (all new `TEST`/`TEST_F` cases), `TISE/CMakeLists.txt` (`ENABLE_COVERAGE`), `TISE/make_and_run.sh`, `TISE/main.cpp` (argc check + cleanup-4 comment), `TISE/README.md`, `docs/SDD.md` (§11.1 only), `docs/planning/resources.md`.

**Created:** none in `TISE/` — no new source/test files; everything extends the existing shared files, per the task breakdown's own merge-conflict-avoidance guidance (§5). (This plan document and its per-task companions under `docs/planning/` are the only new files.)

## Verification (run after each section)
```bash
cd TISE
cmake -S . -B build -DBUILD_TESTING=ON && cmake --build build -j && ctest --test-dir build --output-on-failure   # functional
# then the coverage workflow from Section 0.2 — confirm tise.cpp stays >= 80%
```
Manual sanity check for cleanup 1 specifically: `bash make_and_run.sh` runs without crashing.

## Notes / assumptions carried forward into implementation
- No commits will be made automatically; commit timing is decided separately, after reviewing each section.
- A4's derivative-discontinuity-order derivation is our own working interpretation (the source PDF isn't present in the repo) — flagged above, not silently assumed-verified; it is internally self-consistent (cross-checked against the delta-potential case).
- A2's finite-square-well exact bound-state count will be pinned empirically during implementation (we can actually run the solver), not left as a guess or loose range.
- The 80% coverage target is scoped to `tise.cpp` (where all new code lands), not retroactively applied to `BSpline.cpp`/`time_evolution.cpp`'s pre-existing coverage — that's unrelated to this task set.
