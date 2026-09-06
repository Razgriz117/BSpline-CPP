# TISE Engineer-A — Task A1 Implementation Plan (active slice)

> This plan is scoped to **task A1 only** (boundary-condition asymptote classifier, REQ-F-030). The full five-task master plan (A1–A5 + cleanup) is archived at [`engineer-a-plan.md`](engineer-a-plan.md) in this directory for later rounds.

## Context

`docs/SDD.md` / `docs/planning/tise-task-breakdown.md` assign Engineer A the `CLASSIFY`→`BOUND` branch of the TISE solver's Figure 6/7 flow. Task A1 is the first and most foundational piece: for each unbounded domain side, classify the potential's asymptote into Case 1 (hard wall) / Case 2 (flat or Coulomb, analytic) / Case 3 (irregular — approximate treatment, warn). REQ-NF-010 requires ≥80% coverage maintained on every change; no coverage tooling exists in the repo yet, so a minimal local `gcovr`-based setup is a prerequisite (Section 0 below), not scope creep.

**Case-3 treatment updated per `docs/planning/boundary-condition-case-3-smoothing.md`:** a smooth `sin²`-taper window,
```
Ṽ(x) = W(x)·V(x),   W(x) = 1 for x < R−Δ, sin²(π/2·(x−R)/Δ) for R−Δ≤x≤R, 0 for x>R
```
replaces the flat-truncation-with-discontinuity originally planned for Case 3. This removes the hard discontinuity the original plan explicitly warned about, avoids relying on a fitted asymptotic value for a tail that's irregular by definition, and avoids needing special step-knot treatment at the box edge for Case-3 potentials in a later node-placement task. Trade-off: Δ becomes a new parameter needing a sensible default, and the warning wording changes (no longer claims a "discontinuity" that no longer exists) but is still emitted, since this remains an approximation.

**Environment facts confirmed directly (not assumptions):**
- `BSpline-CPP/` (project root, one level above `1D-QM-Playground/`) is already a git repo, branch `TISE-Generalization`, clean tree, tracks `origin`. No `git init` needed.
- The actual dev environment is missing `GTest`, `muparser`, `nlohmann-json`, and `gcovr`/`lcov` — nothing in `TISE/` can build or test until these are installed (Section 0.1).

---

## Section 0 — Prerequisites (shared, do first)

### 0.1 Install missing dependencies, confirm a clean baseline
```bash
sudo apt-get install libgtest-dev libmuparser-dev nlohmann-json3-dev liblapack-dev libblas-dev libeigen3-dev libyaml-cpp-dev
pip install --user gcovr   # or: sudo apt-get install gcovr lcov
cd TISE && cmake -S . -B build -DBUILD_TESTING=ON && cmake --build build -j && ctest --test-dir build --output-on-failure
```
(`libyaml-cpp-dev` wasn't anticipated when this plan was written — `CMakeLists.txt` picked up a hard `find_package(yaml-cpp REQUIRED)` from a separate, already-landed "Interface (#8)" commit implementing the Controller↔TISE stub CLI, unrelated to Engineer A's work but now a build prerequisite regardless.)

### 0.2 Coverage tooling (`TISE/CMakeLists.txt`)
Insert before line 9 (`option(BUILD_TESTING ...)`), so it lands before `add_subdirectory(tests)` at line 23:
```cmake
option(ENABLE_COVERAGE "Build with --coverage instrumentation (gcov/gcovr)" OFF)
if(ENABLE_COVERAGE)
    add_compile_options(--coverage -O0 -g)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} --coverage")
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} --coverage")
endif()
```
Measurement workflow (rerun after A1 lands):
```bash
cmake -S . -B build-coverage -DBUILD_TESTING=ON -DENABLE_COVERAGE=ON
cmake --build build-coverage -j
find build-coverage -name '*.gcda' -delete
ctest --test-dir build-coverage --output-on-failure
mkdir -p build-coverage/coverage
gcovr --root . --filter '.*tise\.cpp' --exclude '.*/tests/.*' --exclude '.*main\.cpp' \
      --print-summary --html --html-details -o build-coverage/coverage/index.html build-coverage
```
`build-coverage/` matches `.gitignore`'s `build*` pattern. Target is ≥80% on `tise.cpp`'s own reported percentage.

### 0.3 Baseline backfill
`inInterval`/`evaluateFunction` (`tise.cpp:48-111`) have zero direct tests and aren't declared in `tise.hpp`. Promote both to header declarations (pure addition) and add:
- `InIntervalTest.{ClosedBothEnds, HalfOpenVariants, InfinityBounds, ThrowsOnMalformedString}`
- `EvaluateFunctionTest.{PicksCorrectPiece, ThrowsWhenXUncovered}`

This establishes an honest baseline before A1's new code lands.

*(Cleanup tasks 1–4 are deferred — tracked in `engineer-a-plan.md`, not part of this round.)*

---

## Section A1 — Boundary-condition asymptote classifier (REQ-F-030)

### Core types/functions (`tise.hpp`/`tise.cpp`)
```cpp
enum class DomainSide { Left, Right };
struct SpatialDomain { Real xMin; Real xMax; };
enum class AsymptoteCase { HardWall, AnalyticAsymptote, Irregular };       // Case 1 / 2 / 3
enum class AsymptoteSubType { NotApplicable, Flat, Coulomb };              // Case-2 sub-branch

struct AsymptoteClassification
{
    AsymptoteCase    asymptoteCase;
    AsymptoteSubType subType;
    Real             fittedAsymptoticValue;   // V_inf; NaN if HardWall
    Real             powerLawExponent;        // fitted p in V ~ V_inf + C/x^p; NaN if HardWall or Flat
    Real             recommendedTransitionWidth; // Delta for the Case-3 window; NaN unless Irregular
    bool             warningEmitted;
};

struct ConvergenceFit { bool isDivergent; bool isFlat; Real fittedLimit; Real powerLawExponent; };

ConvergenceFit classifySequenceConvergence(const std::vector<Real> &V, Real ratio);

AsymptoteClassification classifyAsymptote(const std::map<std::string, std::string> &potential,
                                           const SpatialDomain &domain,
                                           DomainSide side,
                                           std::ostream &warnOut = std::cerr);
```
`warnOut` defaults to `std::cerr` (mirrors `writeEigenstate`'s existing `std::ostream&` convention); tests inject an `std::ostringstream`.

### Classification algorithm
Sample `V` at 16 points growing geometrically (ratio 4) away from the boundary: `x_k = reference ± scale·4^k`. Feed to `classifySequenceConvergence`:
1. **Flat check:** successive differences `dV[k]` all `≲ 1e-10 + 1e-9·max|V|` → `isFlat=true`.
2. **Divergence check:** `|dV|` in the back half of the window vs. the middle; if it hasn't shrunk ≥2× → `isDivergent=true` (catches `x²`, `x`).
3. **Power-law fit:** else fit exponent `p` from `ln(|dV[k]|/|dV[k+1]|)/ln(ratio)` over a tail window (median of ~5 estimates), extrapolate the limit via a two-point combination.

Mapping: divergent → `HardWall`; flat → `AnalyticAsymptote`/`Flat`; `|p−1| ≤ 0.15` → `AnalyticAsymptote`/`Coulomb`; else → `Irregular`.

### Case 3 handling — updated per the smoothing doc
```cpp
// W(x): sin^2 raised-cosine taper. 1 for x <= R-delta, 0 for x >= R, smooth
// (C^1: value and slope both match) in between. Side-aware: for a Left
// boundary, mirror the taper direction (x >= R+delta -> 1, x <= R -> 0).
Real case3WindowFunction(Real x, Real R, Real delta, DomainSide side);

// V~(x) = W(x) * V(x), per docs/planning/boundary-condition-case-3-smoothing.md.
// Evaluates the caller-supplied potential via evaluateFunction and applies the taper.
Real evaluateWindowedPotential(const std::map<std::string, std::string> &potential,
                                Real x, Real R, Real delta, DomainSide side);
```
Default `Δ`: `0.1 * (domain.xMax - domain.xMin)` (10% of the box size) unless the caller overrides it.

When `classifyAsymptote` detects `Irregular`, it populates `recommendedTransitionWidth` and emits an **updated** warning (no longer claiming a discontinuity, since the taper removes it):
> "Warning: potential asymptote on the `<side>` side is irregular (fitted power-law exponent p=...); the potential will be smoothly tapered to zero over a transition width Δ=... approaching the box boundary (`docs/planning/boundary-condition-case-3-smoothing.md`) rather than truncated with a sharp discontinuity. This remains an approximation — the true asymptotic tail is not analytically known — and continuum normalization will be approximate."

This is a deliberate, reasoned departure from the task-breakdown's literal "warn... this introduces a discontinuity" wording: with the window in place there is no discontinuity to warn about, so asserting one would be inaccurate. The warning's *purpose* (flag an approximation, non-fatal, `stderr`, per SDD §8's physics-warning taxonomy) is preserved; only the specific claim is corrected to match the improved numerics.

Only the *classification* of Case 2's Coulomb sub-branch is in scope — its matching formula is undocumented anywhere in the source material (flagged as a follow-up, unchanged from the original plan).

### Tests to add (`test_tise.cpp`)
- `ClassifyAsymptoteTest.Case1HardWallForQuadraticGrowth` — `{{"[0, inf)", "x*x"}}` → `HardWall`.
- `ClassifyAsymptoteTest.Case1HardWallForLinearGrowth` — `{{"[0, inf)", "x"}}` → `HardWall`.
- `ClassifyAsymptoteTest.Case2FlatForStepPotential` — `{{"[0,5)","0"},{"[5, inf)","10"}}` → `AnalyticAsymptote`/`Flat`, `fittedAsymptoticValue≈10`.
- `ClassifyAsymptoteTest.Case2CoulombForInverseR` — `{{"(0, inf)", "-1/x"}}` → `AnalyticAsymptote`/`Coulomb`, `powerLawExponent≈1`.
- `ClassifyAsymptoteTest.Case3IrregularForPowerLawOneAndHalf` — `{{"(0, inf)", "1/x^1.5"}}` → `Irregular`, `powerLawExponent≈1.5`, `warningEmitted==true`, `recommendedTransitionWidth` populated and positive.
- `ClassifyAsymptoteTest.Case3WarningTextMentionsTaperingNotDiscontinuity` — assert the injected stream mentions "taper"/"transition" and "approximate", and does **not** claim a discontinuity.
- `ClassifyAsymptoteTest.NoWarningForCase1AndCase2` — injected stream stays empty for both.
- `ClassifyAsymptoteTest.LeftSideSymmetric` — `{{"(-inf, inf)", "x*x"}}`, `Left` side → `HardWall`.
- `ClassifySequenceConvergenceTest.{DetectsFlatSequence, DetectsDivergentSequence, FitsCoulombPowerLaw, FitsArbitraryPowerLaw}` — synthetic arrays, no muparser.
- `Case3WindowFunctionTest.{EqualsOneWellInsideBoundary, EqualsZeroAtAndBeyondBoundary, IsContinuousAtBothTransitionEdges, DerivativeVanishesAtBothTransitionEdges, MonotonicWithinTransition}`.
- `EvaluateWindowedPotentialTest.{MatchesRawPotentialWellInsideBoundary, VanishesAtAndBeyondBoundary, MatchesDirectWindowFunctionTimesRawPotential}`.

### Verification
Build+run `TISETests`; all four example potentials (x², flat/step, Coulomb 1/x, 1/r^1.5) classify correctly; Case-3 warning fires with the updated wording; window function is confirmed C¹-continuous at both transition edges numerically (finite-difference derivative check); `gcovr` shows `tise.cpp` ≥80% with the new functions covered.

## File-by-file summary (this round only)
**Modified:** `TISE/tise.hpp`, `TISE/tise.cpp` (new declarations/implementations for A1 + the Section-0.3 backfill), `TISE/tests/test_tise.cpp` (all new test cases above), `TISE/CMakeLists.txt` (`ENABLE_COVERAGE`).

## Notes / assumptions
- No commits made automatically.
- A2–A5 and cleanup tasks 1–4 remain deferred, tracked in `engineer-a-plan.md`.
- Δ's default (10% of box size) is a reasonable starting point, not derived from the smoothing doc itself (which only gives the window's functional form, not a default width) — flagging this as our own engineering choice, adjustable once a consumer (e.g. continuum matching) has concrete accuracy requirements.
