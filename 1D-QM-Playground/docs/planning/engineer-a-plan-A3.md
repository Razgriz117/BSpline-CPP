# TISE Engineer-A — Task A3 Implementation Plan (active slice)

> This plan is scoped to **task A3 only** (well-containment diagnostic). Tasks A1 (boundary-condition asymptote classifier) and A2 (bound/continuum classification) are complete and committed (`15db181`, `725dacc`). The full five-task master plan is archived at [`engineer-a-plan.md`](engineer-a-plan.md); A1's and A2's own plans are at [`engineer-a-plan-A1.md`](engineer-a-plan-A1.md) / [`engineer-a-plan-A2.md`](engineer-a-plan-A2.md). This plan corrects two details in the master plan's A3 section against verified current code (below) and elaborates it to the same rigor as A1/A2.

## Context

`docs/planning/tise-task-breakdown.md` (§2, "A3. Well-containment diagnostic") and `docs/SDD.md` (§5.2.3 Figure 6 `BOUND` node, §6.4 data-validation rules) specify a diagnostic: for a bound state, `ψ'(x_max) ≠ 0` (beyond tolerance) flags that the state is "colliding" with the box wall, so its energy/wavefunction may be inaccurate. This has **no REQ-F-0XX ID** — confirmed by checking the full requirements table in `docs/SDD.md:199-204` and the task-breakdown's own heading (line 74), which conspicuously omits the tag every other Engineer-A task has. It's prose-specified only.

**Environment/repo state (re-verified, no drift since A2 landed):** `HEAD` is `725dacc`, working tree clean, nothing landed since A2. No Section-0-style setup needed this round.

## Ground truth verified directly against current code (not just the master plan's prose)

- `bspline::BSpline::eval(Real x, const Real *fc, int fcSize, int derivativeOrder=0, bool skipFirst=false, int Bsmin=-1, int Bsmax=-1) const` — `TISE/BSpline.hpp:87-93`, impl `BSpline.cpp:483-542`. Exactly matches the master plan's assumed call shape. Out-of-range `x` returns `0.0` silently (no throw); exact equality at either wall is valid/in-range (traced `whichInterval`, `BSpline.cpp:384-431`).
- `eigenstateCoefficients` (`tise.hpp:160-166`/`tise.cpp:384-394`) confirmed zero-padded at both ends (`coeffs[0]=0`, `coeffs[nBSplines-1]=0`), size `nBSplines` — exactly what `bs.eval`'s `fcSize` expects.
- `EigenResult` (`tise.hpp:15-21`) has `.vectors` (column-major eigenvectors, size `dim*dim`) — the correct field to pass as `eigenstateCoefficients`'s `evec` argument.
- **Correction to the master plan's prose:** `engineer-a-plan.md:174` says A3's tests should use "`order=8, L=0, rMin=0.1`" — this doesn't match any real fixture. The two fixtures that actually run the full hydrogenic pipeline at `order=8, L=0` — `SolveEigenTest` (`test_tise.cpp:621-651`) and A2's `ClassifyBoundStatesHydrogenTest` (`test_tise.cpp:736-763`) — both use **`rMin=0.0`**, confirmed by direct read. (`rMin=0.1` only appears on two unrelated `order=4` fixtures that never reach `solveGeneralizedEigenproblem`.) The master plan's own A2 section (`engineer-a-plan.md:154`) correctly cites `SolveEigenTest` as `rMin=0` — its A3 section contradicts this three lines later. **A3's tests use `rMin=0.0`**, matching `SolveEigenTest` exactly (verified bit-exact grid endpoints for both `rMax=2.0` and `rMax=60.0`, so no floating-point boundary risk either way).
- Insertion points, precisely located: `tise.hpp` between `classifyBoundStates`'s declaration (ends line 152) and `analyticHydrogenEnergy` (line 154); `tise.cpp` between `classifyBoundStates`'s impl (ends line 371) and `analyticHydrogenEnergy` (line 373); `test_tise.cpp` between A2's last test (ends line 835) and the `writeEigenstate` banner (line 837). No new `#include`s needed (`<cmath>` already included in `tise.cpp:5`).
- Banner-comment convention has been inconsistent so far (A1 bannered `.cpp` only, A2 bannered `.hpp` only — confirmed by grep, not a deliberate rule). **A3 banners both files**, establishing a cleaner convention going forward without retroactively touching A1/A2's already-committed code.

## Design

### `tise.hpp` (insert after line 152)
```cpp
// === A3: well-containment diagnostic (SDD Sec. 5.2.3, Sec. 6.4) ===
// Result of checkWellContainment: the raw signed derivative plus the
// pass/fail flag derived from it.
struct ContainmentCheck
{
    Real psiPrimeAtBoundary; // psi'(xBoundary), as returned by BSpline::eval
    bool notWellContained;   // true iff |psiPrimeAtBoundary| > tol
};

// Well-containment diagnostic (SDD Sec. 5.2.3 Figure 6 "BOUND" node; Sec. 6.4
// data-validation rule): a bound state confined within the box should decay
// to numerically-zero slope at the outer wall; a nonzero psi'(xBoundary)
// means the state is "colliding" with the wall and its energy/wavefunction
// may be inaccurate. `coeffs` is whatever eigenstateCoefficients produces
// (already zero-padded at both ends) for a single eigenstate; `xBoundary` is
// typically rMax (or rMin for a left-side check). Strict '>': a derivative
// magnitude exactly equal to tol is the marginal case and is not flagged,
// mirroring classifyBoundStates' strict-inequality boundary convention.
ContainmentCheck checkWellContainment(const bspline::BSpline &bs,
                                       const std::vector<Real> &coeffs,
                                       Real xBoundary,
                                       Real tol = 1e-3);
```

### `tise.cpp` (insert after line 371)
```cpp
// === A3: well-containment diagnostic (SDD Sec. 5.2.3, Sec. 6.4) ===
ContainmentCheck checkWellContainment(const bspline::BSpline &bs,
                                       const std::vector<Real> &coeffs,
                                       Real xBoundary,
                                       Real tol)
{
    const int n = static_cast<int>(coeffs.size());
    ContainmentCheck out;
    out.psiPrimeAtBoundary = bs.eval(xBoundary, coeffs.data(), n, 1);
    out.notWellContained = std::abs(out.psiPrimeAtBoundary) > tol;
    return out;
}
```
`static_cast<int>(coeffs.size())` matches the existing call-site style precedent in `writeEigenstate` (`tise.cpp:404,409`) rather than passing `.size()` inline.

### Tolerance/physics reasoning (independently verified, not just trusted from the master plan)
With `L=0` the test expression reduces to the standard hydrogen radial equation (`u(r)=rR(r)` substitution), ground state `E=-0.5`, `u(r) ∝ r·e^{-r}` (decay length 1 in these units) — confirmed `fillBandedMatrices` uses the flat measure (`fUni` weight, `tise.cpp:282,304-305`), i.e. `coeffs` really does reconstruct `u(r)`, not `R(r)`. Using the normalized ground state `u(r)=2r·e^{-r}`:
- `rMax=2.0` (only ~2 decay lengths of margin) → predicted `|ψ'(2.0)| ≈ 0.271`, ~270× above `tol=1e-3` → flags, not borderline.
- `rMax=60.0` (matches `SolveEigenTest`'s converged box) → predicted `|ψ'(60.0)| ≈ 1.0×10⁻²⁴`, ~21 orders of magnitude below `tol` → doesn't flag, not borderline.

Both margins are enormous — this test is not flaky. (Mirrors A2's square-well precedent of citing decay-length margin, e.g. "~27 decay lengths.") **Do not add an energy-accuracy assertion to the small-box fixture** — the ground-state energy visibly deviating from -0.5 at `rMax=2.0` is the point of the test, not a defect.

### Scope: standalone, unwired (matches A1/A2 precedent)
`checkWellContainment` lands as a tested standalone function only — **not** wired into `main.cpp`/`solveTISE`/`writeEigenstate`'s output, even though the SDD implies eventual integration ("should be computed and reported alongside each state"). Neither A1 nor A2 touched those call sites either; wiring diagnostics into actual CLI/report output is a later integration concern (same framing A5's plan already uses for its own call-site wiring).

### Tests added (`test_tise.cpp`, inserted after line 835, before the `writeEigenstate` banner)

GTest mechanics rule out one literal `WellContainmentTest` fixture reused across a small-box and large-box case the way `SolveEigenTest` reuses one `SetUp()` — `rMax` differs per case and `SetUp()` takes no parameters. Following this file's own established habit of separate fixture classes per scenario (e.g. `ClassifyBoundStatesHydrogenTest` vs. `ClassifyBoundStatesSquareWellTest`), and `WriteEigenstateTest`'s precedent of storing `bs`/`coeffs` as members (needed here since the test bodies call `checkWellContainment(bs, coeffs, ...)` directly, unlike `SolveEigenTest` which only needs `result`):

- `WellContainmentSmallBoxTest.FlagsGroundState` — `nNodes=31, order=8, L=0, rMin=0.0, rMax=2.0`, hydrogenic expression → `notWellContained == true`.
- `WellContainmentLargeBoxTest.DoesNotFlagGroundState` — same pattern, `rMax=60.0` (matches `SolveEigenTest`'s converged box) → `notWellContained == false`.
- `WellContainmentLargeBoxTest.DerivativeMatchesDirectBsEvalCall` — cross-check against a direct `bs.eval(..., 1)` call, `EXPECT_NEAR` at `1e-15`.
- `WellContainmentSyntheticTest.DerivativeExactlyAtToleranceIsNotFlagged` — cheap non-physics fixture (no LAPACK/muparser), coefficients ramp from 1 (interior) to 0 (both ends). Pins the exact `>` boundary: passing `tol == |raw|` must NOT flag (marginal case).
- `WellContainmentSyntheticTest.DerivativeJustAboveToleranceIsFlagged` — same fixture, `tol` set just under `|raw|` → must flag.

(Mirrors A2's pairing of synthetic threshold tests alongside physics-fixture tests — A2's `ValueExactlyAtThresholdIsNotBound` is the direct precedent for the last two.)

## Scope / non-goals
- No `CMakeLists.txt`, `main.cpp`, `README.md`, or `docs/SDD.md` changes.
- `docs/planning/engineer-a-plan.md` (the master doc) is left untouched, matching A2's precedent — the `rMin`/banner corrections above live only in this doc.
- No commits made automatically; commit timing is decided after review, matching A1/A2.

## File-by-file summary
**Modified:** `TISE/tise.hpp`, `TISE/tise.cpp` (new struct + function), `TISE/tests/test_tise.cpp` (3 new fixture classes, 5 new `TEST_F` cases).

## Verification
```bash
cd TISE
cmake -S . -B build -DBUILD_TESTING=ON && cmake --build build -j && ctest --test-dir build --output-on-failure
```
```bash
cmake -S . -B build-coverage -DBUILD_TESTING=ON -DENABLE_COVERAGE=ON
cmake --build build-coverage -j
find build-coverage -name '*.gcda' -delete
ctest --test-dir build-coverage --output-on-failure
mkdir -p build-coverage/coverage
gcovr --root . --filter '.*tise\.cpp' --exclude '.*/tests/.*' --exclude '.*main\.cpp' \
      --print-summary --html --html-details -o build-coverage/coverage/index.html build-coverage
```

## Notes / assumptions
- No commits made automatically.
- A3 has no `REQ-F-0XX` ID (verified against the full SDD requirements table) — the code banner cites SDD section numbers instead of inventing one.
- DSBGV's eigenvector sign is arbitrary (no LAPACK guarantee); already handled correctly — `notWellContained` uses `std::abs`, and no test asserts a specific sign on `psiPrimeAtBoundary`.
- A4–A5 and cleanup tasks 1-4 remain deferred, tracked in `engineer-a-plan.md`.
- Wiring this diagnostic into the actual solve/report pipeline (SDD's "reported alongside each state") remains a deferred integration concern, not part of this task.

## Implementation result (2026-08-08)

Implemented via TDD: added the `.hpp` declaration only, added all 5 tests, built and confirmed a **link-time RED** (`undefined reference to tise::checkWellContainment(...)`, all 5 new tests, no other errors — confirmed failing for the right reason, not a typo), then added the `.cpp` implementation and rebuilt to **GREEN**.

All 5 new tests passed on the first try, including the one flagged as an open risk in this doc's design section — the synthetic fixture's `ψ'(xBoundary) ≠ 0` assumption (ramp-to-zero coefficients on an `nNodes=11, order=4` grid) held with no fallback needed. Full suite: `ctest --test-dir build` → 4/4 test binaries pass (`BSplineTests`, `UtilsTests`, `TISETests`, `TimeEvolutionTests`), no regressions.

Coverage (`gcovr`, filtered to `tise.cpp`): **98.6% line coverage** (965/979 lines), up from A2's ~94% baseline — comfortably clears the ≥80% (REQ-NF-010) bar. `checkWellContainment`'s single branch point is fully exercised (small-box/large-box pair hits both `notWellContained` arms; synthetic pair pins the exact `>` boundary).

No changes needed to the design as planned — the physics-margin estimates (`|ψ'(2.0)|≈0.271` flags, `|ψ'(60.0)|≈10⁻²⁴` doesn't) and the `rMin=0.0` correction both held exactly as predicted.
