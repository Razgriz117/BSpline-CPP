# TISE Engineer-A — Task A2 Implementation Plan (active slice)

> This plan is scoped to **task A2 only** (ionization-threshold determination + bound/continuum classification, REQ-F-020). Task A1 (boundary-condition asymptote classifier) is complete and committed (`15db181`). The full five-task master plan is archived at [`engineer-a-plan.md`](engineer-a-plan.md); A1's own plan is at [`engineer-a-plan-A1.md`](engineer-a-plan-A1.md).

## Context

`docs/planning/tise-task-breakdown.md` / REQ-F-020: the solver always diagonalizes the full box and reports **all** sub-threshold eigenvalues as bound states — the count is an output, never a user input. A2 implements the `CLASSIFY` node itself: given an already-solved `EigenResult` and an ionization threshold (a plain `Real`, supplied by the caller — not computed here; that's A1's `classifyAsymptote` job, or 0.0 by convention for a Coulomb-type problem), label each eigenvalue bound or above-threshold.

**Environment/repo state (re-verified, no drift since A1 landed):** the only commit since A1's prerequisite work (`1a52a8a "Interface (#8)"`) is A1's own commit; working tree clean; dependencies, `ENABLE_COVERAGE`, and the full test suite are all already in place and green. No Section-0-style setup needed this round.

## Section A2 — Ionization-threshold / bound-continuum classification (REQ-F-020)

### Types/function (`tise.hpp` — insert between `solveGeneralizedEigenproblem`'s declaration and `analyticHydrogenEnergy`; `tise.cpp` — insert between `solveGeneralizedEigenproblem`'s definition and `analyticHydrogenEnergy`)

```cpp
struct BoundStateClassification
{
    std::vector<bool> isBound; // isBound[i] true iff values[i] < threshold
    int nBound;                // count of true entries in isBound
};

BoundStateClassification classifyBoundStates(const EigenResult &result, Real threshold);
```

Implementation iterates `result.values.size()` (not `result.dim`) — `dim` is a redundant/derived field that's easy to leave stale in a hand-built `EigenResult` (several tests below construct one directly); scanning `values` directly is correct regardless. Plain linear scan, strict `<` (state exactly at threshold is not bound), no assumption that `values` is pre-sorted.

```cpp
BoundStateClassification classifyBoundStates(const EigenResult &result, Real threshold)
{
    BoundStateClassification out;
    out.isBound.resize(result.values.size());
    out.nBound = 0;
    for (std::size_t i = 0; i < result.values.size(); ++i)
    {
        bool bound = result.values[i] < threshold;
        out.isBound[i] = bound;
        if (bound) ++out.nBound;
    }
    return out;
}
```

### Finite-square-well test — physics independently verified (not a "run it and see" placeholder)

Odd-parity bound states of the symmetric well ≡ bound states of `-ψ''/2 + V(x)ψ = Eψ` on `[0,∞)`, `ψ(0)=0`, `V=-V0` for `0≤x<w`, `V=0` for `x≥w` (matches this codebase's Dirichlet-at-origin convention). Verified three independent ways via Python during planning:

1. **Transcendental root-finding** (`k·cot(k·w) = -κ`, `k²+κ²=2V0`) for `V0=2.0, w=5.0` (`k0=2.0`): 3 roots — `k=0.5705→E=-1.8373`, `k=1.1358→E=-1.3549`, `k=1.6846→E=-0.5810`.
2. **Independent cross-check**: dense finite-difference diagonalization of the *full* symmetric well (`x∈[-40,40]`, dx=0.01), classifying eigenvectors by parity — 7 total bound states, 3 odd-parity, energies matching method 1 to ~1e-3.
3. **Box-truncation stability check**: same FD problem on `[0,rMax]` with Dirichlet at both ends, scanning `rMax` from 10 to 30 — count and energies stable to ~1e-5 throughout. Weakest state's decay length is `1/κ=0.928`; `rMax=30` gives ~27 decay lengths of margin beyond the well.

**Chosen parameters:** `V0=2.0`, `w=5.0`, `rMin=0.0`, `rMax=30.0`, `order=8`. **Exact expected count: 3 bound states**, threshold `0.0`.

```cpp
std::map<std::string, std::string> potential = {
    {"[0, 5)",  "-2.0"},
    {"[5, 30]", "0.0"}
};
```

**Implementation result (2026-07-31):** the bound-state **count** was exactly 3 at `nNodes=31` as predicted — `classifyBoundStates` itself is correct. The *energies* converged more slowly than 1e-3 at `nNodes=31` (errors up to ~1e-2), because the potential's step discontinuity at `x=5` isn't yet resolved with degenerate knots (that's task A4's strategic node placement, not yet implemented). A direct convergence scan (`nNodes` = 31/51/81/121/161) confirmed the eigenvalues do converge toward the predicted analytic roots, just slowly; the `ClassifyBoundStatesSquareWellTest` fixture (shared by both the count test and the energy test) was bumped to `nNodes=121` (energy errors ~2e-4, comfortably under `1e-3`). This is exactly the anticipated fallback ("bump `nNodes`, don't touch `V0`/`w`/`rMax`") — just needed for energy precision rather than count correctness, and it costs ~2.7s per test (up from negligible) which is an acceptable one-time trade for a physics-verified regression test.

### Tests to add (`test_tise.cpp`, new banner section between the end of `SolveEigenTest`'s tests and the `writeEigenstate` banner)

Synthetic (plain `TEST`, hand-built `EigenResult`):
- `ClassifyBoundStatesTest.AllBoundWhenAllBelowThreshold` — `values={-5,-3,-1}` → `nBound==3`, all true.
- `ClassifyBoundStatesTest.AllAboveThresholdWhenNoneBelow` — `values={0.5,1.0,2.0}` → `nBound==0`, all false.
- `ClassifyBoundStatesTest.MixedAboveAndBelowThreshold` — `values={-2,-1,0.5,1.5}` → `nBound==2`, pattern `{true,true,false,false}`.
- `ClassifyBoundStatesTest.ValueExactlyAtThresholdIsNotBound` — `values={-1,0,1}`, threshold `0.0` → index 1 is `false` (strict `<`), `nBound==1`.

(Manual per-index `EXPECT_TRUE`/`EXPECT_FALSE`, not gmock's `ElementsAre` — this codebase links only `GTest::gtest`/`gtest_main`, no gmock.)

Hydrogenic structural check — new fixture, `SetUp` copied verbatim from `SolveEigenTest::SetUp`:
- `ClassifyBoundStatesHydrogenTest.GroundStateBoundNBoundLessThanDimAndPrefixStructureHolds` — threshold `0.0`: ground state bound, `nBound < result.dim`, "prefix of trues" property holds.

Square-well fixture (own independent `SetUp`, same grid/basis pattern, potential as above):
- `ClassifyBoundStatesSquareWellTest.ExactlyThreeBoundStatesBelowZero` — `nBound==3`, first 3 true, 4th false.
- `ClassifyBoundStatesSquareWellTest.BoundEnergiesMatchAnalyticOddParityRoots` — `EXPECT_NEAR` each of the 3 energies against `{-1.83728, -1.35493, -0.58099}` at `1e-3`.

### Verification
Build + run `TISETests`; square-well test confirms the physics derivation transfers to the real B-spline solver; rerun `gcovr` (`engineer-a-plan-A1.md` Section 0.2) and confirm `tise.cpp` stays ≥80% (currently ~94%).

## File-by-file summary
**Modified:** `TISE/tise.hpp`, `TISE/tise.cpp` (new struct + function), `TISE/tests/test_tise.cpp` (new test section).

## Notes / assumptions
- No commits made automatically.
- A3–A5 and cleanup tasks 1–4 remain deferred, tracked in `engineer-a-plan.md`.
- The finite-square-well parameters are a fresh, independently-verified derivation (three cross-checking methods).
