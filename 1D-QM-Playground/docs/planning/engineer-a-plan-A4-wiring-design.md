# TISE Engineer-A — A4/A4b Wiring: Design (promoted from the gap inventory)

> Promotes `engineer-a-plan-A4-wiring.md` (2026-08-21 gap inventory, "findings only, not yet a design... do not implement straight from this list") into a full design pass, per that doc's own "Next step" instruction. A1 (`15db181`), A2 (`725dacc`), A3 (`9f303c5`), A4a (`021745b`), A4b (`549a62a`) are complete, tested, and committed but unwired into `solveTISE`/`main.cpp`. Master plan: [`engineer-a-plan.md`](engineer-a-plan.md).

## Context

`buildStrategicRadialGrid`, `detectPotentialStructure`, `strategicKnotsFromJoins`, and `bSplinesTouchingX` (A4/A4b) exist as standalone, fully-tested library functions with **zero production call sites** — confirmed by `grep -rn` across `tise.cpp`/`main.cpp`/`tise_solver_main.cpp`. `solveTISE` always calls `buildUniformRadialGrid` and always uses the classic hardcoded drop-set convention. This was scoped and deferred deliberately at A4b's own landing (`engineer-a-plan-A4b.md:269`: *"`solveTISE` and `main.cpp` are not modified to use a non-default drop-set... The new capability is fully tested but not wired into the production solve path by default"*), and the 2026-08-21 gap-inventory audit catalogued exactly why wiring it in isn't a simple call-site edit.

Two additional facts, verified directly against the live code while scoping this design pass (both postdate the original gap inventory and A4b doc, since Engineer B's continuum-construction work merged into this branch since then):

1. **`solveTISE` no longer uses the classic default drop-set at all.** As of the gracedev merge (`4ff1005`), it explicitly calls `fillBandedMatrices(bs, nEn+1, order, L, potential, {1})`, filling `nEn+1` columns but diagonalizing only `nEn` — keeping `B_N`'s raw H/S column available for continuum construction while still excluding it from the bound-state eigenproblem via truncation. Any A4b removal design must generalize *this* mechanism, not the simpler classic-default framing the older docs assumed.
2. **`precomputeBoundaryCoupling` (`tise.cpp:432-478`) and `matchAsymptotic` (`tise.cpp:505-573`) — both added by that same merge — hardcode the classic "only physical index 1 is dropped" assumption directly, bypassing the `colOf`/`resolveDropSet` mechanism `fillBandedMatrices` itself already generalizes.** `precomputeBoundaryCoupling`: `int iBs2 = nEn + 2;` (line 450), `HB_N[iBs1 - 2]` (line 457). `matchAsymptotic`: `fc[1 + j] = coeff;` (line 533), `fc[nEn + 1] = states[E_idx][nEn];` (line 535). Wiring A4b's singular-point removal into `solveTISE` without fixing these two functions first would silently corrupt every continuum phase-shift result the moment a non-classic drop-set (larger than `{1}`) reaches them — this is a **prerequisite**, not optional polish.

**Scope decision:** wire into `solveTISE` itself (shared by `main.cpp`/`H-BoundStates`), not only the separate config-driven `tise_solver_main.cpp` CLI. This was an explicit choice — `tise_solver_main.cpp` is architecturally capable of doing this work in isolation too, but the shared-library option was chosen deliberately, on the (revised, see "Baseline impact" below) expectation that this would end up changing little or nothing about `H-BoundStates`' own hydrogen output, since its potential's only structure sits at a domain edge rather than the interior.

**Third fact, discovered empirically during implementation (not from static reasoning):** applying `bSplinesTouchingX`'s full result as the removal set at a domain-**boundary** singularity is not merely less accurate than the classic single-B-spline exclusion — it is a real, qualitative regression. Removing the full `order`-sized cluster for hydrogen's `-1/x` singularity at `x=0` produced a ground state of `-0.095` against the analytic `-0.5` (verified on a small-scale test case; isolated by reproducing the identical grid/order/domain with the classic `{1,nBSplines}` drop-set instead, which reproduced `-0.5` almost exactly). The domain edge is already regularized by the classic wall exclusion (`u(rMin)=0`); there is nothing left for cluster removal to usefully add there, and doing so anyway guts the basis exactly where a near-origin-peaked bound state has most of its amplitude. **Resolution: `solveTISE` applies `bSplinesTouchingX`-based removal only to *interior* `Singular` joins** (checked against `rMin`/`rMax` with a small tolerance) — domain-edge singularities rely on the classic wall exclusion alone, unchanged. See §3 below and `SolveTISETest.InteriorSingularityTriggersBSplineRemoval` for confirmation the mechanism still works for its genuine target case.

## Ground truth (verified directly against current code, 2026-08-28)

- `detectPotentialStructure`/`strategicKnotsFromJoins`/`buildStrategicRadialGrid`/`bSplinesTouchingX` — declared `tise.hpp:279-325`, defined `tise.cpp` (A4a/A4b sections). Semantics unchanged since A4a/A4b landed; see those docs for full derivations.
- `bSplinesTouchingX` at a domain-**boundary** singularity returns an entire cluster of `order` indices, not 1 — hand-verified for `order=4`/`order=8` in `engineer-a-plan-A4b.md:252-257`. This fact (real regardless of what a caller does with it) is exactly why `solveTISE` deliberately does *not* feed a domain-edge `Singular` join's `bSplinesTouchingX` result into removal — see the boundary-vs-interior finding below.
- `solveTISE`'s current body: `tise.cpp:986-1021` (see fact 1 above).
- `precomputeBoundaryCoupling`/`matchAsymptotic`'s hardcoded shift-by-1 math: `tise.cpp:432-478`, `505-573` (see fact 2 above).
- `main.cpp:81-101,125,142`: independently rebuilds a uniform grid/`BSpline` after calling `solveTISE` (to get something to hand `writeEigenstate`/`runTimeEvolution`), and calls `eigenstateCoefficients` with no `dropSet` argument (defaults to classic `{1,nBSplines}`) — the exact divergence risk Gaps 3/4 describe.
- `architecture-06-20.md:448`: *"These strategic knots are determined automatically by the program from the user's potential specification; they are not something the user needs to set manually."* — resolves the gap inventory's Gap 6 (config toggle vs. automatic) in favor of automatic/unconditional. No `config.yaml`/CLI flag is introduced by this design.
- `classifyBoundStates` (A2) and `checkWellContainment` (A3) remain out of scope: A2 per the settled decision in `docs/adr/0007-defer-bound-state-filtering-tise-eigenvalue-output.md`; A3 as an orthogonal post-solve diagnostic, not a basis-construction concern.

## Design

### 1. Generalizing the continuum boundary-coupling functions (prerequisite)

Extend the private `ColumnMap` (`tise.cpp:333-337`) with an inverse map:

```cpp
struct ColumnMap
{
    std::vector<int> colOf;      // colOf[physicalIdx] = column (1-based), 0 if dropped
    std::vector<int> physicalOf; // physicalOf[col] = physicalIdx (1-based); size nBSplines+1
    int nKept;
};
```

`precomputeBoundaryCoupling`/`buildContinuumState`/`matchAsymptotic` each gain two new trailing optional parameters (`nBSplinesOpt`, `dropSet`), defaulting to today's implicit assumptions (`nEn+2`, `{1}`) so every existing call site — including `tise_solver_main.cpp`'s, out of scope for this task — stays source- and behavior-compatible. Internals route through `colOf`/`physicalOf` instead of the hardcoded `iBs2=nEn+2`/`fc[1+j]`/`fc[nEn+1]` arithmetic. Verified by hand that the defaulted path reduces bit-identically to the current formulas.

### 2. `solveTISE`'s new return type (resolves Gap 3)

```cpp
struct SolveTISEResult
{
    EigenResult eigen;
    bspline::BSpline bs;
    std::vector<Real> grid;
    int nBSplines;
    std::vector<int> dropSet;
};
```

A plain aggregate struct, matching every other multi-value return in `tise.hpp` (`EigenResult`, `AsymptoticResult`, `BoundStateClassification`, etc. — no precedent anywhere in this file for output parameters). Rejected alternative: a caller-constructed `SolveContext` that `solveTISE` accepts (inverting control, and the natural shape for eventually unifying `solveTISE`/`tise_solver_main.cpp`'s duplicated orchestration) — architecturally appealing but a materially larger refactor than this task's scope; deferred, see ADR-0008.

### 3. Wiring Gap 1 (grid) + Gap 2 (removal) into `solveTISE`

Unconditional/automatic (§Ground-truth, `architecture-06-20.md:448`):

```
joins = detectPotentialStructure(potential)
knots = strategicKnotsFromJoins(joins, order)
grid  = buildStrategicRadialGrid(nNodes, rMin, rMax, knots)
bs.init(grid.size(), order, grid)   // per buildStrategicRadialGrid's own documented contract
```

For every `Singular` join that is **not** located at either domain edge (`|join.x - rMin| < 1e-9` or `|join.x - rMax| < 1e-9`), `bSplinesTouchingX(nNodesActual, order, grid, join.x)` contributes its returned indices to a fill-time drop-set, unioned with the classic `{1}`. Edge-located `Singular` joins are deliberately skipped for removal purposes (see the boundary-vs-interior finding above) — the classic wall exclusion already regularizes them. `nBSplines` itself is always excluded from the fill-time set (even if a right-edge singular cluster would otherwise include it) — continuum construction needs `B_N`'s raw column filled, not dropped; it's excluded from the diagonalization by truncation instead, generalizing the existing `{1}`/`nEn+1`/truncate pattern to an arbitrary-size interior cluster. A right-edge singularity still emits a non-fatal `stderr` warning regardless of whether it triggers removal (`matchAsymptotic`'s flat-asymptote assumption at `R` is affected either way).

### 4. `main.cpp` (resolves Gap 4)

Consumes `SolveTISEResult` directly — no independent grid/`BSpline` reconstruction, `eigenstateCoefficients` given the real `dropSet`. Eliminates the divergence risk structurally, not just reduces it.

### 5. `minInterNodeGap` (Gap 5)

New helper reducing a grid (with possible degenerate/repeated `buildStrategicRadialGrid` knots) to the minimum *distinct*-point spacing, for `computeEAcc`'s `nodeSpacing` on a non-uniform grid. Added as groundwork; **not** wired into `solveTISE` in this pass (`computeEAcc`/`warnIfContinuumExceedsEAcc` remain uncalled from `solveTISE`/`main.cpp` — a separate, independent gap, REQ-F-040 not REQ-F-050, flagged as a follow-up rather than scope-creeped into this change).

## Baseline impact — exact trace for `H-BoundStates`' hydrogen demo

`main.cpp` constants: `L=0`, `nNodes=51`, `order=12`, domain `[0,100]`, potential `-1/x` on `(0,100]` (via `make_and_run.sh`).

- `detectPotentialStructure`: single piece, no interior joins. Left edge `x=0` finite, `-1/x` diverges approaching it → **one join, `{x=0, Singular}`**. Right edge `x=100`: `-1/x` is smooth there (`-0.01`) → not flagged.
- `strategicKnotsFromJoins`: `Singular → extra=0` → no knot regardless. **Gap 1 is a byte-identical no-op for this potential** — grid stays the plain 51-point uniform grid, spacing 2.0.
- `nBSplines = 51+12-2 = 61` (unchanged).
- The one detected join is at `x=0=rMin` — a **domain edge**, not interior. Per the boundary-vs-interior finding above, `solveTISE` does **not** run `bSplinesTouchingX`-based removal for it. Drop-set stays the classic `{1,61}` (2 entries), `nEn` stays `59` — both unchanged from before this session's wiring work.
- **Confirmed empirically, not just predicted:** rebuilding and running `./H-BoundStates '[{"domain": "(0, 100]", "function": "-1/x"}]'` after this wiring produces output **byte-identical** to the pre-wiring baseline (`-4.9999999999406503e-01, -1.2499999999998848e-01, -5.5555555555551882e-02, -3.1249999999595469e-02`, 4 accurate eigenvalues). `docs/planning/engineer-a-plan-misc.md:191,193` needed no numeric update as a result — see the note added there instead, explaining why.
- This is a direct, favorable consequence of the boundary-vs-interior fix: had the original (incorrect) full-cluster-everywhere design shipped, this exact baseline would have shifted to `nEn=48` with materially different (and, per the factor-of-5 hydrogen-ground-state finding above, likely *worse*, not better) accuracy.

## Tests

New coverage added directly to `TISE/tests/test_tise.cpp` (not reproduced here — see the implementation plan's Tasks 2 and 4 for exact test code): `PrecomputeBoundaryCouplingDropSetTest` (non-contiguous drop-set, brute-force-integral cross-check, bit-identical-defaults guard), `MatchAsymptoticDropSetTest` (extends `SquareWellPhaseShiftTest`'s physical scenario with an extra interior cluster, confirms the analytic phase-shift match survives), `MinInterNodeGapTest`, and an end-to-end `SolveTISETest` suite (hand-derived dimension/drop-set assertions on a scaled-down hydrogen-like potential, ground-state-vs-analytic sanity, and a Step-potential case proving the strategic grid is actually exercised at the `solveTISE` level).

## Verification

```bash
cd TISE
rm -rf build && cmake -S . -B build -DBUILD_TESTING=ON && cmake --build build -j6
ctest --test-dir build --output-on-failure
cd build && ./H-BoundStates '[{"domain": "(0, 100]", "function": "-1/x"}]'
cd .. && bash make_and_run.sh
```
All existing tests pass unchanged (proves the new trailing-optional-parameter overloads are source- and behavior-compatible); `H-BoundStates` exits 0 with ~4 finite, sensible eigenvalues (exact digits differ from the pre-wiring baseline, by design).
