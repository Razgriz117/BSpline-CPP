# TISE Engineer-A — Task A4a Implementation Plan (active slice)

> Scoped to **task A4a only** (strategic node placement, required core — detect step/stitched-kink/singular joins, build a grid with degenerate knots). A4b (generalizing the hardcoded "drop `B_1`/`B_N`" convention) is an explicit stretch goal, attempted only once A4a is solid — **not** part of this plan. A1 (`15db181`), A2 (`725dacc`), A3 (`9f303c5`) are complete and committed. Master plan: [`engineer-a-plan.md`](engineer-a-plan.md) §"Section A4". This plan elaborates that section to the same rigor as A1-A3, fixing two gaps its original sketch didn't fully work out (below).

## Context

`docs/planning/tise-task-breakdown.md` (§"A4. Strategic node placement (REQ-F-050)") and `docs/SDD.md` (REQ-F-050: *"place B-spline knots using uniform spacing augmented by strategic potential-driven nodes as the default collocation scheme"*) specify: detect potential structure (steps, stitched-derivative kinks, singularities) from the piecewise `potential` map, and insert degenerate knots at step/kink joins so the B-spline basis can represent the resulting wavefunction discontinuity in `ψ''`/`ψ'''`. Singular-potential B-spline *removal* (the table's 4th row) is detection-only here — actually removing B-splines is A4b, out of scope.

**Done when** (task-breakdown's exact wording): *"a potential with a domain-boundary discontinuity (e.g. the existing 'particle in a box with a rectangular barrier' example from the config schema spec) produces a grid with degenerate knots at the boundary, and eigenvalue accuracy improves (smaller `eigenvalueError`-style residual) relative to a purely uniform grid at the same node count."*

Per the user's established cadence for this task set (one task implemented at a time, paused for review after each — see A1/A2/A3), this is the next task after A3.

**Environment/repo state (re-verified 2026-08-14, no drift since A3 landed):** `HEAD` is `9f303c5`, working tree clean, nothing landed since A3. No Section-0-style setup needed.

**Prior-session note:** A design for this task was sketched in chat on 2026-08-08 (session `cd3de9fa`) after a design subagent hit a 64,000-output-token API error and terminated early. That sketch got the shape right (reuse `classifySequenceConvergence`, Richardson-extrapolated one-sided value/slope jumps, `order-3`/`order-4` multiplicity formula) but was never written to a durable plan doc, never had its insertion points/line numbers verified against current code, and — as this plan's ground-truth pass below found — missed a real correctness bug (see "Gap 3" below). This document supersedes that sketch; nothing here is copied from it without independent verification against current code.

## Ground truth verified directly against current code (not just prose)

- **`BSpline::init(int numberOfNodes, int order, const std::vector<Real> &gridIn)`** — `BSpline.cpp:175-256`. Two facts that directly drive this task's design:
  - `nBSplines_ = numberOfNodes + order - 2` (`BSpline.cpp:187`) — a degenerate knot (the same x value repeated in `gridIn`) increases `numberOfNodes` and therefore `nBSplines_` by exactly the number of repeats. This is *how* degenerate knots relax continuity: more coincident knots at one x → more independent B-splines whose support boundary sits there → less forced continuity.
  - **Confirmed silent-truncation risk**: `numberOfNodes` and `gridIn.size()` are independent parameters (`BSpline.cpp:190-193`: only checks `gridIn.size() >= numberOfNodes`, then `BSpline.cpp:212-216` only reads `gridIn[0..numberOfNodes-1]`). If a caller builds a longer grid (with repeated knots for multiplicity) but passes the *original* `nNodes` as `numberOfNodes`, the extra knots are silently dropped — no error, no warning. `buildStrategicRadialGrid`'s contract must make the caller pass its own returned `.size()`, documented prominently, with a dedicated test.
  - **Monotonicity check** (`BSpline.cpp:196-202`) only rejects strict decreases (`gridIn[i] < gridIn[i-1]`) — repeated values pass. No changes needed to `BSpline.cpp` itself; repeated knots are already handled correctly by the existing math, just never previously exercised by a non-uniform caller.
- **`inInterval`** (`tise.cpp:50-88`) — regex `^\s*([\[\(])\s*(bound)\s*,\s*(bound)\s*([\]\)])\s*$` (case-insensitive; bound = signed decimal or `inf`/`infinity`), plus a `parseBound` lambda handling the infinity cases and `std::stod` for finite ones. This regex+parsing logic is what gets factored out into a new file-local `parseInterval` helper (see Design), reused by `detectPotentialStructure` to get numeric piece bounds sorted correctly — `std::map`'s native iteration order is lexicographic on the domain *string* (e.g. `"[10,...)"` sorts before `"[5,...)"`), which would silently misorder joins if pieces were processed in raw map-iteration order.
- **`classifySequenceConvergence`** (`tise.hpp:52-54` already anticipates this: *"Shared numeric core for asymptote classification (REQ-F-030) and, later, singular-point detection for strategic node placement (REQ-F-050)"*). Read in full (`tise.cpp:118-153`+): flatness pre-check (`maxAbsDV <= 1e-10 + 1e-9*maxAbsV` → `isFlat`), then a `shrinkFactor = |dV[N-2]| / max(|dV[mid]|, 1e-300)` divergence check (`shrinkFactor >= 0.5` → `isDivergent`). Critically, **this divergence test is direction-agnostic** — it only compares finite-difference magnitudes across the array `V[0..N-1]` in the order given; it has no built-in assumption about whether the samples move toward or away from a reference point. That's what makes reusing it for an *inward*-shrinking probe (toward a potential singularity) instead of A1's *outward*-growing probe (toward infinity) valid: for a genuine singularity, V grows as samples approach it, so `dV` grows near the end of the array exactly the way A1's "doesn't converge" case does — same `shrinkFactor >= 0.5` logic correctly fires. Hand-verified by tracing the Coulomb-origin test case below.
- **`classifyAsymptote`** (`tise.cpp:207-274`) — the sampling pattern to mirror: `kNumSamples=16`, `kRatio=4.0`, `scale = max(|reference|, 1.0)`, samples at `reference + sign*scale*ratio^k` for `k=0..15`.
- **`tise.cpp` is 445 lines.** `checkWellContainment` (A3) ends at line 384; `analyticHydrogenEnergy` starts at line 386. **A4a's new code inserts in that gap** (after line 384, before line 386), continuing the "insert after the previous task, before `analyticHydrogenEnergy`" convention A2/A3 both used. Separately, the `parseInterval` refactor modifies `inInterval` **in place** at its current location (lines 50-88), not at the end-of-file insertion point.
- **`tise.hpp` is 204 lines.** `checkWellContainment`'s declaration ends at line 175; `analyticHydrogenEnergy` is declared at line 178. New A4a declarations insert in that gap.
- **`tests/test_tise.cpp` is 1039 lines.** A3's last test (`WellContainmentSyntheticTest.DerivativeJustAboveToleranceIsFlagged`) ends at line 973; the `writeEigenstate` banner starts at line 975. **A4a's new tests insert in that gap.**
- **`Real` = `double`** (`BSpline.hpp:12`, aliased in `tise.hpp:12`). `inInterval`/`parseInterval` keep using raw `double` (unchanged, behavior-preserving refactor of pre-existing code); new A4a-specific code uses `Real`, matching A1-A3 convention for new physics-domain code.
- **No new `#include`s needed anywhere.** `tise.cpp` already includes `<algorithm>`, `<cmath>`, `<regex>` (lines 3, 5, 13) — covers `std::sort`, `std::lower_bound`, `std::min`/`std::max`, `std::abs`, `std::isfinite`, `std::pow`. `tests/test_tise.cpp` already includes `<algorithm>` (line 5) — covers `std::count`, `std::is_sorted` used in the new tests.
- **muparser `^` power syntax confirmed supported** — already used in existing fixtures, e.g. `tests/test_tise.cpp:316,328`: `{{"(0, inf)", "1/x^1.5"}}`. Needed for the Coulomb-origin test's `"-1/x + 1/x^2"`.
- **Box+barrier example confirmed sourced from `docs/TDSE-original-design/2026-06-28-config-yaml-schema-design.md:101-104`** ("Particle in a box with a rectangular barrier"): `{"[0, 5)": "0"}, {"[5, 6]": "10"}, {"(6, 10]": "0"}` — domain `[0,10]`. This is the task-breakdown's own literal "Done when" example; used for both the `Step` detection tests and the accuracy test.

## Two real gaps found and fixed during this planning pass

**Gap 1 (silent truncation, from a prior chat sketch — see "Prior-session note" above):** `buildStrategicRadialGrid`'s returned grid must be passed to `BSpline::init` with `numberOfNodes = returned.size()`, not the original `nNodes` argument, or the extra knots get silently dropped (see "Ground truth" above). Fixed via prominent header-comment contract + a dedicated `ReturnedSizeAccountsForAllInsertions` test, and by having the accuracy test itself pass `strategicGrid.size()` (not `nCoarse`) to `bs.init`.

**Gap 2 (knot-coincidence, also from the prior sketch):** a strategic knot's x isn't guaranteed to already sit on the uniform base grid. Fixed: `buildStrategicRadialGrid` splices in a new grid point via `std::lower_bound` + `insert` when the join x isn't already present (tolerance `1e-12` for "already present"), *then* adds `extraMultiplicity` more copies on top.

**Gap 3 (found independently during this pass, not present in the prior chat sketch — a real correctness bug, not just a documentation issue):** the prior sketch's plan to "mirror A1's `classifyAsymptote` sampling exactly, just inverted" is unsafe as literally stated. A1's `scale = max(|reference|, 1.0)` is fine for its own use (samples move *outward* toward infinity, so there's no upper bound to worry about) — but for A4a's *inward*-shrinking probe, the largest sample offset (`h0`, at `k=0`) must stay inside whichever piece is being sampled, or `evaluateFunction` throws (`tise.cpp:109-112`, "Function domain does not cover x"). Concretely: for the box+barrier potential's join at `x=6` (between `[5,6]` and `(6,10]`), checking the *left* piece `[5,6]` (width 1) with the naive `scale = max(|6|,1) = 6` steps as far as `x = 6-6 = 0` — three pieces away from where it should be sampling, hard outside `[5,6]`, guaranteed exception. **Fix:** clamp `h0` by the sampled piece's own width: `h0 = isfinite(width) ? min(scale, 0.25*width) : scale`. Verified safe for all four canonical test cases by hand (see Design below); verified this clamp doesn't change behavior for the wide-piece Coulomb-origin case (width 100 ≫ scale 1, clamp is a no-op there).

## Design

### `parseInterval` refactor (`tise.cpp`, modifies existing lines 50-88 in place — behavior-preserving)

Extracts `inInterval`'s regex/bound-parsing into a file-local (anonymous-namespace) helper that returns numeric bounds instead of just a bool, so `detectPotentialStructure` can reuse it without duplicating a nontrivial regex. Not exposed in `tise.hpp` — no caller outside this file needs the numeric bounds directly.

```cpp
// given a string, e.g. [0, 20), return true if input x falls in the bounds and false otherwise
// Bounds may be finite numbers or (+/-)inf/infinity (case-insensitive); brackets
// select inclusive ([, ]) vs. exclusive ((, )) endpoints.
//
// Regex/bound-parsing logic factored into parseInterval() (below) so A4a's
// detectPotentialStructure can reuse it to get numeric piece bounds without
// duplicating this nontrivial regex. File-local only -- not exposed in
// tise.hpp, since no caller outside this file needs the numeric bounds
// directly.
namespace
{
struct ParsedInterval
{
    double lower;
    double upper;
    bool lowerInclusive;
    bool upperInclusive;
};

ParsedInterval parseInterval(const std::string& interval)
{
    static const std::regex re(
        R"(^\s*([\[\(])\s*(-?(?:\d+(?:\.\d+)?)|[+-]?(?:inf|infinity))\s*,\s*(-?(?:\d+(?:\.\d+)?)|[+-]?(?:inf|infinity))\s*([\]\)])\s*$)",
        std::regex_constants::icase
    );

    std::smatch m;
    if (!std::regex_match(interval, m, re)) {
        throw std::runtime_error("Invalid interval: " + interval);
    }

    auto parseBound = [](const std::string& s) {
        std::string t;
        t.reserve(s.size());
        for (char c : s)
            t += std::tolower(static_cast<unsigned char>(c));

        if (t == "inf" || t == "+inf" ||
            t == "infinity" || t == "+infinity")
            return std::numeric_limits<double>::infinity();

        if (t == "-inf" || t == "-infinity")
            return -std::numeric_limits<double>::infinity();

        return std::stod(s);
    };

    ParsedInterval out;
    out.lowerInclusive = m[1] == "[";
    out.upperInclusive = m[4] == "]";
    out.lower = parseBound(m[2]);
    out.upper = parseBound(m[3]);
    return out;
}
} // namespace

bool inInterval(double x, const std::string& interval)
{
    const ParsedInterval iv = parseInterval(interval);
    bool leftOK = iv.lowerInclusive ? (x >= iv.lower) : (x > iv.lower);
    bool rightOK = iv.upperInclusive ? (x <= iv.upper) : (x < iv.upper);
    return leftOK && rightOK;
}
```

Behavior-preserving: same regex, same `parseBound`, same comparison logic — `inInterval`'s existing tests are unaffected and need no changes.

### `tise.hpp` (insert after line 175, before line 177's `analyticHydrogenEnergy` comment)

```cpp
// === A4a: strategic node placement (REQ-F-050), required core only ===
// A4b (generalizing the hardcoded "drop B_1/B_N" convention so singular-
// potential B-spline removal is more than detection-only) is deferred --
// see docs/planning/engineer-a-plan.md Section A4. Singular joins are
// detected here but not remediated by this task.

// One of the four required-treatment categories from the task's source
// table (docs/planning/tise-task-breakdown.md, "A4. Strategic node
// placement"): Step (V itself jumps), StitchedKink (V continuous, V' jumps),
// Singular (V diverges, e.g. 1/r), or Continuous (no special treatment).
enum class JoinType
{
    Continuous,
    Step,
    StitchedKink,
    Singular
};

// A detected join (interior piece-to-piece boundary) or flagged domain edge
// in the potential's piecewise structure, at position x.
struct DetectedJoin
{
    Real x;
    JoinType type;
};

// A single knot location that should receive `extraMultiplicity` additional
// degenerate copies beyond the ordinary grid point already there (or
// spliced in if not already present -- see buildStrategicRadialGrid).
struct StrategicKnot
{
    Real x;
    int extraMultiplicity;
};

// Detect step/stitched-kink/singular structure in a piecewise potential's
// domain joins and finite outer domain edges. Interior joins (two pieces
// meet) are always emitted, even Continuous ones -- so a Continuous join's
// presence in the output means "detected, not actionable", not "silently
// dropped". Finite domain edges (a piece's own outer boundary with no
// neighbor, e.g. x=0 in {"(0,100]": "-1/x+1/x^2"}) are only emitted if
// singular there -- an ordinary box wall isn't actionable for this task.
// See docs/planning/engineer-a-plan-A4.md for the detection algorithm and
// threshold derivations.
std::vector<DetectedJoin> detectPotentialStructure(const std::map<std::string, std::string> &potential);

// Convert detected joins into strategic knots, per the multiplicity formula
// derived in docs/planning/engineer-a-plan-A4.md: Step -> order-3,
// StitchedKink -> order-4 (both clamped at 0 -- a case this codebase's
// order values, 4 and 8, only actually reach for StitchedKink at order=4).
// Singular/Continuous produce no knot (Singular's remediation -- B-spline
// removal -- is A4b, out of scope here).
std::vector<StrategicKnot> strategicKnotsFromJoins(const std::vector<DetectedJoin> &joins, int order);

// Build a radial grid combining a uniform base (buildUniformRadialGrid) with
// degenerate knots at each StrategicKnot's location: extraMultiplicity
// repeated copies, splicing in a brand-new grid point first if `x` doesn't
// already coincide (within 1e-12) with a uniform grid point. Returns the
// FULL combined grid.
//
// CONTRACT: callers MUST pass the returned vector's own .size() as
// BSpline::init's `numberOfNodes` argument, not the original nNodes --
// BSpline::init only reads the first numberOfNodes entries of the grid it's
// given and silently drops the rest (verified against BSpline.cpp:190-216).
std::vector<Real> buildStrategicRadialGrid(int nNodes, Real rMin, Real rMax,
                                            const std::vector<StrategicKnot> &knots);
```

### `tise.cpp` (insert after line 384, before line 386's `analyticHydrogenEnergy`)

```cpp
// === A4a: strategic node placement (REQ-F-050), required core only ===
namespace
{
// Shrinking-offset probe toward x0 from one side, reusing A1's numeric core
// (classifySequenceConvergence) but inverted: offsets shrink toward x0
// instead of growing toward infinity -- valid because that function's
// divergence test only compares finite-difference magnitudes across the
// given array, with no built-in assumption about sample direction (see
// docs/planning/engineer-a-plan-A4.md). `sign` is -1 to sample the piece
// left of x0, +1 for the piece right of x0; `pieceWidth` is that piece's
// own width, used to clamp the starting offset so the probe can never step
// outside a narrow piece (e.g. the box+barrier's width-1 barrier slab) --
// A1's own scale formula alone isn't safe for this inward direction.
bool isSingularApproaching(const std::map<std::string, std::string> &potential,
                            Real x0, Real sign, Real pieceWidth)
{
    constexpr int kNumSamples = 16;
    constexpr Real kRatio = 4.0;
    const Real scale = std::max(std::abs(x0), 1.0);
    const Real h0 = std::isfinite(pieceWidth) ? std::min(scale, 0.25 * pieceWidth) : scale;

    std::vector<Real> V(kNumSamples);
    for (int k = 0; k < kNumSamples; ++k)
        V[k] = evaluateFunction(potential, x0 + sign * h0 / std::pow(kRatio, k));

    return classifySequenceConvergence(V, kRatio).isDivergent;
}

// One-sided value estimate at x0, approaching from `sign` (-1 = left piece,
// +1 = right piece). Richardson-extrapolated over step sizes h and h/2 to
// cancel the O(h) bias from necessarily sampling strictly inside the piece
// (open piece boundaries, e.g. "[0,5)", never include x0=5 itself):
// V(x0+sign*h) = V0 + sign*h*V0' + O(h^2), so 2*V(h/2) - V(h) = V0 + O(h^2).
Real oneSidedValue(const std::map<std::string, std::string> &potential,
                    Real x0, Real sign, Real h)
{
    const Real gH  = evaluateFunction(potential, x0 + sign * h);
    const Real gH2 = evaluateFunction(potential, x0 + sign * h / 2.0);
    return 2.0 * gH2 - gH;
}

// One-sided slope (dV/dx) estimate at x0, same side convention. Builds a
// standard one-sided 2nd-order finite-difference derivative,
// f'(a) ~ [-3f(a) + 4f(a+h) - f(a+2h)] / (2h), centered at a = x0+sign*h
// (not x0 itself, since x0 may not be in this piece's domain) -- that
// off-centering is an O(h) bias, which the same Richardson trick (evaluate
// at h and h/2, then 2*D(h/2) - D(h)) cancels, leaving O(h^2) error --
// negligible at h~1e-4 against the StitchedKink threshold below.
Real oneSidedSlope(const std::map<std::string, std::string> &potential,
                    Real x0, Real sign, Real h)
{
    auto forwardDeriv = [&](Real step) {
        const Real g1 = evaluateFunction(potential, x0 + sign * step);
        const Real g2 = evaluateFunction(potential, x0 + sign * 2.0 * step);
        const Real g3 = evaluateFunction(potential, x0 + sign * 3.0 * step);
        return sign * (-3.0 * g1 + 4.0 * g2 - g3) / (2.0 * step);
    };
    const Real dH  = forwardDeriv(h);
    const Real dH2 = forwardDeriv(h / 2.0);
    return 2.0 * dH2 - dH;
}
} // namespace

std::vector<DetectedJoin> detectPotentialStructure(const std::map<std::string, std::string> &potential)
{
    struct PieceInfo { Real lower, upper; };
    std::vector<PieceInfo> pieces;
    pieces.reserve(potential.size());
    for (const auto &entry : potential)
    {
        const ParsedInterval iv = parseInterval(entry.first);
        pieces.push_back({iv.lower, iv.upper});
    }
    // Sort by numeric lower bound -- std::map's native order is lexicographic
    // on the domain *string* (e.g. "[10,...)" sorts before "[5,...)"), which
    // would silently misorder joins if used directly.
    std::sort(pieces.begin(), pieces.end(),
              [](const PieceInfo &a, const PieceInfo &b) { return a.lower < b.lower; });

    constexpr Real kValueJumpTol = 1e-6; // O(1)-O(10) potential magnitudes in
                                          // this codebase; a genuine Step is
                                          // an O(1) jump, ~6 orders of margin.
    constexpr Real kSlopeJumpTol = 1e-4; // one-sided slope estimate's O(h^2)
                                          // truncation error is ~1e-8 at
                                          // h=1e-4, ~4 orders below this.
    constexpr Real kFdStep = 1e-4;

    auto clampedStep = [&](Real width) {
        return std::isfinite(width) ? std::min(kFdStep, 0.1 * width) : kFdStep;
    };

    auto classifyJoin = [&](Real x0, Real leftWidth, Real rightWidth) {
        if (isSingularApproaching(potential, x0, -1.0, leftWidth) ||
            isSingularApproaching(potential, x0, +1.0, rightWidth))
            return JoinType::Singular;

        const Real hL = clampedStep(leftWidth);
        const Real hR = clampedStep(rightWidth);
        const Real vL = oneSidedValue(potential, x0, -1.0, hL);
        const Real vR = oneSidedValue(potential, x0, +1.0, hR);
        if (std::abs(vL - vR) > kValueJumpTol)
            return JoinType::Step;

        const Real sL = oneSidedSlope(potential, x0, -1.0, hL);
        const Real sR = oneSidedSlope(potential, x0, +1.0, hR);
        if (std::abs(sL - sR) > kSlopeJumpTol)
            return JoinType::StitchedKink;

        return JoinType::Continuous;
    };

    std::vector<DetectedJoin> joins;

    for (std::size_t i = 0; i + 1 < pieces.size(); ++i)
    {
        const Real x0 = pieces[i].upper;
        const Real leftWidth = pieces[i].upper - pieces[i].lower;
        const Real rightWidth = pieces[i + 1].upper - pieces[i + 1].lower;
        joins.push_back({x0, classifyJoin(x0, leftWidth, rightWidth)});
    }

    // Global domain edges: only flagged if singular (an ordinary box wall
    // isn't interesting to strategic node placement). Also probes each
    // piece's own endpoints, not just inter-piece joins, so a single-piece
    // potential like {"(0, inf)": "-1/x + 1/x^2"} still gets its x=0
    // singularity flagged even though it has no "join" at all.
    if (!pieces.empty())
    {
        const PieceInfo &first = pieces.front();
        if (std::isfinite(first.lower))
        {
            const Real width = first.upper - first.lower;
            if (isSingularApproaching(potential, first.lower, +1.0, width))
                joins.push_back({first.lower, JoinType::Singular});
        }
        const PieceInfo &last = pieces.back();
        if (std::isfinite(last.upper))
        {
            const Real width = last.upper - last.lower;
            if (isSingularApproaching(potential, last.upper, -1.0, width))
                joins.push_back({last.upper, JoinType::Singular});
        }
    }

    return joins;
}

std::vector<StrategicKnot> strategicKnotsFromJoins(const std::vector<DetectedJoin> &joins, int order)
{
    std::vector<StrategicKnot> knots;
    knots.reserve(joins.size());
    for (const auto &j : joins)
    {
        int extra = 0;
        switch (j.type)
        {
            case JoinType::Step:         extra = std::max(0, order - 3); break;
            case JoinType::StitchedKink: extra = std::max(0, order - 4); break;
            case JoinType::Singular:
            case JoinType::Continuous:   extra = 0; break;
        }
        if (extra > 0)
            knots.push_back({j.x, extra});
    }
    return knots;
}

std::vector<Real> buildStrategicRadialGrid(int nNodes, Real rMin, Real rMax,
                                            const std::vector<StrategicKnot> &knots)
{
    std::vector<Real> grid = buildUniformRadialGrid(nNodes, rMin, rMax);

    std::vector<StrategicKnot> sorted = knots;
    std::sort(sorted.begin(), sorted.end(),
              [](const StrategicKnot &a, const StrategicKnot &b) { return a.x < b.x; });

    for (const auto &knot : sorted)
    {
        auto it = std::lower_bound(grid.begin(), grid.end(), knot.x);
        const bool alreadyPresent = (it != grid.end()) && (std::abs(*it - knot.x) < 1e-12);

        if (!alreadyPresent)
            it = grid.insert(it, knot.x); // splice a new base point at knot.x

        grid.insert(it + 1, knot.extraMultiplicity, knot.x); // extra degenerate copies
    }

    return grid;
}
```

**Multiplicity formula derivation** (unchanged from the master plan, re-verified independently): for a B-spline basis of order `k`, a knot of *total* multiplicity `m` gives continuity `C^{k-1-m}`. To make ψ's derivative of order `n` the first discontinuous one, need continuity `C^{n-1}`, i.e. `m = k-n`. An ordinary grid point already has base multiplicity 1, so *extra* multiplicity is `m-1 = k-n-1`. Step forces a ψ'' jump (`n=2`, from `ψ'' = 2m(V-E)ψ`) → `extra = order-3`. StitchedKink forces a ψ''' jump (`n=3`, from differentiating once more) → `extra = order-4`. Cross-check: a delta potential forces a ψ' jump (`n=1`) → `extra = order-2`, one less than full multiplicity (continuity `C^0`), matching the standard delta-matching condition — confirms the formula independent of the (missing-from-repo) source PDF.

### Hand-traced verification of all 4 canonical `detectPotentialStructure` cases

- **`StepAtBothBarrierEdges`** `{{"[0,5)","0"},{"[5,6]","10"},{"(6,10]","0"}}`: join x=5 (pieces "0"/"10", widths 5/1): neither side singular (constant functions → `isFlat`, not divergent); `vL≈0`, `vR≈10`, jump=10 ≫ 1e-6 → `Step`. Join x=6 (pieces "10"/"0", widths 1/4): same reasoning, jump=10 → `Step`. Domain edges: x=0 into width-5 piece "0" (not singular), x=10 into width-4 piece "0" (not singular) → neither edge flagged. Result: exactly `{x=5,Step}, {x=6,Step}` — matches.
- **`StitchedKinkAtCornerJoin`** `{{"[0,5)","x"},{"[5,10]","10-x"}}`: join x=5 (widths 5/5): not singular; `vL≈5`, `vR≈10-5=5`, jump≈0 → not Step; `sL≈1` (d(x)/dx), `sR≈-1` (d(10-x)/dx), jump=2 ≫ 1e-4 → `StitchedKink`. Edges: x=0 into "x" (not singular), x=10 into "10-x" (value 0, not singular) → not flagged. Result: exactly `{x=5,StitchedKink}` — matches.
- **`SingularAtCoulombOrigin`** `{{"(0,100]","-1/x+1/x^2"}}`: single piece, no interior joins. Left edge x=0: `scale=max(0,1)=1`, width=100, `h0=min(1,25)=1`; samples `V[k]=V(1/4^k)` → dominated by `1/x^2`, grows explosively as k→15 (x→~2.3e-10, V→~1.9e19) → `dV` grows near the end → `shrinkFactor≫0.5` → `isDivergent=true` → `Singular`. Right edge x=100: approaching from the left (sign=-1) into the same piece, `-1/x+1/x^2` is finite and smooth at x=100 → not singular → not flagged. Result: exactly `{x=0,Singular}` — matches.
- **`ContinuousJoinNotFlagged`** `{{"[0,5)","x*x"},{"[5,10]","x*x"}}`: join x=5: not singular (x² finite everywhere near 5); `vL≈vR≈25`, jump≈0 → not Step; `sL≈sR≈10` (d(x²)/dx=2x=10), jump≈0 → not StitchedKink → `Continuous`. Edges: x=0 and x=10 into "x*x", neither singular → not flagged. Result: exactly `{x=5,Continuous}` — one join present (detected, not actionable), matches the test's intent.

### Scope: standalone, unwired (matches A1/A2/A3 precedent)

`detectPotentialStructure`/`strategicKnotsFromJoins`/`buildStrategicRadialGrid` land as tested standalone functions only. `solveTISE` (`tise.cpp:428-443`) is **not** modified to call them by default — it keeps calling `buildUniformRadialGrid` exactly as today. REQ-F-050's "as the default collocation scheme" phrasing is aspirational for the eventual wired system; wiring it into the actual solve pipeline's call site is a later integration concern, same framing A1's asymptote classifier and A3's containment diagnostic already used for their own call-site wiring. The "Done when" criterion's demonstration of `BSpline::init` accepting a strategic grid happens inside the accuracy test itself (which does call `bs.init` with both grids), not in production code.

### Tests to add (`tests/test_tise.cpp`, insert after line 973, before the `writeEigenstate` banner at line 975)

```cpp
// ---------------------------------------------------------------------------
// detectPotentialStructure
// ---------------------------------------------------------------------------

TEST(DetectPotentialStructureTest, StepAtBothBarrierEdges)
{
    std::map<std::string, std::string> potential = {
        {"[0,5)", "0"}, {"[5,6]", "10"}, {"(6,10]", "0"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 2u);
    std::sort(joins.begin(), joins.end(),
              [](const tise::DetectedJoin &a, const tise::DetectedJoin &b) { return a.x < b.x; });
    EXPECT_DOUBLE_EQ(joins[0].x, 5.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Step);
    EXPECT_DOUBLE_EQ(joins[1].x, 6.0);
    EXPECT_EQ(joins[1].type, tise::JoinType::Step);
}

TEST(DetectPotentialStructureTest, StitchedKinkAtCornerJoin)
{
    std::map<std::string, std::string> potential = {
        {"[0,5)", "x"}, {"[5,10]", "10-x"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 5.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::StitchedKink);
}

TEST(DetectPotentialStructureTest, SingularAtCoulombOrigin)
{
    std::map<std::string, std::string> potential = {
        {"(0, 100]", "-1/x + 1/x^2"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 0.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Singular);
}

TEST(DetectPotentialStructureTest, ContinuousJoinNotFlagged)
{
    std::map<std::string, std::string> potential = {
        {"[0,5)", "x*x"}, {"[5,10]", "x*x"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 5.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Continuous);
}

TEST(DetectPotentialStructureTest, SingularAtInteriorJoin)
{
    // 1/(1-x) diverges approaching x=1 from the left; exercises
    // classifyJoin's own Singular branch (distinct from the domain-edge
    // singularity path SingularAtCoulombOrigin covers).
    std::map<std::string, std::string> potential = {
        {"[0,1)", "1/(1-x)"}, {"[1,2]", "0"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 1.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Singular);
}

TEST(DetectPotentialStructureTest, SingularAtRightDomainEdge)
{
    // 1/(100-x) diverges approaching x=100 from the left; exercises the
    // last-piece (right) global-edge singularity push, distinct from
    // SingularAtCoulombOrigin's left-edge case.
    std::map<std::string, std::string> potential = {
        {"[0, 100)", "1/(100-x)"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 100.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Singular);
}

TEST(DetectPotentialStructureTest, HandlesInteriorJoinWithInfiniteWidthPiece)
{
    // Right piece [5, inf) has infinite width, exercising
    // isSingularApproaching's isfinite(pieceWidth)==false fallback (falls
    // back to the plain scale=max(|x0|,1) probe, matching A1's own formula,
    // since there's no finite piece boundary to clamp against). -1/x^2 is
    // smooth at x=5 (only singular at x=0, outside this piece), so the
    // join is a genuine Step (0 vs -1/25), not Singular -- confirming the
    // fallback doesn't spuriously flag a non-singular infinite-width piece.
    std::map<std::string, std::string> potential = {
        {"[0,5)", "0"}, {"[5, inf)", "-1/x^2"}
    };
    auto joins = tise::detectPotentialStructure(potential);
    ASSERT_EQ(joins.size(), 1u);
    EXPECT_DOUBLE_EQ(joins[0].x, 5.0);
    EXPECT_EQ(joins[0].type, tise::JoinType::Step);
}

// ---------------------------------------------------------------------------
// strategicKnotsFromJoins
// ---------------------------------------------------------------------------

TEST(StrategicKnotsFromJoinsTest, MultiplicityFormulaPerOrder)
{
    for (int order : {4, 6, 8})
    {
        std::vector<tise::DetectedJoin> joins = {
            {1.0, tise::JoinType::Step}, {2.0, tise::JoinType::StitchedKink}
        };
        auto knots = tise::strategicKnotsFromJoins(joins, order);

        auto findKnotAt = [&](double x) -> const tise::StrategicKnot * {
            for (const auto &k : knots)
                if (k.x == x) return &k;
            return nullptr;
        };

        const tise::StrategicKnot *step = findKnotAt(1.0);
        ASSERT_NE(step, nullptr);
        EXPECT_EQ(step->extraMultiplicity, order - 3);

        const tise::StrategicKnot *kink = findKnotAt(2.0);
        if (order - 4 > 0)
        {
            ASSERT_NE(kink, nullptr);
            EXPECT_EQ(kink->extraMultiplicity, order - 4);
        }
        else
        {
            EXPECT_EQ(kink, nullptr);
        }
    }
}

TEST(StrategicKnotsFromJoinsTest, StitchedKinkClampsToZeroAtOrderFour)
{
    std::vector<tise::DetectedJoin> joins = {{3.0, tise::JoinType::StitchedKink}};
    auto knots = tise::strategicKnotsFromJoins(joins, 4);
    EXPECT_TRUE(knots.empty());
}

TEST(StrategicKnotsFromJoinsTest, SingularAndContinuousProduceNoKnots)
{
    std::vector<tise::DetectedJoin> joins = {
        {1.0, tise::JoinType::Singular}, {2.0, tise::JoinType::Continuous}
    };
    auto knots = tise::strategicKnotsFromJoins(joins, 8);
    EXPECT_TRUE(knots.empty());
}

// ---------------------------------------------------------------------------
// buildStrategicRadialGrid
// ---------------------------------------------------------------------------

TEST(BuildStrategicRadialGridTest, ExistingPointGetsExtraMultiplicityOnly)
{
    // nNodes=11 on [0,10] -> spacing 1.0; x=5.0 already on the grid.
    std::vector<tise::StrategicKnot> knots = {{5.0, 3}};
    auto grid = tise::buildStrategicRadialGrid(11, 0.0, 10.0, knots);
    EXPECT_EQ(grid.size(), 11u + 3u);
    EXPECT_EQ(std::count(grid.begin(), grid.end(), 5.0), 1 + 3);
}

TEST(BuildStrategicRadialGridTest, NewPointInsertionSplicesAndAddsMultiplicity)
{
    // nNodes=11 on [0,10] -> spacing 1.0; x=5.5 is NOT on the grid.
    std::vector<tise::StrategicKnot> knots = {{5.5, 2}};
    auto grid = tise::buildStrategicRadialGrid(11, 0.0, 10.0, knots);
    EXPECT_EQ(grid.size(), 11u + 1u + 2u);
    EXPECT_EQ(std::count(grid.begin(), grid.end(), 5.5), 1 + 2);
    EXPECT_TRUE(std::is_sorted(grid.begin(), grid.end()));
}

TEST(BuildStrategicRadialGridTest, RemainsNonDecreasing)
{
    std::vector<tise::StrategicKnot> knots = {{2.0, 4}, {5.5, 3}, {8.0, 1}};
    auto grid = tise::buildStrategicRadialGrid(11, 0.0, 10.0, knots);
    EXPECT_TRUE(std::is_sorted(grid.begin(), grid.end()));
}

TEST(BuildStrategicRadialGridTest, ReturnedSizeAccountsForAllInsertions)
{
    // One on-grid knot (5.0, +3) and one off-grid knot (5.5, +2): total
    // extra entries = 3 (multiplicity only) + (1 splice + 2 multiplicity).
    std::vector<tise::StrategicKnot> knots = {{5.0, 3}, {5.5, 2}};
    auto grid = tise::buildStrategicRadialGrid(11, 0.0, 10.0, knots);
    EXPECT_EQ(grid.size(), 11u + 3u + (1u + 2u));
}

// ---------------------------------------------------------------------------
// Strategic node placement -- literal "Done when" criterion
// ---------------------------------------------------------------------------

TEST(StrategicNodePlacementAccuracyTest, ImprovesOverUniformGridForBoxBarrier)
{
    // Particle in a box [0,10] with a rectangular barrier on [5,6]
    // (docs/TDSE-original-design/2026-06-28-config-yaml-schema-design.md's
    // own worked example -- the task-breakdown's literal "Done when" case).
    // Ground-state eigenvalue only; L is unused by fillBandedMatrices' own
    // potential evaluation (kept only for eigenvalueError's hydrogen-
    // specific analytic comparison, not used here), passed as 0.
    std::map<std::string, std::string> potential = {
        {"[0,5)", "0"}, {"[5,6]", "10"}, {"(6,10]", "0"}
    };
    const int order = 8;
    const int L = 0;
    const double rMin = 0.0, rMax = 10.0;

    auto groundStateEnergy = [&](const std::vector<double> &grid, int nNodesForInit) {
        bspline::BSpline bs;
        int info = bs.init(nNodesForInit, order, grid);
        EXPECT_EQ(info, 0);
        int nEn = bs.getNBSplines() - 2;
        auto [H, S] = tise::fillBandedMatrices(bs, nEn, order, L, potential);
        auto result = tise::solveGeneralizedEigenproblem(std::move(H), std::move(S), nEn, order);
        return result.values[0];
    };

    // Fine uniform reference: converged ground truth.
    const int nFine = 201;
    double eFine = groundStateEnergy(tise::buildUniformRadialGrid(nFine, rMin, rMax), nFine);

    // Coarse uniform: 21 nodes, spacing 0.5 -- x=5 and x=6 already fall
    // exactly on the grid, deliberately, so the strategic comparison below
    // adds only knot degeneracy, no new distinct spatial points -- the
    // cleanest possible reading of the task's "same node count" wording.
    const int nCoarse = 21;
    double eCoarseUniform = groundStateEnergy(tise::buildUniformRadialGrid(nCoarse, rMin, rMax), nCoarse);

    // Coarse strategic: same 21-node uniform base, degenerate knots added
    // at the two Step joins (x=5, x=6).
    auto joins = tise::detectPotentialStructure(potential);
    auto knots = tise::strategicKnotsFromJoins(joins, order);
    auto strategicGrid = tise::buildStrategicRadialGrid(nCoarse, rMin, rMax, knots);
    // Per buildStrategicRadialGrid's contract: pass the returned grid's OWN
    // size, not nCoarse -- see Gap 1 above.
    double eCoarseStrategic = groundStateEnergy(strategicGrid, static_cast<int>(strategicGrid.size()));

    double errUniform = std::abs(eCoarseUniform - eFine);
    double errStrategic = std::abs(eCoarseStrategic - eFine);
    EXPECT_LT(errStrategic, errUniform);
}
```

**Flagged risk (can't verify further without executing code):** the exact node counts (`nFine=201`, `nCoarse=21`) are chosen from physical reasoning (spacing 0.5 puts both barrier edges exactly on the coarse grid; 201 should be amply converged for a box of size 10) but not run. **First implementation step for this test specifically: build and run it before anything else.** If `errStrategic < errUniform` doesn't hold at these exact counts, the fallback is adjusting `nFine`/`nCoarse` (e.g. increasing `nFine` for better reference convergence), not abandoning the test — record whatever is found in this plan doc, same precedent as A2's `nNodes` 31→121 correction and A3's synthetic-fixture risk.

## Scope / non-goals

- No `CMakeLists.txt`, `main.cpp`, `README.md`, or `docs/SDD.md` changes.
- `docs/planning/engineer-a-plan.md` (the master doc) is left untouched, matching A1/A2/A3's precedent.
- A4b (generalizing the `B_1`/`B_N` drop convention) is explicitly out of scope — not attempted, not stubbed.
- `solveTISE`'s call site is not modified to use strategic grids by default (see "Scope: standalone, unwired" above).
- `evaluateFunction`'s existing by-value `std::map` parameter (`tise.cpp:90`, copies the whole potential map on every call) is a pre-existing signature this task doesn't touch — the new detection code calls it many times (16 samples × up to 2 sides per join, plus the value/slope probes), which is a trivial cost for the small potentials used here, not a concern worth fixing in this task's scope.
- No commits made automatically; commit timing is decided after review, matching A1/A2/A3.

## File-by-file summary

**Created:** `docs/planning/engineer-a-plan-A4.md` (this file).
**Modified:**
- `TISE/tise.hpp` — new `JoinType`/`DetectedJoin`/`StrategicKnot` + 3 function declarations (after line 175).
- `TISE/tise.cpp` — `parseInterval` refactor (in place, lines 50-88) + new A4a block: 3 file-local helpers (`isSingularApproaching`, `oneSidedValue`, `oneSidedSlope`) + 3 public functions (after line 384).
- `TISE/tests/test_tise.cpp` — 15 new `TEST`/`TEST_F` cases across 4 groups (after line 973).

## Implementation status (post-execution)

This plan has been implemented and verified. Corrections were made during implementation, recorded here per this task set's established precedent (A2's `nNodes` 31→121 correction, A3's synthetic-fixture risk check):

1. **Compile fix:** `clampedStep`'s lambda (`[]` → `[&]`) needed to capture `kFdStep` even though it's `constexpr` — passing it to `std::min` (which takes `const T&`) is an odr-use, which GCC correctly rejected without a capture. Fixed in both this doc and the actual code; the design/reasoning is otherwise unchanged.
2. **Two coverage gaps found and closed after the first green build:** `gcovr` showed `classifyJoin`'s own `Singular` return branch (an interior join where one side itself diverges, as opposed to a domain-edge singularity) and the right/last-edge singular push were both unexercised — the original 12 tests only ever hit *domain-edge* singularities (`SingularAtCoulombOrigin`, left edge only). Added `DetectPotentialStructureTest.SingularAtInteriorJoin` (`1/(1-x)` diverging at the interior join `x=1`) and `DetectPotentialStructureTest.SingularAtRightDomainEdge` (`1/(100-x)` diverging at the right domain edge `x=100`) — both hand-verified to actually diverge (sampled values checked numerically before writing the assertions) and both pass.
3. **The remaining flagged gap (`isfinite(pieceWidth)==false` clamp branch) closed on request:** added `DetectPotentialStructureTest.HandlesInteriorJoinWithInfiniteWidthPiece` — an interior join `{"[0,5)":"0"}, {"[5, inf)":"-1/x^2"}}` where the right piece's width is genuinely infinite, forcing `isSingularApproaching`'s fallback (`h0 = scale`, no width clamp) to execute. Chose a *non*-singular right side (`-1/x^2` is smooth at `x=5`, only singular at `x=0`, which is outside this piece) so the test also confirms the fallback doesn't spuriously flag an infinite-width piece as singular — the join correctly resolves to `Step` (value jump `0` vs. `-1/25`). Note: this specific branch is on a single ternary line that's always "executed" regardless of which arm runs, so it was never visible as a "missing line" in gcovr's line-coverage report (the metric this project's ≥80% target actually tracks) — it was identified by manual code reading, not by the coverage tool. Test count is 15, not 12.

Verification results: all 85 tests in `TISETests` pass (`ctest --test-dir build`), including the flagged-risk accuracy test at the proposed node counts (`nFine=201`, `nCoarse=21`) with no adjustment needed. `inInterval`'s 4 existing tests pass unchanged, confirming the `parseInterval` refactor is behavior-preserving. Coverage: `tise.cpp` line coverage is 95% (well above the 80% target); the only remaining uncovered lines in `tise.cpp` are pre-existing (the `DSBGV` failure throw, and the entirely-untested `solveTISE` wrapper, neither touched by this task) plus closing braces gcov attributes oddly — no uncovered lines remain inside any A4a code. A clean rebuild produces zero compiler warnings.

No commit made yet — paused for review, per this task set's established cadence.

## Verification

```bash
cd TISE
cmake -S . -B build -DBUILD_TESTING=ON && cmake --build build -j && ctest --test-dir build --output-on-failure
```
Confirm the literal "Done when" criterion: `StrategicNodePlacementAccuracyTest.ImprovesOverUniformGridForBoxBarrier` passes. Confirm all 4 `DetectPotentialStructureTest` cases match the hand-traced results above exactly (not just pass/fail — check the actual `JoinType` values if anything is unexpected). Confirm `inInterval`'s existing tests still pass unchanged (behavior-preserving refactor check).

```bash
cmake -S . -B build-coverage -DBUILD_TESTING=ON -DENABLE_COVERAGE=ON
cmake --build build-coverage -j
find build-coverage -name '*.gcda' -delete
ctest --test-dir build-coverage --output-on-failure
mkdir -p build-coverage/coverage
gcovr --root . --filter '.*tise\.cpp' --exclude '.*/tests/.*' --exclude '.*main\.cpp' \
      --print-summary --html --html-details -o build-coverage/coverage/index.html build-coverage
```
Confirm `tise.cpp` stays ≥80% (currently ~98.6% after A3). Confirm all four `JoinType` branches are exercised (the `DetectPotentialStructureTest` cases cover this) and both the on-grid and off-grid paths in `buildStrategicRadialGrid` are hit (the `BuildStrategicRadialGridTest` cases cover this). If gcovr reveals an uncovered branch not already listed above, add a targeted test for it then, following the same reasoning style as this plan rather than a placeholder. (The `isfinite(pieceWidth)==false` clamp branch, only reachable for an infinite-width piece, was exactly this kind of gap — closed by `HandlesInteriorJoinWithInfiniteWidthPiece`, see "Implementation status" below.)

## Notes / assumptions carried into implementation

- Interior-join detection assumes adjacent pieces' domains touch exactly (`pieces[i].upper == pieces[i+1].lower`), using `pieces[i].upper` as the join x. This mirrors `evaluateFunction`'s own existing looseness (no validation that pieces tile the domain without gaps/overlaps) — not a new assumption introduced by this task.
- `parseInterval`/`ParsedInterval` are file-local (anonymous namespace), matching the design decision to not expose numeric-bound parsing in the public header. No dedicated test file for `parseInterval` itself — it's exercised thoroughly by proxy through `inInterval`'s existing tests (unchanged) and indirectly through all 4 `DetectPotentialStructureTest` cases.
- The singularity-check clamp (Gap 3) is the one piece of this design with no precedent in either the master plan or the prior chat sketch — flagged here explicitly as new reasoning from this planning pass, not carried over from anywhere else.
- `classifyJoin`'s singularity check runs before the value/slope checks and short-circuits on the first divergent side (`||`) — a join singular on either side is `Singular` regardless of the other side's behavior, matching the task table's framing of "singular potentials" as a distinct, higher-priority category from ordinary jumps.
