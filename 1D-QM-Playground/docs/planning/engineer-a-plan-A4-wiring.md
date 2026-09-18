# TISE Engineer-A — A4/A4b Wiring Gaps (findings only, not yet a design)

> Not an implementation plan yet — this is a gap inventory produced by an audit (2026-08-21) checking the current codebase against `docs/planning/architecture-06-20.md`'s "Strategic node placement" decision and REQ-F-050. A4a (`021745b`) and A4b (`549a62a`) are complete, tested, and committed — see [`engineer-a-plan-A4.md`](engineer-a-plan-A4.md) / [`engineer-a-plan-A4b.md`](engineer-a-plan-A4b.md). What follows is everything standing between that standalone library capability and REQ-F-050's "as the default collocation scheme" actually being true of a real solve. Per this task set's established cadence ([[feedback-plan-storage-and-cadence]]), promote this into a full A1-A5-style plan doc (Context/Design/Tests/Verification) when this is picked up as its own task — do not implement straight from this list.

## Context

`architecture-06-20.md` (§"Strategic node placement", stakeholder decision) and `docs/SDD.md` REQ-F-050 both specify: *"Use a uniform B-spline basis as the default, augmented only by the strategic nodes dictated by the potential structure."* A4a/A4b implement the detection and construction machinery this requires:

- `detectPotentialStructure` / `strategicKnotsFromJoins` / `buildStrategicRadialGrid` (`tise.hpp:226-256`, `tise.cpp`) — build a strategic (non-uniform) grid from a piecewise potential.
- `bSplinesTouchingX` (`tise.hpp:272`) plus the generalized `dropSet` parameter on `fillBandedMatrices`/`eigenstateCoefficients` (`tise.hpp:141-143`, `289-293`) — remove B-splines near a singular point.

All of it is unit-tested (98 tests total in `TISETests` as of A4b) and includes end-to-end accuracy/correctness demonstrations. **None of it is called from `solveTISE` or `main.cpp`.** The gaps below are what "wiring it in" actually has to resolve — several are not simple call-site edits, because the current API doesn't have a place to put the extra information a caller would need.

## Gaps

**Gap 1 — `solveTISE` never calls the A4a grid-building pipeline.**
`tise.cpp:740`, inside `solveTISE`, calls `buildUniformRadialGrid(nNodes, rMin, rMax)` unconditionally — never `detectPotentialStructure`/`strategicKnotsFromJoins`/`buildStrategicRadialGrid`, even though `solveTISE` already has the `potential` map and `order` it would need to run that pipeline.

**Gap 2 — `solveTISE` never calls the A4b removal pipeline.**
`tise.cpp:751` calls `fillBandedMatrices(bs, nEn, order, L, potential)` with no `dropSet` argument, so it always gets the classic `{1, nBSplines}` default (`resolveDropSet`'s fallback). A `Singular` join detected by `detectPotentialStructure` is never turned into a `bSplinesTouchingX` call and fed back in — the table's 4th row (remove B-splines at a singular point) stays detection-only in production, exactly the gap A4b's own scope note flagged as *"needed for the table's 4th row to be more than detection-only"* (`engineer-a-plan.md:189`).

**Gap 3 — `solveTISE` doesn't expose the grid it actually used, and `main.cpp` independently reconstructs one.** This is the structural blocker, not just a missing call.
`solveTISE`'s signature (`tise.hpp:305`) returns only `EigenResult` — eigenvalues/eigenvectors/dim, no grid, no `nBSplines`, no `dropSet`. Meanwhile `main.cpp:82-88` calls `buildUniformRadialGrid` **again**, independently, and builds its own `BSpline` object from it, to get something to hand to `writeEigenstate` (`main.cpp:114`) and `runTimeEvolution` (`main.cpp:123-127`). Today this duplication is harmless because both builds are the same uniform grid. If `solveTISE` starts building a strategic (non-uniform) grid internally without `main.cpp` changing, `main.cpp`'s independently-rebuilt `BSpline` would silently diverge from the one `solveTISE` actually solved against — eigenvector coefficients would then be decoded against the wrong basis. Not a crash; silently wrong wavefunctions/energies downstream. Wiring Gap 1 safely requires giving callers a way to get back the grid (and, per Gap 4, the drop-set) `solveTISE` used — a signature/return-type change, not a one-line edit.

**Gap 4 — `main.cpp:113`'s `eigenstateCoefficients` call never passes a `dropSet`.**
Same failure mode as Gap 3, specifically for A4b: if `solveTISE` starts using a non-default `dropSet` internally (Gap 2), `main.cpp:113`'s `tise::eigenstateCoefficients(er.vectors, iEn, nEn, nBSplines)` call — no `dropSet` argument, so it defaults to `{1, nBSplines}` — would misattribute coefficients to the wrong physical B-spline indices whenever the real drop-set differs from the classic one.

**Gap 5 — no non-uniform-aware `nodeSpacing` for A5's `computeEAcc`.**
`tise.hpp:318-322`'s own header comment already flags this: `computeEAcc(nodeSpacing, mass)` takes a single scalar spacing, but "for a non-uniform grid (see A4's `buildStrategicRadialGrid`) the physically-correct value to pass is the *minimum* inter-node gap, not an average." No helper exists yet that reduces a `std::vector<Real>` grid (with its degenerate/repeated knots) down to that minimum gap. Irrelevant while `solveTISE` only builds uniform grids; becomes a real blocker the moment Gap 1 lands, since A5 is also unwired into `solveTISE` today regardless (see [[project-engineer-a-status]]).

**Gap 6 — no config/CLI surface to reason about.**
`main.cpp`'s only runtime input is the potential JSON (`parsePiecewise`, `main.cpp:44-54`); grid parameters (`BS_NNODS`, `BS_ORDER`, `BS_GRMIN`, `BS_GRMAX`, `main.cpp:29-32`) are compile-time constants. There's no flag distinguishing "use strategic placement" from "use plain uniform." The architecture doc's own framing (*"these strategic knots are determined automatically by the program... not something the user needs to set manually"*) suggests wiring should make it unconditional rather than opt-in, but that's a design call this inventory doesn't resolve — flagging it so it isn't silently assumed either way.

**Gap 7 — REQ-F-050 and ADR-0002 both describe strategic placement as already the default.**
`docs/SDD.md`'s REQ-F-050 text and `docs/adr/0002-defer-wkb-collocation.md` (*"strategic placement... always applied automatically regardless of this decision"*) both read as present-tense claims about shipped behavior. Neither is true of the current solve path (Gaps 1-2). `engineer-a-plan-A4.md:391-392` already flagged this as "aspirational for the eventual wired system" at A4a's own landing — this gap is just making the same point explicit in one place scoped to the wiring task specifically, so whoever picks this up doesn't have to re-derive it.

**Gap 8 — adjacent, independent doc defect found in the same audit (not itself a wiring gap).**
`architecture-06-20.md:427` and `docs/SDD.md:1437` both claim the mixed exponential+linear (Bachau Appendix A.1) grid "is the scheme in the current implementation." No exponential-linear grid builder exists anywhere in `TISE/` — only `buildUniformRadialGrid` is wired in anywhere. Unrelated to strategic placement specifically, but in the same neighborhood (grid construction) and worth fixing whenever these docs next get touched.

## Non-goals of this doc

- No code, tests, or insertion points proposed here — that's the next plan doc's job.
- Doesn't decide Gap 6's design question (unconditional vs. config-gated) or exactly how `solveTISE`'s return type/signature should change to close Gap 3/4 (e.g., a new struct bundling `EigenResult` + grid + `nBSplines` + resolved `dropSet`, vs. output parameters, vs. a new `SolveContext` the caller builds first and passes in).

## Next step

When this becomes its own implementation task, write `engineer-a-plan-A4-wiring.md`'s full design pass the way A1-A5 did: resolve Gaps 3/4/6's open questions first (they gate everything else), then wire Gap 1 → Gap 2 → Gap 5, updating `main.cpp` and `solveTISE`'s tests together so the "same grid used throughout" invariant has a regression test, not just documentation. Gap 8 can be fixed independently, any time, with no dependency on the rest.
