# TISE Engineer-A — Task A5 Implementation Plan (final task in the A1–A5 series)

> This plan is scoped to **task A5 only** (`E_acc` continuum-accuracy warning). Tasks A1 (asymptote classifier), A2 (bound/continuum classification), A3 (well-containment diagnostic), A4a+A4b (strategic node placement) are complete and committed (`15db181`, `725dacc`, `9f303c5`, `021745b`, `549a62a`). The five-task master plan is archived at [`engineer-a-plan.md`](engineer-a-plan.md) (§A5 sketch at lines 230–243); per-task docs for the prior four are at `engineer-a-plan-A1.md` through `engineer-a-plan-A4b.md`. This plan elaborates the master's A5 section to the same rigor, re-verified against current code (git working tree confirmed clean on branch `TISE-Generalization` before starting).

## Context

REQ-F-040 (`docs/SDD.md` §3.1) requires the solver to accept a continuum energy range and warn when it exceeds what the basis can accurately represent. The task is split into a "warning half" (A5, this task) and a "user-input half" — accepting `[E_threshold, E_max]` as real config/CLI input — which is explicitly deferred: partly to Engineer B's B2 task (energy-grid loop, which will consume the range as plain parameters, no config wiring either) and partly to the project's excluded Phase 1–3 Controller↔TISE interface work (`tise-task-breakdown.md` §6).

The physics: a B-spline basis with node spacing `Δx_node` cannot accurately represent states whose energy is high enough that their half de Broglie wavelength is comparable to or smaller than that spacing. This sets an accuracy ceiling `E_acc ≈ π²/(2·m·Δx_node²)`. A5 implements exactly two standalone functions — `computeEAcc` and `warnIfContinuumExceedsEAcc` — with no config wiring and no call site (that integration happens later, once Engineer B's B2 exists). This is the last of Engineer A's five assigned tasks; after this, only the unrelated §0.4 cleanup tasks (main.cpp argc guard, README rewrite, doc reconciliation — none touched here) remain outstanding from the master plan.

**Source docs, verified verbatim:**
- `docs/planning/tise-task-breakdown.md:103-111` — the task definition itself ("What"/"Source"/"Touches"/"Done when").
- `docs/planning/architecture-06-20.md:395-407` — the physics derivation and design decision ("Basis accuracy limit" / "Decision"), duplicated verbatim in `docs/SDD.md:1401-1417` (§12 Appendix B, "Open Design Questions" → "Continuum range" entry, ending "**Status: resolved.** Formalized as REQ-F-040.").
- `docs/SDD.md:936` (§6.4 data-validation rule) and `docs/SDD.md:1025-1030` (§8 warning taxonomy: this is a "physics warning" — computation completes, result may be unreliable; destination is `stderr`, consistently with every other C++ warning in this codebase).

## Ground truth verified directly against current code

- **Git log confirms A4a/A4b are committed** (`021745b`, `549a62a`), working tree clean, nothing staged/unstaged. (This corrects an earlier assumption from a 6-day-old memory snapshot that A4 was still just a chat sketch — it has since been fully planned and implemented.)
- **Exact current file lengths** (`wc -l`, cross-checked with direct `Read`): `TISE/tise.hpp` is **307 lines** (`solveTISE` declared at line 305, blank line 306, `} // namespace tise` at line 307 — the last line); `TISE/tise.cpp` is **755 lines** (`solveTISE`'s body closes at line 753, blank line 754, `} // namespace tise` at line 755); `TISE/tests/test_tise.cpp` is **1522 lines** (last line is the closing `}` of `SingularPotentialBSplineRemovalTest.RemovingBSplinesNearSingularityProducesFiniteSolution`). A5's code is the very last thing in each file.
- **A reusable π constant already exists**: `tise.cpp:201-204`, `namespace { constexpr Real kPi = 3.14159265358979323846; }`, internal linkage inside `namespace tise`, currently used only by `case3WindowFunction`, but in scope for anything placed later in the file — which A5 is. No `M_PI`/`<numbers>` used anywhere in `tise.hpp`/`tise.cpp` today; `kPi` is the file's own established convention, so A5's `.cpp` implementation reuses it rather than introducing a second pi constant.
- **`M_PI` in the test file is safe, verified empirically, not assumed.** `CMakeLists.txt:7` sets `CMAKE_CXX_EXTENSIONS OFF` (strict `-std=c++17`, not GNU extensions), which can matter for glibc's feature-test-gated `M_PI`. Compiled a throwaway `<cmath>`-only program with this exact flag combination on this environment's actual toolchain (GCC 11.4.0/Ubuntu 22.04) — `M_PI` resolves correctly. `test_tise.cpp` already includes `<cmath>` (line 6) and `<sstream>` (line 7); no new includes needed anywhere.
- **Test file has zero uses of `tise::Real`** — 42 existing scalar test locals all use plain `double`. New tests follow that 100%-consistent convention rather than `Real`.
- **`checkWellContainment` (A3, `tise.cpp:454-465`) confirmed as the style template** for a small standalone function: no doc comment on the `.cpp` definition (doc comment lives entirely on the `.hpp` declaration), no input validation (trusts the caller), tight body. Default arguments (e.g. `tol = 1e-3`) are declared once in the `.hpp` and **not** repeated in the `.cpp` definition — a hard C++ rule, already followed throughout this file.
- **No CMakeLists.txt changes needed.** `tise_lib`/`test_tise` are single-source-file targets (`tise.cpp`, `test_tise.cpp`); any code landing inside those files is picked up automatically. `ENABLE_COVERAGE` is already wired (landed with A1).
- **`TISE/build/` and `TISE/build-coverage/` already exist and are configured** from prior rounds — no fresh `cmake -S -B` configure needed, just rebuild.
- **Existing on-disk coverage artifact is stale — do not cite it as current.** `build-coverage/coverage/index.tise.cpp....html` reports 95.0% line coverage for `tise.cpp`, but its timestamp (2026-08-14 10:00) predates the A4b commit (2026-08-14 13:07) by ~3 hours — it's an A4a-era snapshot, not current. Treat the first fresh `gcovr` run after implementing A5 as the real baseline (see Verification).
- **Repo-wide grep confirms a clean slate**: zero hits for `EAcc`, `E_acc`, `eAcc`, `computeEAcc`, `warnIfContinuumExceedsEAcc`, or `nodeSpacing` anywhere under `TISE/` (excluding build dirs).

## Gaps found during this planning pass

Per the task's own scope framing plus direct inspection of adjacent code, the following are genuinely out of scope, unresolved, or worth flagging for whoever picks up work after A5 — not silently assumed away:

1. **`E_max`/`E_threshold` stay plain function parameters — no config/CLI schema wiring.** Explicit in the source task (`tise-task-breakdown.md:105`, "no config schema wiring, that's out of scope"). The `tise.continuum` YAML block these would eventually bind to already exists in the SDD's schema docs (`SDD.md:856-864`: `enabled`, `E_threshold`, `E_max`, `n_energies`, `n_pts`) but nothing in this task reads it.
2. **No call site — functions land standalone and unwired**, matching A1–A4's own precedent (none of those are called from `main.cpp` either). The task breakdown explicitly defers wiring to whoever integrates Engineer B's B2 (energy-grid loop) task; coordinating that call site is out of scope here.
3. **The `warnings.json` sidecar (SDD §8) is out of scope — A5 only writes to the injectable `warnOut` stream (stderr by default)**, matching the "Done when" criterion exactly. Concretely verified: `TISE/tise_solver_main.cpp:168-175` (`writeCoreOutputs`) already unconditionally writes an empty `warnings.json` (`[]`), with its own comment: *"No real physics warnings exist until a later phase... ready to be populated once warnings exist."* That file's header also states outright that it's a Phase-1 contract stub and *"nothing here should be treated as a foundation to build real numerics on top of"* — confirming there's a real (if currently inert) future destination for this warning's text, but wiring it there is not this task. Practical implication: whoever wires that later will likely lift `warnIfContinuumExceedsEAcc`'s message string into that file's `"message"` field, which is a mild argument for the message being clean, self-contained prose rather than something terser assuming a stderr-only reader — reflected in the wording chosen below.
4. **No "minimum inter-node gap" helper exists yet for non-uniform grids.** The master plan's own A5 sketch (`engineer-a-plan.md:236`) already flags that once A4 lands, the physically-correct `nodeSpacing` to pass is the grid's *minimum* inter-node gap, not a naive average. A4a/A4b are now committed, so this is no longer hypothetical — but confirmed directly, `buildStrategicRadialGrid` (`tise.hpp:255-256`) returns a plain `std::vector<Real>` with no accompanying min-spacing utility anywhere in the codebase. `computeEAcc` intentionally takes a scalar, not a grid, so this reduction stays a future call-site concern — flagged here as a concrete unimplemented follow-up rather than something A4 already solved.
5. **The source physics relation is asymptotic (≳/≲), implemented here as an exact closed-form (=).** `architecture-06-20.md`/SDD's derivation is an order-of-magnitude scaling bound throughout, but the task's "Done when" requires a specific testable value. `computeEAcc` returns the relation's leading-order coefficient exactly and `warnIfContinuumExceedsEAcc` treats it as a hard threshold for strict `>` — a deliberate reading of an inherently fuzzy physical bound as a crisp one, named explicitly (including in the `.hpp` doc comment) rather than silently glossed over.
6. **No input validation on `nodeSpacing`/`mass`** (non-positive values silently yield `inf`/`nan`/negative results rather than throwing). Deliberate, matching `checkWellContainment`'s/`classifyAsymptote`'s precedent of trusting caller-supplied physical parameters. Lower-risk than it might sound: unlike A4b's `resolveDropSet` (where a bad index risked silent out-of-bounds memory access), a bad `computeEAcc` input produces only a nonsensical-but-harmless floating-point value, not memory corruption.
7. **No canonical "atomic units" (`mass=1`) constant exists anywhere on the TISE side.** `main.cpp:23`'s `PARTICLE_MASS = 3.75` belongs to the unrelated Part-2 TDSE/wavepacket code path; the TISE solver's own formulas (`analyticHydrogenEnergy`, `radialPotential`) assume `mass=1` implicitly via hardcoded literals, never as a named constant anywhere. Not a defect A5 introduces (the task explicitly wants `mass` as a plain parameter), but whoever wires the future B2 call site needs to know to pass `1.0` — nothing in the code documents that convention today.
8. **`computeEAcc(Real nodeSpacing, Real mass)` has two adjacent same-typed parameters** — a transposition-order footgun (swapped arguments compile silently, produce a wrong-but-plausible number). Not changing the signature, since it's fixed by the master plan and by the agreed test call shape (`computeEAcc(0.5, 1.0)`); flagged as an accepted, named micro-risk rather than an invisible one.

*(Minor aside, not a numbered gap: the task breakdown's own "Source" citation, "SDD §5.2.3/§6.4/§8," is slightly loose — confirmed by direct read that §5.2.3 never actually mentions `E_acc`; the real citations are §6.4, §8, §6.1 (config schema), and Appendix B, which isn't in that citation list at all. Noted, not corrected in the source doc — matches A1–A4's existing precedent of flagging such things without editing `tise-task-breakdown.md` itself.)*

## Design

### `tise.hpp` — insert after line 305 (`solveTISE`'s declaration), before blank line 306 / closing `} // namespace tise` at line 307

```cpp
// === A5: E_acc continuum-accuracy warning (REQ-F-040, warning half) ===
// Basis accuracy ceiling (docs/planning/architecture-06-20.md "Continuum
// range" -> "Basis accuracy limit"; duplicated docs/SDD.md Appendix B
// "Continuum range" entry; see also SDD Sec. 6.4, Sec. 8). States whose
// half de Broglie wavelength is commensurate with or smaller than the
// B-spline node spacing cannot be accurately represented:
//   lambda/2 = pi/k <~ dx_node  ==>  k >~ pi/dx_node  ==>  E >~ pi^2 / (2 m dx_node^2)
// The source relation is an asymptotic ("<~"/">~") scaling bound, not an
// exact equality; this function takes its leading-order coefficient as
// the concrete E_acc threshold, per this task's closed-form "Done when"
// criterion -- treat the result as an order-of-magnitude ceiling, not a
// razor-sharp cutoff. `nodeSpacing` is a single scalar spacing; for a
// non-uniform grid (see A4's buildStrategicRadialGrid) the
// physically-correct value to pass is the *minimum* inter-node gap, not
// an average -- that reduction is a future call-site concern, not
// performed here (see docs/planning/engineer-a-plan-A5.md, Gaps).
Real computeEAcc(Real nodeSpacing, Real mass);

// Warns (to warnOut, default stderr -- SDD Sec. 8's "physics warning"
// class: computation completed, but a result may be unreliable) iff the
// requested continuum ceiling eMax exceeds the basis accuracy ceiling
// eAcc (from computeEAcc), and returns true iff it did so. Strict '>':
// eMax exactly equal to eAcc is the marginal case and does not warn,
// mirroring checkWellContainment's/classifyBoundStates' existing
// strict-inequality boundary convention. Standalone: no dependency on
// Engineer B's energy-grid loop (B2) -- wiring this at that call site is
// future integration work, matching A1-A4's "standalone, unwired"
// precedent.
bool warnIfContinuumExceedsEAcc(Real eMax, Real eAcc, std::ostream &warnOut = std::cerr);
```

### `tise.cpp` — insert after line 753 (`solveTISE`'s closing brace) and blank line 754, before closing `} // namespace tise` at line 755

```cpp
// === A5: E_acc continuum-accuracy warning (REQ-F-040, warning half) ===
Real computeEAcc(Real nodeSpacing, Real mass)
{
    return kPi * kPi / (2.0 * mass * nodeSpacing * nodeSpacing);
}

bool warnIfContinuumExceedsEAcc(Real eMax, Real eAcc, std::ostream &warnOut)
{
    const bool exceeds = eMax > eAcc;
    if (exceeds)
        warnOut << "Warning: requested continuum E_max=" << eMax
                << " exceeds the basis accuracy ceiling E_acc=" << eAcc
                << " (set by the B-spline node spacing); results at energies "
                << "above E_acc are unreliable.\n";
    return exceeds;
}
```

Warning text deliberately echoes the source doc's own phrase — *"the program should issue a warning that results are **unreliable**"* (`architecture-06-20.md:407`, duplicated `SDD.md:1407`) — and prints `E_max`/`E_acc` labels matching SDD prose and the `tise.continuum.E_max` config field name, for traceability (see Gap 3 above for why this text may get reused verbatim later).

### `tests/test_tise.cpp` — append after line 1522 (EOF)

```cpp

// ---------------------------------------------------------------------------
// computeEAcc
// ---------------------------------------------------------------------------

TEST(ComputeEAccTest, MatchesClosedForm)
{
    EXPECT_NEAR(tise::computeEAcc(0.5, 1.0), 2.0 * M_PI * M_PI, 1e-12);
}

TEST(ComputeEAccTest, ScalesInverselyWithMassAndSpacingSquared)
{
    double base = tise::computeEAcc(0.5, 1.0);

    // Doubling nodeSpacing (mass fixed): E_acc ~ 1/dx^2 -> quarters.
    EXPECT_NEAR(tise::computeEAcc(1.0, 1.0), base / 4.0, 1e-12);

    // Doubling mass (nodeSpacing fixed): E_acc ~ 1/m -> halves.
    EXPECT_NEAR(tise::computeEAcc(0.5, 2.0), base / 2.0, 1e-12);
}

// ---------------------------------------------------------------------------
// warnIfContinuumExceedsEAcc
// ---------------------------------------------------------------------------

TEST(WarnIfContinuumExceedsEAccTest, FiresWhenEMaxExceeds)
{
    std::ostringstream warn;
    EXPECT_TRUE(tise::warnIfContinuumExceedsEAcc(10.0, 5.0, warn));
    EXPECT_FALSE(warn.str().empty());
}

TEST(WarnIfContinuumExceedsEAccTest, SilentWhenEMaxBelow)
{
    std::ostringstream warn;
    EXPECT_FALSE(tise::warnIfContinuumExceedsEAcc(3.0, 5.0, warn));
    EXPECT_TRUE(warn.str().empty());
}

TEST(WarnIfContinuumExceedsEAccTest, BoundaryEqualDoesNotFire)
{
    std::ostringstream warn;
    EXPECT_FALSE(tise::warnIfContinuumExceedsEAcc(5.0, 5.0, warn));
    EXPECT_TRUE(warn.str().empty());
}

TEST(WarnIfContinuumExceedsEAccTest, WarningTextMentionsUnreliableResults)
{
    std::ostringstream warn;
    tise::warnIfContinuumExceedsEAcc(10.0, 5.0, warn);
    EXPECT_NE(warn.str().find("unreliable"), std::string::npos);
}
```

**Why 2 `ComputeEAccTest` cases are sufficient, and why a 6th `WarnIfContinuumExceedsEAccTest` was added beyond the master plan's named 5:** `MatchesClosedForm` alone can't pin the exponent on `mass`, because it calls with `mass=1.0` — any power of 1 is indistinguishable, so a buggy implementation using `mass²` or omitting `mass` entirely would still pass. `ScalesInverselyWithMassAndSpacingSquared` closes that hole by varying `mass` and `nodeSpacing` independently; together they pin the full functional form of this single straight-line arithmetic expression, leaving no further edge case. `WarningTextMentionsUnreliableResults` is added (beyond the master plan's 5) because otherwise the "echoes SDD's language for traceability" design choice is just a comment nobody enforces — mirrors A1's own precedent of splitting "did it fire" from "does the text say the right thing" into separate tests.

## Scope / non-goals

- No config/CLI wiring for `E_max`/`E_threshold` (Gap 1).
- No call site — not invoked from `main.cpp`, `tise_solver_main.cpp`, or anywhere else (Gap 2).
- No `warnings.json` sidecar population (Gap 3).
- No minimum-inter-node-gap grid utility (Gap 4) — future work for whoever wires a non-uniform-grid call site.
- No input validation on `computeEAcc`'s parameters (Gap 6).
- No changes to `CMakeLists.txt`, `main.cpp`, `README.md`, `docs/SDD.md`, or `docs/planning/resources.md` — none are needed for this task.

## File-by-file summary

**Modified:** `TISE/tise.hpp`, `TISE/tise.cpp`, `TISE/tests/test_tise.cpp` (all additions only, appended at each file's tail — no existing code touched).
**Created:** `docs/planning/engineer-a-plan-A5.md` (this document).
**Not touched:** everything else, including `CMakeLists.txt` (no changes needed — see Ground truth above).

## Verification

TDD flow, matching A3's own documented history (link-time RED, not compile-time RED):
```bash
cd TISE
# 1. Add .hpp declarations + all 6 new tests, leave .cpp unimplemented.
cmake --build build -j   # expect: link failure, "undefined reference to tise::computeEAcc(...)" /
                          # "...warnIfContinuumExceedsEAcc(...)" -- confirms failing for the right reason.
# 2. Add the .cpp implementation.
cmake --build build -j && ctest --test-dir build --output-on-failure   # expect: all tests green, including the 6 new ones.

# 3. Fresh coverage measurement (do not reuse the stale on-disk artifact):
cmake -S . -B build-coverage -DBUILD_TESTING=ON -DENABLE_COVERAGE=ON
cmake --build build-coverage -j
find build-coverage -name '*.gcda' -delete
ctest --test-dir build-coverage --output-on-failure
gcovr --root . --filter '.*tise\.cpp' --exclude '.*/tests/.*' --exclude '.*main\.cpp' \
      --print-summary --html --html-details -o build-coverage/coverage/index.html build-coverage
# Confirm tise.cpp stays >= 80% (REQ-NF-010); both new functions should show as fully covered
# given all 4 logical branches (exceeds/doesn't, plus computeEAcc's single straight-line path) are exercised.
```

## Notes / assumptions carried into implementation

- No commits will be made automatically — commit timing is decided separately, after review (standing instruction for this task set).
- The existing `build-coverage/` coverage artifact predates A4b and must not be cited as a current baseline (Ground truth above) — the post-implementation section below will record a freshly measured number.
- Warning-message wording and the 6th test are reasoned engineering additions beyond the master plan's literal text, documented above with rationale rather than silently introduced.
- This is the last of Engineer A's five assigned tasks. The master plan's separate §0.4 cleanup items (an unrelated `main.cpp` argc crash, README rewrite, SDD/resources.md doc reconciliation) remain outstanding but are out of scope for "task A5" as requested — noted here only so they aren't mistaken for part of this task's gaps.

## Implementation status (post-execution)

Implemented exactly as designed above — no corrections needed, the design's line numbers, code, and test bodies were used verbatim. TDD flow followed and verified at each step, matching A3's link-time-RED precedent:

1. **RED**: added the `.hpp` declarations and all 6 tests, left `.cpp` unimplemented, built. Got a link failure with `undefined reference to tise::computeEAcc(double, double)` and `undefined reference to tise::warnIfContinuumExceedsEAcc(double, double, std::ostream&)` — both new symbols, referenced from exactly the new tests, no other errors. Confirmed failing for the right reason.
2. **GREEN**: added the `.cpp` implementation, rebuilt clean. `ctest --test-dir build` — all 4 suites pass (`BSplineTests`, `UtilsTests`, `TISETests`, `TimeEvolutionTests`), 0 failures. Ran the 6 new tests individually via `--gtest_filter="ComputeEAccTest.*:WarnIfContinuumExceedsEAccTest.*"` for an explicit per-test confirmation, not just an aggregate pass count — all 6 `OK`. Full `TISETests` suite is now at 106 cases (100 pre-existing + 6 new).
3. **Coverage** (fresh `gcovr` run, not the stale on-disk artifact flagged in Ground truth above): **`tise.cpp` at 96% line coverage (384/400 executable lines)** — the master plan's actual REQ-NF-010 gate ("`tise.cpp`'s own reported percentage"), comfortably clearing the ≥80% bar. (The blended tise.cpp+test_tise.cpp figure `gcovr --print-summary` reports is 98.9%/1489 of 1505 — that number conflates the test file's trivially-100%-covered own lines and is not the gate this project uses; 96% is the correct one to cite.) Confirmed via the line-level `--txt` report that **both new functions are fully covered — neither appears anywhere in the "Missing" line list**. The 16 lines gcovr does report missing in `tise.cpp` (`434,452,632,689,738,740,742-746,748-749,751-753`) are all at or below the pre-A5 line 753, i.e. entirely within/before `solveTISE`'s body (which no test calls directly — tests exercise `fillBandedMatrices`/`solveGeneralizedEigenproblem` directly instead) plus four other scattered pre-existing lines — none of this is new or touched by A5; coverage was not regressed.
4. **Diff scope confirmed minimal and exactly as planned**: `git diff --stat` shows only `tise.hpp` (+30), `tise.cpp` (+17), `tests/test_tise.cpp` (+52) — 99 insertions, **0 deletions**, no existing code touched, nothing outside the three planned files. Working tree otherwise clean; this plan doc is the only new file.

No commit has been made — per standing instruction, commit timing is decided separately after review.
