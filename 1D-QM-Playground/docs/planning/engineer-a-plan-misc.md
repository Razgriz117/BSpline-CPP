# TISE Engineer-A — Miscellaneous Cleanup Tasks Implementation Plan

> This plan is scoped to the **four cleanup tasks** from `docs/planning/tise-task-breakdown.md` §4 ("Cleanup Tasks"), reproduced and elaborated as §0.4 of the master plan ([`engineer-a-plan.md`](engineer-a-plan.md) lines 62–86). Tasks A1–A5 (boundary-condition asymptote classifier, bound/continuum classification, well-containment diagnostic, strategic node placement, continuum-accuracy warning) are all complete and committed (`15db181`, `725dacc`, `9f303c5`, `021745b`, `549a62a`, `655e9c5`) — see [`engineer-a-plan-A1.md`](engineer-a-plan-A1.md) through [`engineer-a-plan-A5.md`](engineer-a-plan-A5.md). These four cleanup items are the only work `engineer-a-plan.md` still lists as outstanding. Per this task set's cadence (plan + implementation in one pass for one already-well-scoped bundle, matching A5's precedent), all four are designed and implemented together here.

## Context

`docs/planning/tise-task-breakdown.md` §4, in full (lines 190–198):

> Small, real gaps surfaced during research. Pick up whichever fits naturally alongside the adjacent workstream — none require interface/CLI work.
>
> 1. **Fix `TISE/make_and_run.sh`.** It currently invokes `./build/H-BoundStates` with zero arguments; `main.cpp:60` unconditionally dereferences `argv[1]` via `parsePiecewise(argv[1])`, so this will crash. Update the script to pass a valid JSON potential array argument matching the actual CLI contract.
> 2. **Update `TISE/README.md`.** It still describes the project as hydrogenic-only with angular momentum set via a `constexpr int L` — it doesn't mention the muparser-based general piecewise-potential mechanism that's actually implemented. Update it to describe the current JSON-argv/muparser mechanism.
> 3. **Reconcile the expression-parser decision record.** `docs/SDD.md` §11.1 and `docs/planning/resources.md` both describe the expression-parser choice (FunctionParser vs. NFParam) as "not yet chosen definitively" — but `TISE/CMakeLists.txt:16-18` has already committed to **muparser**, which isn't even one of the two shortlisted candidates. Update SDD §11.1 / `resources.md` to record muparser as the adopted choice, or open a short ADR if there's genuine appetite to revisit.
> 4. **Flag (do not fix) the unconditional time-evolution call.** `main.cpp:121-133` always calls `tevol::runTimeEvolution(...)` after solving the TISE problem, regardless of whether a bound/continuum-only run was intended — there is no `run_tise`/`run_tdse` gating, since no config-flag plumbing exists in `main.cpp` at all yet. A proper fix requires the (out-of-scope) config-driven run-gating from the interface phases. Document this as a known limitation; do not attempt to fix it in this phase.

**Assignment status:** these four items are explicitly **not** assigned to either engineer by name in the source doc ("Pick up whichever fits naturally alongside the adjacent workstream") — `engineer-a-plan.md:7` is what claims them for this workstream: *"We are picking up Engineer A's five tasks (A1–A5) ... plus the four cleanup tasks (§4) not assigned to either engineer."* None carry a `REQ-F`/`REQ-NF` ID — they're framed as "small, real gaps surfaced during research," not requirements.

**Confirmed via file mtimes: none of the four have been touched yet.** `main.cpp`, `make_and_run.sh`, and `README.md` are all still at their original (pre-A1) timestamp, and `docs/SDD.md`/`docs/planning/resources.md` are unchanged since the plan proposing their edit was written.

**Explicitly out of scope — a related but separately-sourced item is deliberately excluded:** `engineer-a-plan-A4-wiring.md` (a gap-inventory doc from a later audit, unrelated to `tise-task-breakdown.md`) separately flags a "Gap 8" — `docs/SDD.md:1437` and `docs/planning/architecture-06-20.md:427` both claim a mixed exponential+linear grid "is the scheme in the current implementation," which is false (only `buildUniformRadialGrid`/`buildStrategicRadialGrid` exist anywhere in `TISE/`; confirmed by grep — zero hits for "exponential"/"Bachau" in `tise.cpp`/`tise.hpp`). This is the same *flavor* of fix as Cleanup 3 below (stale doc vs. actual code), but it did not originate from `tise-task-breakdown.md` §4 — it surfaced from a separate, later audit, and is explicitly out of scope for "the miscellaneous tasks from the tise-task-breakdown.md doc." Left untouched here; available to pick up independently any time (the wiring doc itself notes Gap 8 has no dependency on anything else).

## Ground truth verified directly against current code

- **`main.cpp:56-60`**, confirmed by direct read: `int main(int argc, char *argv[])` has no `argc` check anywhere in the file; line 60 is `auto potential = parsePiecewise(argv[1]);`. When `argc==1`, `argv[1]` equals `argv[argc]`, which the C++ standard guarantees is a null pointer; constructing a `std::string` (via `parsePiecewise`'s `const std::string&` parameter) from a null `char*` is undefined behavior — segfaults in practice. This call is outside both of the file's two `try`/`catch` blocks.
- **`make_and_run.sh`, full current content** (2 lines): `cmake -S . -B build && cmake --build build -j6` then `cd build && ./H-BoundStates ; cd -` — zero arguments passed, directly triggering the above.
- **`L`'s actual role, verified against `tise.cpp:357-367`** (not just assumed from the master plan's prose):
  ```cpp
  bspline::D2DFun fUni = [](double, const double *) { return 1.0; };
  // Piecewise potential supplied by the caller, evaluated per-x via muparser.
  // `L` is no longer used to select the potential here; it is retained for
  // eigenvalueError()'s comparison against the analytic hydrogen spectrum.
  bspline::D2DFun fPot = [potential](double x, const double *) {
      return evaluateFunction(potential, x);
  };
  // Previous hardcoded radial hydrogen-like potential, kept for reference:
  // bspline::D2DFun fPot = [L](double x, const double *) {
  //     return radialPotential(x, L);
  // };
  ```
  The code's own comment confirms the master plan's claim directly: `L` (passed into `fillBandedMatrices` and used only by `analyticHydrogenEnergy`/`eigenvalueError`) no longer shapes the Hamiltonian's potential term at all — that comes entirely from the caller-supplied piecewise `potential` map via `evaluateFunction`. Safe to state as fact in the README.
- **`README.md`, current structure** (191 lines): ToC lists Code Structure / Changing the Angular Momentum / Dependencies / Building with CMake / Building and Running Tests / Building with `g++` / Output Files. Intro describes the potential as the fixed hydrogenic form \(V_\ell(r) = \ell(\ell+1)/2r^2 - 1/r\). Dependencies already lists muparser, but mislabeled "[Function Parser Library](https://beltoforion.de/en/muparser/)" — a name that collides confusingly with the *different*, unchosen "FunctionParser" library at `warp.povusers.org` (see Cleanup 3). No "Coverage" or "Known Limitations" section exists.
- **`docs/SDD.md` §11.1, line 1109, exact current text:** *"An expression parser for `potential.function`/`tdse.field.expression` ([§6.1], [§6.4]) — candidates from `docs/planning/resources.md` are [FunctionParser](http://warp.povusers.org/FunctionParser/) and [NFParam](https://github.com/nativeformat/NFParam); not yet chosen definitively."* Confirmed via grep: zero hits for "muparser" anywhere in `docs/SDD.md`.
- **`docs/planning/resources.md`, full current content** (5 lines, bare URLs, no annotations):
  ```
  http://warp.povusers.org/FunctionParser/
  https://github.com/nativeformat/NFParam
  https://www.boost.org/doc/libs/1_77_0/libs/math/doc/html/interpolation.html
  https://github.com/nlohmann/json
  https://github.com/jbeder/yaml-cpp
  ```
- **`main.cpp:118-133`, confirmed unconditional:**
  ```cpp
  try
  {
      tevol::runTimeEvolution(bs, er, nBSplines,
                              BS_GRMIN, BS_GRMAX,
                              TIME_STEPS, DT,
                              INITIAL_POSITION, M_OMEGA, HBAR,
                              TIMESTEPS_DIR);
  }
  catch (const std::exception &e) { ... }
  ```
  No `if` gates this. `main.cpp` has zero YAML/config parsing (its includes are `cmath, fstream, iomanip, iostream, sstream, string, map, nlohmann/json.hpp, BSpline.hpp, tise.hpp, time_evolution.hpp` — no `yaml-cpp`), confirmed by `CMakeLists.txt`: `H-BoundStates` is not linked against `yaml-cpp` (only the separate `tise_solver` stub target is). **However, `config.yaml` (project root) already defines the semantically-right flag**, currently unused by anything on this path:
  ```yaml
  run:
    run_tise:     true
    run_tdse:     false
  ```
  (`config.yaml` lines 6-8; `docs/SDD.md` line 808 documents `run_tdse` as "Whether to invoke the TDSE solver.") This makes the Cleanup-4 comment more concrete than a vague "future config flag" — there's a specific, already-named flag sitting unused.
- **`CMakeLists.txt` confirms muparser is a hard, already-committed dependency:** `pkg_check_modules(MUPARSER REQUIRED IMPORTED_TARGET muparser)`, linked into `tise_lib`.
- **Coverage impact: none.** The project's `gcovr` invocation (`docs/planning/engineer-a-plan.md` §0.2) filters to `'.*tise\.cpp'` and explicitly excludes `'.*main\.cpp'` — none of this task's changes touch `tise.cpp`/`tise.hpp`, and `main()` isn't reachable from GTest regardless, so none of it could move the coverage number either way.

## Design

### Cleanup 1 — `TISE/make_and_run.sh` + `main.cpp` argc guard

**`TISE/make_and_run.sh`, full replacement:**
```bash
cmake -S . -B build && cmake --build build -j6
cd build && ./H-BoundStates '[{"domain": "(0, 100]", "function": "-1/x"}]' ; cd -
```

**Deviation from the master plan's literal suggested argument, found and corrected during implementation:** `engineer-a-plan.md:68` originally proposed `'[{"domain": "(0, 100]", "function": "-1/x + 1/x^2"}]'` (the \(\ell=1\) hydrogenic potential). Running it against the actually-built binary showed this produces a technically-non-crashing but empty-looking demo — **"Number of Accurate Eigenvalues : 0"**, no `EigenState_*` files written — because `main.cpp`'s `constexpr int L` defaults to `0`, and `eigenvalueError` compares against the \(\ell=0\) analytic formula regardless of what potential was actually supplied (verified: \(\ell=1\) ground state \(\approx -0.125\) vs. the \(\ell=0\) analytic comparison value \(-0.5\), so the very first eigenvalue already fails `ERROR_THRESHOLD` and the printing loop `break`s immediately). Switching to `-1/x` (the \(\ell=0\) potential, matching the actual default `L`) was verified directly to print 4 accurate eigenvalues matching the closed-form hydrogen spectrum (\(-1/2, -1/8, -1/18, -1/32\)) and write 4 `EigenState_XXX` files — a genuinely working first-run demo instead of a silent no-op. This is the same class of judgment call as Case 3's smoothing-window upgrade in the A1 plan: a small, reasoned deviation from the literal source text, made because direct testing showed the literal text didn't serve its own stated purpose (a working demonstration), and documented here rather than silently substituted.

**`TISE/main.cpp`, insert immediately after line 57 (`{`), before the existing `argv[1]` comment/parse:**
```cpp
    if (argc < 2)
    {
        std::cerr << "Usage: " << argv[0] << " '<JSON array of {\"domain\":...,\"function\":...} pieces>'\n";
        return EXIT_FAILURE;
    }
```
`EXIT_FAILURE` is already used four times elsewhere in this file with no explicit `<cstdlib>` include (transitively available), so this needs no new include.

**Non-goal:** no `try`/`catch` added around malformed-JSON input (`nlohmann::json::parse` still throws uncaught if `argv[1]` isn't valid JSON) — that's a distinct failure mode from the one this task names (`argc` only), and adding it would be scope creep beyond the task's literal text.

### Cleanup 2 — `TISE/README.md`

Section-by-section edits (not a full rewrite):

1. **Table of Contents:** rename "Changing the Angular Momentum" → "Specifying the Potential"; add "Coverage" (after "Building with `g++`") and "Known Limitations" (after "Output Files") entries.
2. **Intro paragraph:** reframe from "hydrogen-like atom" / fixed potential to a general piecewise-potential solver, keeping the hydrogenic case as an example and noting the analytic-error column is only meaningful when the supplied potential actually is hydrogenic.
3. **"Code Structure" → `main.cpp` bullet:** describe argv-driven potential parsing rather than "Defines the angular-momentum parameter `L`"; add a bullet noting the Part-2 time-evolution call (cross-referencing Known Limitations).
4. **Replace "Changing the Angular Momentum" with "Specifying the Potential":** documents the `{"domain": "<interval>", "function": "<muparser expression in x>"}` DSL using `main.cpp`'s own comment example, gives the hydrogenic-reproduction example (`-1/x + 1/x^2` for \(\ell=1\)), and states that `L` now only feeds the analytic-comparison formula (per the verified code comment above) — not the potential's shape — and that mismatched `L` makes the printed error column meaningless (not validated against the potential).
5. **Dependencies:** relabel "Function Parser Library" → "muparser" with a short description, removing the FunctionParser-name collision.
6. **"Running the executable" and the manual `g++` build's run step:** both updated to pass the JSON arg (`-1/x`, per the Cleanup-1 deviation above — matches `main.cpp`'s default `L=0` so the example genuinely produces accurate eigenvalues, not just "doesn't crash"), since after Cleanup 1 running with none now prints a Usage message and exits.
7. **New "Coverage" section:** documents the existing `ENABLE_COVERAGE` CMake option + `gcovr` workflow (already established in `engineer-a-plan.md` §0.2, landed with task A1); states the ≥80% REQ-NF-010 gate is tracked on `tise.cpp`'s own reported line coverage, not a blended total.
8. **New "Known Limitations" section:** one bullet, for Cleanup 4 (exact text below).

### Cleanup 3 — reconcile the expression-parser record

**`docs/SDD.md` line 1109, replace the clause:**
> — candidates from `docs/planning/resources.md` are [FunctionParser](http://warp.povusers.org/FunctionParser/) and [NFParam](https://github.com/nativeformat/NFParam); not yet chosen definitively.

**with:**
> — **adopted: [muparser](https://beltoforion.de/en/muparser/)**, wired via `pkg_check_modules(MUPARSER REQUIRED IMPORTED_TARGET muparser)` in `TISE/CMakeLists.txt` and used by `tise::evaluateFunction` (`TISE/tise.cpp`). The candidates originally shortlisted from `docs/planning/resources.md`, [FunctionParser](http://warp.povusers.org/FunctionParser/) and [NFParam](https://github.com/nativeformat/NFParam), were superseded by this choice.

A docs update (not a new ADR) matches this repo's own precedent: all 4 existing ADRs record "defer X" decisions, whereas this is the opposite — an already-adopted, already-in-use choice with nothing left to decide.

**`docs/planning/resources.md`, insert a new line grouped with the other expression-parser candidates** (after the existing FunctionParser/NFParam lines, before Boost.Math), in the file's existing terse one-URL-per-line style:
```
https://beltoforion.de/en/muparser/  # adopted -- TISE/CMakeLists.txt, tise::evaluateFunction (see docs/SDD.md Sec. 11.1)
```

### Cleanup 4 — flag-only comment on the unconditional time-evolution call

**`TISE/main.cpp`, insert between the existing "Project Part 2" banner comment (line 120) and the `try` block (line 121):**
```cpp
    // NOTE (known limitation -- tise-task-breakdown.md Sec. 4 item 4): this call
    // is unconditional -- time evolution always runs after the TISE solve above,
    // regardless of intent. config.yaml already defines a matching run.run_tdse
    // flag (docs/SDD.md Sec. 6.1), but main.cpp has no YAML parsing of its own and
    // never reads it; gating this call on that flag needs the (out-of-scope)
    // config-driven Controller<->TISE plumbing from the interface phases. Flagged
    // here, not fixed, per the cleanup task's own scope.
```

**README "Known Limitations" bullet (added as part of Cleanup 2's edits):**
```
* **Time evolution always runs.** After solving the bound-state problem, `main.cpp` unconditionally calls `tevol::runTimeEvolution(...)` — there is currently no way to request a TISE-only run. `config.yaml` already defines a `run.run_tdse` flag for exactly this purpose (see `docs/SDD.md` §6.1), but `main.cpp` does not parse `config.yaml` (it has no YAML dependency at all) and so never reads it. Wiring this up is deferred to the Controller↔TISE configuration plumbing (`docs/planning/tise-task-breakdown.md` §4 item 4).
```

**Non-goal:** no gating logic added — matches the task's explicit "flag, do not fix" instruction.

## Scope / non-goals

- No changes to `TISE/tise.hpp`, `TISE/tise.cpp`, `TISE/tests/test_tise.cpp`, or `TISE/CMakeLists.txt` — none of the four cleanup items touch solver logic.
- No changes to `docs/planning/engineer-a-plan.md` (the master doc) — left untouched, matching A1–A5 precedent.
- No fix for Gap 8 (`docs/SDD.md:1437` / `docs/planning/architecture-06-20.md:427`'s stale exponential-linear-grid claim) — deliberately excluded, see Context above.
- No `try`/`catch` around malformed-JSON `argv[1]` input — a distinct failure mode from the one Cleanup 1 names.
- No gating logic for the time-evolution call — Cleanup 4 is explicitly comment-only.
- No new ADR for the muparser decision — a docs update fits this repo's existing precedent better (see Cleanup 3).

## File-by-file summary

**Modified:** `TISE/make_and_run.sh`, `TISE/main.cpp` (argc guard + Cleanup-4 comment), `TISE/README.md`, `docs/SDD.md` (line 1109 only), `docs/planning/resources.md`.
**Created:** `docs/planning/engineer-a-plan-misc.md` (this document).
**Not touched:** `TISE/tise.hpp`, `TISE/tise.cpp`, `TISE/tests/test_tise.cpp`, `TISE/CMakeLists.txt`, `docs/planning/engineer-a-plan.md`, `docs/planning/architecture-06-20.md`.

## Verification

```bash
cd TISE
cmake -S . -B build -DBUILD_TESTING=ON && cmake --build build -j && ctest --test-dir build --output-on-failure
```
Confirms no regressions — same test suites/counts as before this change, since nothing here touches `tise.cpp`/`tise.hpp` logic.

Manual checks (`main()` isn't reachable from GTest, so these are CLI checks, not new unit tests):
```bash
cd TISE/build
./H-BoundStates                                                          # expect: Usage message to stderr, exit 1, no crash
./H-BoundStates '[{"domain": "(0, 100]", "function": "-1/x + 1/x^2"}]'   # expect: normal solve output, no crash
cd .. && bash make_and_run.sh                                             # expect: builds, runs to completion without crashing
```

## Notes / assumptions carried forward into implementation

- No commits will be made automatically; commit timing is decided separately, after review (standing instruction for this task set).
- The "Function Parser Library" → "muparser" README relabel and the analytic-error-column caveat in "Specifying the Potential" are reasoned engineering additions slightly beyond the master plan's literal text, documented here with rationale rather than silently introduced — both directly serve the same doc-accuracy goal as the plan's own cleanup items, in sections already being touched.
- Gap 8 (stale exponential-linear-grid claim) remains available as independent future work; flagged, not fixed, here.
- This closes out every item `engineer-a-plan.md` lists as outstanding for this workstream — after this, no known undone work remains from the master plan (modulo Gap 8 and the separate, not-yet-designed A4-wiring task, both explicitly out of scope here).

## Implementation status (post-execution)

Implemented exactly as designed above, with one deviation found and corrected during verification (Cleanup 1's demo-potential swap, documented in its Design section).

1. **All 5 files edited** as planned: `TISE/make_and_run.sh`, `TISE/main.cpp` (argc guard + Cleanup-4 comment), `TISE/README.md` (full section-by-section rewrite), `docs/SDD.md` (line 1109), `docs/planning/resources.md` (one added line).
2. **Build verified clean:** `cmake -S . -B build -DBUILD_TESTING=ON && cmake --build build -j` — only `H-BoundStates` needed rebuilding (only `main.cpp` changed among compiled sources); no warnings/errors.
3. **No regressions:** `ctest --test-dir build --output-on-failure` — all 4 suites pass, 0 failures (`BSplineTests`, `UtilsTests`, `TISETests`, `TimeEvolutionTests`), exactly as before this change (nothing here touches `tise.cpp`/`tise.hpp`/`test_tise.cpp`).
4. **Crash fix confirmed directly:** `./H-BoundStates` with no args now prints the Usage message to stderr and exits with code 1 — no segfault, no crash.
5. **Normal path confirmed:** `./H-BoundStates '[{"domain": "(0, 100]", "function": "-1/x"}]'` prints 4 accurate eigenvalues (\(-4.9999999999406503\text{e-}01\), \(-1.2499999999998848\text{e-}01\), \(-5.5555555555551882\text{e-}02\), \(-3.1249999999595469\text{e-}02\) — matching \(-1/2, -1/8, -1/18, -1/32\) to the `1e-10` `ERROR_THRESHOLD`) and writes 4 `EigenState_XXX` files.
6. **`bash make_and_run.sh` run end-to-end:** builds cleanly and reproduces the same 4-eigenvalue output above with no crash.
7. **Deviation found by testing, not assumed:** the master plan's originally-proposed demo potential (`-1/x + 1/x^2`, i.e. \(\ell=1\)) was actually run against the built binary before finalizing — it printed **"Number of Accurate Eigenvalues : 0"** and wrote no eigenstate files, because `main.cpp`'s `constexpr int L` stays at its default `0` and `eigenvalueError` doesn't validate `L` against the supplied potential (exactly the caveat documented in the new "Specifying the Potential" README section). Switched the demo argument (in `make_and_run.sh`, the README's "Running the executable" section, and the manual-`g++`-build run step — but *not* the README's own worked \(\ell=1\) illustration, which intentionally keeps `-1/x + 1/x^2` to demonstrate the general-\(\ell\) case) to `-1/x`, verified per item 5 above to produce genuinely meaningful output.
8. **Diff scope:** exactly the 5 files listed in File-by-file summary; no changes to `tise.hpp`/`tise.cpp`/`test_tise.cpp`/`CMakeLists.txt`/`engineer-a-plan.md`/`architecture-06-20.md`.

No commit has been made — per standing instruction, commit timing is decided separately after review.

**Re-verified 2026-08-28, after wiring strategic node placement and A4b singular-join B-spline removal into `solveTISE`** (`docs/planning/engineer-a-plan-A4-wiring-design.md`): item 5's exact 4 eigenvalues above are still current — a fresh `./H-BoundStates '[{"domain": "(0, 100]", "function": "-1/x"}]'` run after that wiring reproduces them **byte-identical**. This was expected, not incidental: `-1/x`'s only detected structure is its Coulomb singularity at `x=0`, which coincides with the domain's own left wall, and `solveTISE` deliberately does not run `bSplinesTouchingX`-based B-spline removal for domain-edge singularities (only interior ones) — see that doc's "Baseline impact" section for why. `bash make_and_run.sh` re-run end-to-end too, same result.
