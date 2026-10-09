# PR #18 ("Phy5606 projects") Review

**Reviewed:** 2026-10-02
**PR:** [Razgriz117/BSpline-CPP#18](https://github.com/Razgriz117/BSpline-CPP/pull/18), head `PHY5606_Projects` @ `433e11a`, base `main` @ `73229c8`
**Scope:** 6 commits, 29 files changed, +2820 / -197 lines. Not yet merged as of this review.

## 1. Summary / Verdict

**Do not merge as-is.** Two independently-reproduced correctness regressions sit in the solver's physics core, not its edges: the generalized eigensolver that replaced the old LAPACK `DSBGV` call does not reliably detect a degenerate overlap matrix the way `DSBGV` did, and the new delta-potential feature has no guard against the one input shape (colliding/duplicate delta positions) that produces exactly that degeneracy. Four smaller, also-reproduced bugs sit around them: an unguarded division that turns `n_pts_eigenstate: 1` into a silently-all-`NaN` "successful" run, two Python crashes in the new plotting code reachable from ordinary config typos, and a validation gap where a config that passes the Python pre-flight check is rejected by the C++ solver anyway. None of the six have test coverage, so CI would not catch any of them today.

This PR is the implementation of the delta-potential feature that `docs/adr/0014-defer-delta-potential-join-detection.md` (on `main`, pre-PR) marked deferred — it adds `potential_deltas` as a config section (not, as ADR-0014's deferred path assumed, as a `JoinType` detected from `potential` expressions), and it swaps the solver's eigenvalue backend from LAPACK's banded `dsbgv_` to `Eigen::GeneralizedSelfAdjointEigenSolver`. Both changes are real engineering improvements — `dsbgv_`'s hardcoded `-I/usr/include/eigen3` dependency was Debian-only, and the delta feature is physically sound where it has been tested (a single delta reproduces the exact textbook energy to 10 significant figures, verified below) — but the error-handling contract from the code they replaced was not fully carried over.

The rest of the PR (new `eigenstates.dat`/`continuum_states.dat`/`potential.dat` output tables, the `plot_spectrum_overview` figure, build-portability fixes, the student guide) is sound in its happy path but duplicates a fair amount of code across the table writers and the CMake dependency blocks, which is where two of the six bugs above actually live.

## 2. Methodology

The PR was reviewed in two passes.

**Pass 1 — broad coverage.** Seven independent review passes, each briefed on a distinct angle, read the full diff against the actual PR worktree (checked out as a local `pr-18` branch) and reported findings back for consolidation:

1. A **removed-behavior audit** — diffing every error-handling path that existed before the PR against what replaced it, specifically hunting for a guard that used to fire and silently doesn't anymore.
2. A **language-pitfall sweep** — C++ signed/unsigned, `thread_local`/`static` lifetime, dangling-reference/use-after-move, float-equality; Python mutable-default/late-binding-closure/truthiness; CMake cache/`FORCE`/`find_package` ordering.
3. A **cross-file consistency tracer** — for each new parameter/field threaded through multiple functions, checking every call site agrees on its meaning.
4. A **line-by-line diff scan** with its own build of the solver, run against hand-written configs to empirically trigger the hypothesized failure modes rather than just reasoning about them.
5. A **root-cause/architecture scan** — symptom-vs-root-cause framing for duplication and signature growth, cross-checked against the project's own ADRs.
6. Three **cleanup-focused scans** (simplification, reuse/duplication, efficiency) covering the quality side separate from correctness.

**Pass 2 — independent re-verification (this document's author).** Rather than transcribing Pass 1's claims, every finding that ended up in Section 3 below was re-derived first-hand:

- The PR worktree's `TISE/` CMake project was configured and built from source (`cmake -S TISE -B TISE-build -DTISE_FETCH_DEPENDENCIES=OFF -DCMAKE_BUILD_TYPE=Release && cmake --build TISE-build --target tise_solver`, system Eigen/yaml-cpp/muparser — no network fetch needed), producing a real `tise_solver` binary.
- Minimal YAML configs (reproduced verbatim in Section 5) were run against that binary to trigger each C++-side bug directly, rather than relying on a report of having done so.
- A standalone ~40-line Eigen3 program was written and compiled against the system's `/usr/include/eigen3` headers to characterize `Eigen::GeneralizedSelfAdjointEigenSolver`'s error-detection behavior in isolation, across six matrix pathologies. **This independent check corrected Finding 3.1 below**: an initial Pass-1 report claimed an exactly-singular overlap matrix passes the solver's check silently; direct testing showed that specific case is actually caught, while several more realistic near-degenerate cases are not. The finding below reflects the corrected, independently-measured version, not the original claim.
- The two Python findings (Section 3.3, 3.4) were reproduced in isolated interpreter sessions against the actual functions imported from the PR's `analysis.py`/`controller.py`, not synthesized from memory of what the code does.
- The Eigen-triangle and stale-comment findings (Section 4) were checked against the installed Eigen headers and the PR's actual comment text via `grep`, not asserted from the diff alone.

Every command and its literal output used to reach a conclusion below is reproduced in Section 5 so the results can be checked independently.

## 3. Correctness Findings

### 3.1 The new eigensolver's error check does not reliably detect a degenerate overlap matrix

- **File:** `1D-QM-Playground/TISE/tise.cpp:1359-1365` (`solveGeneralizedEigenproblem`)
- **Status:** Open, not fixed

**Symptom.** The old code called LAPACK's `dsbgv_` and threw `"DSBGV failed with info=..."` when the overlap matrix `S` was not positive-definite (confirmed absent from this PR's `main`-branch base; see Section 5.1). The new code instead checks:

```cpp
Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>
    solver(A, B, Eigen::ComputeEigenvectors | Eigen::Ax_lBx);

if (solver.info() != Eigen::Success)
    throw std::runtime_error(
        "Eigen GeneralizedSelfAdjointEigenSolver failed; the overlap matrix "
        "is most likely not positive definite");
```

The comment immediately above this code (`tise.cpp:1256-1258`) explicitly asserts this preserves the old guarantee: *"both properties the continuum construction downstream relies on, and both matching what the previous LAPACK implementation guaranteed."* That assertion is false for a meaningful subset of degenerate inputs.

**Root cause.** `GeneralizedSelfAdjointEigenSolver::compute()` internally runs `Eigen::LLT<MatrixXd> cholB(matB)` to reduce the problem to a standard eigenproblem, but it never inspects `cholB.info()` — the one place a non-positive-definite `B` would actually be flagged. `solver.info()` instead reports only whether the *separate* tridiagonalized QR iteration (on the already-transformed standard eigenproblem) converged — a condition essentially orthogonal to whether `B` was positive-definite in the first place.

**Independent verification.** A standalone program (`eigen_repro_test2.cpp`, Section 5.2) was compiled against `/usr/include/eigen3` and run against six synthetic `(A, B)` pairs with `A = I`:

| Case | `B` | `solver.info()==Success`? | Returned eigenvalues |
|---|---|---|---|
| 1. Indefinite | `diag(1, -1, 1)` | **true** | `1 1 1` — wrong |
| 2. Exactly singular | `diag(1, 0, 1)` | false (caught) | `-nan -nan -nan` |
| 3. Near-singular, still PD | `diag(1, 1e-14, 1)` | **true** | `1 1 1e+14` — spurious |
| 4. Tiny negative (roundoff-scale) | `diag(1, -1e-10, 1)` | **true** | `1 1 1e+20` — spurious |
| 5. Rank-deficient (duplicate row) | `[[1,1,0],[1,1,0],[0,0,1]]` | **true** | `0.38 1 2.62` — plausible-looking, still wrong |
| 6. Same, with fp-scale noise | as above + `1e-16` noise | **true** | `0.38 1 2.62` — plausible-looking, still wrong |

A direct `Eigen::LLT` Cholesky factorization of the same indefinite matrix from case 1 correctly reports `info() == NumericalIssue`, confirming the information needed is computed somewhere in Eigen's call graph — `GeneralizedSelfAdjointEigenSolver` simply never surfaces it.

The most consequential row is **5/6**: a rank-deficient `B` arising from two linearly dependent basis functions (physically, the situation Finding 3.2 below can create) passes the check and returns finite, unremarkable-looking numbers — not `NaN`, not a crash, not even an unusually large value. There is nothing a caller or a human skimming output could use to notice anything went wrong.

**Why this matters now, specifically.** This PR's own `eigenstates.dat`/`eigenvalues.dat` are, per the PR's own `TISE/README.md` update, "the ones to reach for if you are writing your own analysis." There is no `isfinite`/sanity check anywhere downstream in `tise.cpp`/`tise_solver_main.cpp` on the solve results before they are written to those files.

**Recommended fix.** After `GeneralizedSelfAdjointEigenSolver::compute()`, explicitly check `Eigen::LLT<MatrixXd>(B).info() != Eigen::Success` (or reuse the Cholesky Eigen's solver already computes internally, if a non-breaking way to access it exists in the version in use) before trusting `solver.info()`. This restores the actual guarantee the comment at `tise.cpp:1256` already claims is in place. A supplementary `isfinite` check on `result.values`/`result.vectors` before any file write would be a cheap second line of defense against whatever this specific fix misses.

### 3.2 Colliding delta-potential positions degenerate the basis with no validation, and the resulting crash depends on exactly this reliability gap

- **File:** `1D-QM-Playground/TISE/tise.cpp:1929-1952` (`buildStrategicGridAndDropSet`) and `~1610-1619` (`buildStrategicRadialGrid`'s knot-insertion logic)
- **Status:** Open, not fixed

**Symptom.** Two `potential_deltas` entries at the same (or very close) `x`, or a delta landing exactly on a potential's own `Step`/kink join, crash `tise_solver`.

**Root cause.** Each interior `Singular` join and each `DeltaTerm` independently pushes `order - 2` extra knot copies at its own `x` into the same `knots` vector (`tise.cpp:1938`, `:1952`). `buildStrategicRadialGrid`'s insertion logic sums multiplicities from distinct sources that land at the same `x` (within `1e-12`) rather than capping the total at what a single degenerate point actually needs. `validateDeltaTerms` (`tise.cpp:1887`, called at `:1949`) only checks that each delta's `x` lies strictly inside `(rMin, rMax)` — it never checks for a duplicate delta location or a collision with a detected potential join. The same gap exists on the Python side: `controller.py`'s `validate_potential_deltas` (around line 202) also only checks domain bounds.

Once a knot's total multiplicity reaches `order`, the corresponding B-spline basis function collapses to zero-measure support, which is exactly the rank-deficiency shape exercised in Finding 3.1's cases 5/6 above.

**Independent verification.** Built from source and run directly (full commands and output in Section 5.3):

```
$ ./tise_solver --config duplicate-delta.yaml --output-dir out-dup
tise_solver: Eigen GeneralizedSelfAdjointEigenSolver failed; the overlap matrix is most likely not positive definite
exit=1
```

against a config whose only change from a working control run (`order: 8`, `potential_deltas: [{x: 0.0, strength: -1.0}]`, which independently reproduces the exact analytic hydrogen-delta-well energy `E_0 = -0.5` to `-5.0000000000002454e-01`) is a second, duplicate `{x: 0.0, strength: -1.0}` entry.

In this specific repro, the resulting degeneracy happened to be clean enough to land in Finding 3.1's "case 2" bucket (caught, with a — misleadingly worded, since the real cause is the duplicate knot, not an inherently bad overlap matrix — error). Per Finding 3.1, a less exact collision (e.g. two deltas at very close but not bit-identical `x`, or a delta a few floating-point units away from an existing join) is not guaranteed to land in that bucket, and could instead silently return the finite-looking wrong numbers of cases 5/6. Whether a given misconfiguration crashes or silently corrupts output is accidental, not a designed property of the code.

**Recommended fix.** In `validateDeltaTerms` (and its Python-side counterpart, `validate_potential_deltas`), reject any delta whose `x` coincides (within some explicit tolerance, consistent with the `1e-12` the knot-insertion logic already uses) with another delta's `x` or with a detected potential join location. This is a cheap, purely additive validation change — it does not touch the solve path at all.

### 3.3 `n_pts_eigenstate`/`n_pts == 1` silently produces an all-`NaN` output file reported as a successful run

- **Files:** `1D-QM-Playground/TISE/tise.cpp:1751` (`writeEigenstateTable`), `:1808` (`writeContinuumTable`), `:1829` (`writePotential`) — all three share the formula below
- **Status:** Open, not fixed

**Symptom.** A config with `tise.n_pts_eigenstate: 1` (or `tise.continuum.n_pts: 1`) exits 0 and writes a grid-point file whose every value is `-nan`.

**Root cause.** All three new table writers compute each grid point as:

```cpp
x = rMin + (rMax - rMin) * static_cast<Real>(ix - 1) / static_cast<Real>(npts - 1)
```

for `ix` in `[1, npts]`. When `npts == 1`, the denominator `npts - 1` is `0`, making every `x` (and therefore every dependent value) `0.0 / 0.0 = NaN`. Neither `tise_solver_main.cpp` (reads the field with no range check) nor `controller.py`'s config validator rejects `npts <= 1`.

**Independent verification.** Reproduced directly against the built binary (full output in Section 5.4):

```
$ ./tise_solver --config npts1.yaml --output-dir out-npts1
tise_solver: warning: 1 of 91 computed states are below E=0.0 (...)
exit=0

$ tail -1 out-npts1/eigenstates.dat
-nan   -nan   -nan   -nan   -nan   -nan   -nan   -nan   ...  (all 91 columns)
```

The only diagnostic the run emits is an unrelated informational warning about bound-state count — nothing flags the actual problem, and the process exit code is `0`.

**Recommended fix.** Validate `n_pts_eigenstate >= 2` and `tise.continuum.n_pts >= 2` wherever the other `tise.*` numeric fields are already validated in `controller.py`, and add a defensive check in `tise_solver_main.cpp` (or centrally, per Section 4.2's recommended shared grid helper) so a config bypassing the Python layer still fails cleanly rather than writing corrupt output.

### 3.4 `plot_spectrum_overview` crashes on an uncaught `ValueError` for an empty eigenstate

- **File:** `1D-QM-Playground/analysis.py:600`
- **Status:** Open, not fixed

**Symptom.** If any plotted eigenstate's sample-point list is empty, `analysis.py`'s `main()` surfaces a raw Python traceback instead of the clean `AnalysisError`-wrapped message every other failure mode gets.

**Root cause.**

```python
peak = max((max(abs(p.phi) for p in st) for _, _, st in drawn), default=1.0) or 1.0
```

The *outer* `max()` has `default=1.0`, guarding the case where `drawn` itself is empty (already handled by an earlier `if not drawn: return`). The *inner* `max(abs(p.phi) for p in st)` has no such default — Python's `max()` raises `ValueError: max() iterable argument is empty` on an empty iterable with no `default` given. Finding 3.3 is one concrete way to reach this: `n_pts_eigenstate <= 0` produces a zero-row `eigenstate_NNN.dat`, which parses to an empty `st`.

**Independent verification.**

```
$ python3 -c "
drawn = [(0, -0.5, [])]
peak = max((max(abs(p.phi) for p in st) for _, _, st in drawn), default=1.0) or 1.0
"
ValueError: max() iterable argument is empty
```

reproduced in isolation against the exact expression from the PR's `analysis.py:600`. The resulting `ValueError` is not a subclass of `AnalysisError`, so `main()`'s `except AnalysisError` clause does not catch it.

**Recommended fix.** Add `default=0.0` (or skip empty-`st` entries from `drawn` upstream, which is arguably the more correct fix — an eigenstate with zero sample points shouldn't be plotted at all) to the inner `max()` call.

### 3.5 `potential_deltas: null` passes Python pre-flight validation but is rejected by the C++ solver

- **Files:** `1D-QM-Playground/controller.py:202` (`validate_potential_deltas`) vs. `1D-QM-Playground/TISE/tise_solver_main.cpp:134-139` (`parseDeltaConfig`)
- **Status:** Open, not fixed

**Symptom.** A config with the key present but empty (`potential_deltas:` with nothing after the colon, or explicit `potential_deltas: null`) passes `controller.py`'s config validation — whose own docstring states its purpose is *"catching them here gives the error before a subprocess is launched"* — and then fails when the `tise_solver` subprocess actually runs.

**Root cause.** PyYAML parses both an *absent* `potential_deltas` key and a key present with an explicit `null` value to the same Python `None`. `validate_potential_deltas` treats any `None` as "no deltas" and returns cleanly (`controller.py:202-203`). `tise_solver_main.cpp`'s `parseDeltaConfig`, using yaml-cpp, distinguishes them: a key with an explicit null value is a *defined* node (`IsDefined() == true`), so the function's `if (!deltasNode) return deltas;` early-out does not fire; execution falls through to `if (!deltasNode.IsSequence())`, which is true for a null node, and throws.

**Independent verification.** Both sides reproduced directly:

```
$ python3 -c "
import controller
controller.validate_potential_deltas(None, (-20.0, 20.0))
print('validation PASSED')
"
validation PASSED

$ ./tise_solver --config null-delta.yaml --output-dir out-null   # potential_deltas: <empty>
tise_solver: 'potential_deltas' must be a YAML list of {x: ..., strength: ...} mappings
exit=1
```

**Recommended fix.** Make `validate_potential_deltas`'s notion of "absent" match yaml-cpp's: either reject an explicit null the same way the C++ side does (closing the gap by tightening Python), or have `parseDeltaConfig` treat a defined-but-null node the same as an absent one (closing it by loosening C++). Either is a one-line fix; the important part is making the two sides agree, since the whole point of the Python-side check is to make the C++ side's behavior predictable in advance.

### 3.6 `visualization: null` crashes the new `spectrum_overview` toggle read (pre-existing pattern, newly most exposed)

- **File:** `1D-QM-Playground/analysis.py:676`
- **Status:** Open, not fixed. Lower confidence/priority than 3.1-3.5 — this pattern pre-dates the PR.

**Symptom.** A config with `visualization:` present but with no value (parses to `None`) crashes with `AttributeError` instead of a clean error.

**Root cause.** `cfg.get("visualization", {}).get("spectrum_overview", True)` — `dict.get`'s default substitutes only when the key is *absent*. When the key is present with value `None` (which is exactly how YAML parses `visualization:` with nothing after it), `.get("visualization", {})` correctly returns `{}`'s *non-use* — it returns `None`, because the key *was* found — and the chained `.get(...)` on `None` raises.

```
$ python3 -c "
import yaml
cfg = yaml.safe_load('visualization:\n')
cfg.get('visualization', {}).get('spectrum_overview', True)
"
AttributeError: 'NoneType' object has no attribute 'get'
```

This exact pattern already existed pre-PR for the `eigenstates`/`bound_states_squared`/`phase_shifts` toggles (confirmed via `git show <base>:analysis.py`) — it is not newly introduced. It is flagged here because the new `spectrum_overview` line runs *first* and, unlike its siblings, defaults to `True` (opt-out rather than opt-in), making it the earliest and most easily triggered occurrence of a pre-existing bug.

**Recommended fix.** Replace `cfg.get("visualization", {})` with `(cfg.get("visualization") or {})` at all four call sites (the three pre-existing plus the new one), or resolve `visualization_cfg = cfg.get("visualization") or {}` once and reuse it.

## 4. Code-Quality Findings (not correctness bugs)

These do not block merge on their own but are worth fixing alongside Section 3, since two of them (4.2, 4.3) directly border the bugs above.

### 4.1 The matrix-fill loop writes a triangle Eigen never reads

`tise.cpp:1346-1349` writes both `A(i,j)`/`A(j,i)` and `B(i,j)`/`B(j,i)` for every band entry. Checked against the installed Eigen3 headers:

```
$ grep -n "selfadjointView" /usr/include/eigen3/Eigen/src/Eigenvalues/GeneralizedSelfAdjointEigenSolver.h
184:    MatrixType matC = matA.template selfadjointView<Lower>();
```

`GeneralizedSelfAdjointEigenSolver` unconditionally reads only the lower triangle — there is no `Upper`/`Lower` option in its constructor API at all. The `A(i,j) = H[k]` / `B(i,j) = S[k]` (upper-triangle) writes are therefore dead code in their entirety, not merely redundant on the diagonal. This could mislead a future maintainer into thinking a full symmetric fill is required. Simplify to writing only the lower triangle.

### 4.2 Stale comments describe the removed LAPACK call

`tise.hpp:540` ("Solve H c = E S c via LAPACK DSBGV.") and `tise_solver_main.cpp:469-471` ("solveGeneralizedEigenproblem's internal LAPACK call overwrites its own copies, not these...") both still describe `dsbgv_`'s in-place-overwrite semantics. Confirmed via `grep -n DSBGV tise.hpp tise_solver_main.cpp tise.cpp` that no code calls LAPACK anymore — only these two comments (plus one accurate historical reference at `tise.cpp:1256,1326`) still mention it. The new code never touches `H`/`S` (it builds separate dense `A`/`B`); the comments' *conclusion* ("H/S remain valid") still holds, but the stated *mechanism* is fabricated — and it's fabricated in exactly the function Finding 3.1 shows has a real, non-obvious pitfall, which makes a wrong mental model here more costly than usual.

### 4.3 The GoogleTest CMake fetch block is duplicated and has already diverged

`TISE/CMakeLists.txt` and `student-base-code/CMakeLists.txt` (the latter newly adding this block in this PR) carry a near-verbatim ~20-line GoogleTest `FetchContent` block. Confirmed:

```
$ grep -n FETCH_DEPENDENCIES TISE/CMakeLists.txt        # TISE_FETCH_DEPENDENCIES
$ grep -n FETCH_DEPENDENCIES student-base-code/CMakeLists.txt  # FETCH_DEPENDENCIES
```

The two copies already use different names for what should be the same toggle, within the same PR that introduced the duplication. Recommend factoring into a shared CMake module/macro so a future GTest version bump or CRT-linkage fix can't be applied to one copy and forgotten in the other.

### 4.4 Table-writer formatting and grid-point formula duplicated 5-8x

The `setprecision(16)`/`setw(24)` output idiom and the `(rMax-rMin)*(ix-1)/(npts-1)` grid formula are each independently copy-pasted across `writeEigenstate`, `writeEigenstateTable`, `writeContinuumInfo`, `writeContinuumTable`, `writePotential`, `writeEigenvalues`, `writeEigenvectors`, and `writeBandedMatrix`. This is also the literal site of Finding 3.3 — because the formula is duplicated rather than centralized, a fix applied to one call site is easy to miss at the others. A shared `linspace`-style helper (taking `npts` and asserting `npts >= 2`, which would also directly fix 3.3) would remove both the duplication and the missing-validation problem in one change.

### 4.5 Efficiency: redundant basis evaluation in two hot loops (reported by Pass 1, not independently re-measured)

- `fillBandedMatrices`' delta-term loop (`tise.cpp:745-747`, confirmed present by direct reading) recomputes `bs.eval(d.x, iBs2, 0)` once per inner `iBs1` iteration even though it does not depend on `iBs1` — roughly doubles this loop's basis-evaluation count for configs with delta terms.
- `writeEigenstateTable`/`writeContinuumTable` re-derive the full B-spline basis-value window at each grid point once per state/energy rather than once per point — O(`npts·nStates·order`) where O(`npts·order + npts·nStates`) would suffice.

Neither is a correctness bug; both are cheap, mechanical hoisting fixes if addressed.

### 4.6 No test exercises any of the bugs in Section 3

```
$ grep -c "potential_deltas\|DeltaTerm" tests/test_controller_unit.py tests/test_controller_integration.py
tests/test_controller_unit.py:0
tests/test_controller_integration.py:0
```

Confirmed: zero references. Equally, no C++ `SolveEigenTest` constructs a non-PD/singular `S`, no test sets `n_pts_eigenstate <= 1`, and `plot_spectrum_overview` has no dedicated test at all. None of Section 3's findings would be caught by the existing CI (`.github/workflows/ci.yml` runs the same `ctest`/`pytest` suites).

## 5. Reproduction

All commands below were run against a clean worktree of the PR's head commit (`433e11a`), checked out separately so the primary working tree (on `TISE-Generalization`) was never touched:

```bash
git worktree add <scratch-dir>/pr18-worktree pr-18   # pr-18 tracks PR #18's head, 433e11a
```

### 5.1 Confirming the base lacks the delta feature and the LAPACK call

```bash
$ git show 73229c87:1D-QM-Playground/TISE/tise.hpp | grep -c "DeltaTerm\|potential_deltas"
0
$ git show 73229c87:1D-QM-Playground/TISE/tise.cpp | grep -n "dsbgv_\|DSBGV"
24:    void dsbgv_(char *jobz, char *uplo,
1228:    dsbgv_(&jobz, &uplo,
1238:        throw std::runtime_error("DSBGV failed with info=" + std::to_string(info));
```

### 5.2 Eigen eigensolver behavior (Finding 3.1)

```cpp
// eigen_repro_test2.cpp
#include <Eigen/Eigenvalues>
#include <iostream>
void test(const std::string &label, Eigen::MatrixXd A, Eigen::MatrixXd B) {
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXd>
        solver(A, B, Eigen::ComputeEigenvectors | Eigen::Ax_lBx);
    std::cout << label << ": info==Success? " << (solver.info()==Eigen::Success)
               << "  eigenvalues: " << solver.eigenvalues().transpose() << "\n";
}
int main() {
    using Eigen::MatrixXd;
    MatrixXd I3 = MatrixXd::Identity(3,3);
    { MatrixXd B(3,3); B.setZero(); B(0,0)=1; B(1,1)=-1;    B(2,2)=1; test("1. Indefinite",      I3, B); }
    { MatrixXd B(3,3); B.setZero(); B(0,0)=1; B(1,1)=0;     B(2,2)=1; test("2. Exactly singular", I3, B); }
    { MatrixXd B(3,3); B.setZero(); B(0,0)=1; B(1,1)=1e-14; B(2,2)=1; test("3. Near-singular PD",  I3, B); }
    { MatrixXd B(3,3); B.setZero(); B(0,0)=1; B(1,1)=-1e-10;B(2,2)=1; test("4. Tiny-negative",     I3, B); }
    { MatrixXd B(3,3); B << 1,1,0, 1,1,0, 0,0,1;                      test("5. Rank-deficient",   I3, B); }
}
```

```bash
$ g++ -std=c++17 -I/usr/include/eigen3 eigen_repro_test2.cpp -o eigen_repro_test2 && ./eigen_repro_test2
1. Indefinite        : info==Success? 1  eigenvalues: 1 1 1
2. Exactly singular   : info==Success? 0  eigenvalues: -nan -nan -nan
3. Near-singular PD   : info==Success? 1  eigenvalues: 1 1 1e+14
4. Tiny-negative      : info==Success? 1  eigenvalues: 1 1 1e+20
5. Rank-deficient     : info==Success? 1  eigenvalues: 0.381966 1 2.61803
```

### 5.3 Building the real solver and reproducing the delta-collision crash (Finding 3.2)

```bash
$ cmake -S TISE -B TISE-build -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DTISE_FETCH_DEPENDENCIES=OFF
-- muparser: found via pkg-config
-- yaml-cpp: found on system
-- Configuring done

$ cmake --build TISE-build -j --target tise_solver
[100%] Built target tise_solver
```

`control.yaml` (single delta — sanity check it works):
```yaml
bspline: {n_nodes: 81, order: 8, domain: [-20.0, 20.0]}
potential: [{'domain': '[-20, 20]', 'function': '0'}]
potential_deltas: [{x: 0.0, strength: -1.0}]
tise: {n_pts_eigenstate: 1001, continuum: {enabled: false}}
```
```bash
$ ./tise_solver --config control.yaml --output-dir out-control && head -2 out-control/eigenvalues.dat
0  -5.0000000000002454e-01     # exact analytic E_0 = -0.5 for a unit delta well
```

`duplicate-delta.yaml` (same, with the one `potential_deltas` entry duplicated):
```bash
$ ./tise_solver --config duplicate-delta.yaml --output-dir out-dup
tise_solver: Eigen GeneralizedSelfAdjointEigenSolver failed; the overlap matrix is most likely not positive definite
exit=1
```

### 5.4 `n_pts_eigenstate: 1` (Finding 3.3)

Same config as `control.yaml` above with `n_pts_eigenstate: 1`:
```bash
$ ./tise_solver --config npts1.yaml --output-dir out-npts1
tise_solver: warning: 1 of 91 computed states are below E=0.0 (...)
exit=0

$ tail -1 out-npts1/eigenstates.dat
-nan -nan -nan -nan -nan -nan -nan -nan ...   # every one of 91 columns
```

### 5.5 `plot_spectrum_overview` empty-state crash (Finding 3.4)

```bash
$ python3 -c "
drawn = [(0, -0.5, [])]
peak = max((max(abs(p.phi) for p in st) for _, _, st in drawn), default=1.0) or 1.0
"
ValueError: max() iterable argument is empty
```

### 5.6 `potential_deltas: null` mismatch (Finding 3.5)

```bash
$ python3 -c "
import controller
controller.validate_potential_deltas(None, (-20.0, 20.0))
print('validation PASSED')
"
validation PASSED
```
`null-delta.yaml` (same as `control.yaml`, with `potential_deltas:` left empty):
```bash
$ ./tise_solver --config null-delta.yaml --output-dir out-null
tise_solver: 'potential_deltas' must be a YAML list of {x: ..., strength: ...} mappings
exit=1
```

### 5.7 `visualization: null` crash (Finding 3.6)

```bash
$ python3 -c "
import yaml
cfg = yaml.safe_load('visualization:\n')
cfg.get('visualization', {}).get('spectrum_overview', True)
"
AttributeError: 'NoneType' object has no attribute 'get'
```

### Environment

- `cmake` 3.x / `g++` 11.4.0 (system), Eigen3 (system, via `pkg-config`), yaml-cpp (system), muparser (system via `pkg-config`) — no `FetchContent` network access was needed.
- Python 3 with the project's existing `yaml` module; `controller`/`analysis` imported directly from the PR worktree's `1D-QM-Playground/` directory.

## 6. Recommendations / Implementation Plan

Ordered by what blocks merge:

1. **Fix Finding 3.1** — add an explicit `Eigen::LLT` positive-definiteness check (or equivalent) after the generalized eigensolver call, independent of `solver.info()`. This is the root dependency for 3.2's severity: once it's fixed, a delta collision will reliably fail loudly instead of sometimes failing loudly and sometimes failing silently.
2. **Fix Finding 3.2** — reject colliding/duplicate delta positions (and delta-on-join collisions) in `validateDeltaTerms` and `validate_potential_deltas`, with a tolerance consistent with the knot-insertion logic's existing `1e-12`.
3. **Fix Finding 3.3** — validate `n_pts_eigenstate >= 2` / `continuum.n_pts >= 2` in `controller.py`'s config validation, and defensively in `tise_solver_main.cpp`.
4. **Fix Findings 3.4-3.6** — each is a one-line change in `analysis.py`/`controller.py`; low risk, no architectural decision needed.
5. **Add regression tests for all of the above** (Finding 4.6) before merge — a non-PD-overlap `SolveEigenTest`, an `n_pts_eigenstate <= 1` config test, a duplicate-delta config test, and the two Python crash repros as unit tests — so this class of bug has a tripwire going forward.
6. Section 4's quality items (4.1-4.5) are not blockers; fold them in opportunistically while touching the same files for items 1-5, since 4.2 and 4.4 directly border fixes already being made.

## 7. Appendix: Full Finding List (ReportFindings output, for reference)

The structured finding list reported via the review tooling, for traceability against this narrative:

| # | File | Category | Severity (this doc's ordering) |
|---|---|---|---|
| 1 | `TISE/tise.cpp:1362` | correctness | Critical — Section 3.1 |
| 2 | `TISE/tise.cpp:1938` | correctness | Critical — Section 3.2 |
| 3 | `TISE/tise.cpp:1751` | correctness | High — Section 3.3 |
| 4 | `analysis.py:600` | correctness | High — Section 3.4 |
| 5 | `controller.py:202` | correctness | Medium — Section 3.5 |
| 6 | `analysis.py:676` | correctness | Low — Section 3.6 |
| 7 | `TISE/tise.hpp:540` | simplification | Section 4.2 |
| 8 | `TISE/tise.cpp:745` | efficiency | Section 4.5 |
| 9 | `TISE/tise.cpp:1751` | efficiency | Section 4.5 |
| 10 | `TISE/tise.cpp:1346` | simplification | Section 4.1 |
| 11 | `TISE/CMakeLists.txt:211` | reuse | Section 4.3 |
| 12 | `TISE/tise.cpp:1700` | reuse | Section 4.4 |
| 13 | `TISE/tests/test_tise.cpp` | test-coverage | Section 4.6 |
