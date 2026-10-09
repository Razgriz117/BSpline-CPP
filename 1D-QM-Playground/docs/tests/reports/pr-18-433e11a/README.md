# Verification of `docs/planning/PR-18-phy5606-projects-review.md`

**Subject:** the thirteen findings in the PR #18 review (`PHY5606_Projects` @ `433e11a`, base `main` @ `73229c8`).
**Method:** same as the known-solution reports — a clean checkout of `433e11a` was built in an isolated container (`cmake -DBUILD_TESTING=ON -DTISE_FETCH_DEPENDENCIES=OFF`, system Eigen 3.4.0 / yaml-cpp / muparser / LAPACK; all 4 ctest suites pass), every claimed reproduction was re-run from the YAML or program listed in `repro/`, and where the review reasoned about behaviour without measuring it, the measurement was added. Nothing below was taken from the review's own transcript.
**Date:** 2026-10-09.

## Verdict table

| # | Finding | Review says | Verified? | What the evidence shows |
|---|---|---|---|---|
| 3.1 | `GeneralizedSelfAdjointEigenSolver::info()` does not detect non-PD overlap | Critical | **Confirmed**, with two corrections | Eigen's `info()` reports only QR convergence; `LLT(B).info()` and LAPACK `dsbgv` (`info = n+2`) both flag the same matrices. Review's cases 3 and 7 ("spurious") are **not** detection failures — $10^{14}$ is the exact eigenvalue of that problem. And the proposed `LLT` fix **does not catch the real-solver corruption in 3.2** (see below). |
| 3.2 | Colliding delta positions degenerate the basis, unvalidated | Critical | **Confirmed and worse than stated** | Exact duplicate: caught (exit 1). Two deltas $10^{-9}$–$10^{-13}$ apart: **exit 0, $E_0 = -2\times10^{4}$ … $-2\times10^{12}$** (exact $-2.0$). Cholesky of $S$ *succeeds* in every one of these cases; the damage comes from $\mathrm{cond}(S) = 10^{13}$–$10^{17}$, not indefiniteness. |
| 3.3 | `n_pts_eigenstate: 1` → all-NaN file, exit 0 | High | **Confirmed**; e2e consequence added | `tise_solver` exits 0; `eigenstates.dat` and every `eigenstate_NNN.dat` are `-nan`. Through `controller.py` the run then fails in `analysis.py`'s reader with a clean message — so the pipeline is not silent, only the solver is. `n_pts_eigenstate: 0` writes zero-row files and reaches 3.4. |
| 3.4 | `plot_spectrum_overview` raises bare `ValueError` on an empty state | High | **Confirmed e2e** | `python3 controller.py --config npts0_full.yaml` ends in a raw traceback at `analysis.py:600`, not an `AnalysisError`. |
| 3.5 | `potential_deltas: null` passes Python, fails C++ | Medium | **Confirmed e2e** | `validate_config()` passes; `tise_solver` exits 1 with the "must be a YAML list" message; `controller.py` surfaces it as a solver-stage failure. Both `potential_deltas:` (empty) and explicit `null` behave identically. |
| 3.6 | `visualization:` (null) → `AttributeError` | Low, pre-existing | **Confirmed e2e**; pre-existence confirmed | Raw `AttributeError` traceback at `analysis.py:676` through the full pipeline. Same pattern at lines 678/681/683 predates the PR. |
| 4.1 | Upper-triangle writes are dead | quality | **Confirmed** | Eigen 3.4.0 header: `selfadjointView<Lower>()` at lines 184, 197, 210 — all three problem types. |
| 4.2 | Stale LAPACK comments | quality | **Confirmed** | `tise.hpp:540` and `tise_solver_main.cpp:469–471` as quoted; no `dsbgv_` symbol remains. |
| 4.3 | Duplicated, divergent GTest CMake block | quality | **Confirmed** | `TISE_FETCH_DEPENDENCIES` vs `FETCH_DEPENDENCIES`, same block. |
| 4.4 | Grid formula / formatting duplicated | quality | **Confirmed** | `(npts - 1)` formula at 5 sites (1245, 1701, 1751, 1808, 1829); `setprecision(16)` at 9. |
| 4.5 | Redundant `bs.eval` in delta loop | quality | **Confirmed by reading** (not timed) | `tise.cpp:746–747`: `bs.eval(d.x, iBs2, 0)` inside the `iBs1` loop. Cheap to hoist; negligible at $n\sim100$. |
| 4.6 | No test exercises the six bugs | test-coverage | **Confirmed for the bugs; overstated for the feature** | Zero Python references to `potential_deltas`, true. But `TISE/tests/test_tise.cpp` has **seven** `DeltaPotentialTest` cases (closed-form ground state, double-delta doublet to $10^{-9}$, out-of-domain rejection, …). What is missing is specifically: duplicate/near-duplicate deltas, non-PD or ill-conditioned $S$, `n_pts ≤ 1`, and `plot_spectrum_overview`. |

**Overall:** the review's six correctness findings are all genuine and all reproduce end-to-end. The most important adjustment is to its **fix plan**: item 1 (an `LLT` positive-definiteness check) would not have caught a single one of the silent-corruption cases that make finding 3.2 critical, because the overlap matrix in those cases is positive definite to Cholesky and merely ill-conditioned. Delta-position validation (item 2) is the fix that matters, and a conditioning guard on $S$ is the right second line of defence — not an `isfinite` check, which also passes on every corrupted output here.

---

## 1. The delta feature itself is physically sound (sanity baseline)

Before testing the failure modes, the feature's happy path was checked against closed forms (atomic units, $m=\hbar=1$):

| Configuration | solver $E_0$ | exact | rel. error |
|---|---|---|---|
| single $-\delta(x)$ on $[-20,20]$, order 8, 81 nodes | −0.50000000000025 | $-g^2/2 = -0.5$ | 5e-13 |
| same, state 2 | 0.012337005501363 | box $n=2$: $4\pi^2/2\cdot40^2 = 0.012337005501$ (odd state, unaffected by the delta) | 2e-14 |
| two deltas at $x=\pm1$ (even) | −0.61478253628788 | $\kappa = 1+e^{-2\kappa}$ → −0.614782536288 | 1e-13 |
| two deltas at $x=\pm1$ (odd) | −0.31745478527345 | $\kappa = 1-e^{-2\kappa}$ → −0.317454785274 | 2e-13 |

The knot-multiplicity treatment ($\text{order}-2$ extra copies at each delta, giving a $C^0$ basis there) is the same device that fixed the interior-singularity test in iteration 3, and it is exact here. The review's "10 significant figures" claim for the single delta is correct (13, in fact).

## 2. Finding 3.1 — eigensolver error detection

### 2.1 What was measured

`repro/eig.cpp` runs the review's six synthetic $(I, B)$ pairs plus one more through three independent checks: Eigen's `GeneralizedSelfAdjointEigenSolver::info()`, Eigen's `LLT<MatrixXd>(B).info()`, and LAPACK's `dsbgv_` (the routine the PR removed), on the same matrices:

| case | $B$ | GSAES `info==Success` | `LLT(B)` Success | `dsbgv info` | eigenvalues returned |
|---|---|---|---|---|---|
| 1 indefinite | diag(1,−1,1) | **true** | false | 5 (= n+2: not PD) | 1 1 1 |
| 2 exactly singular | diag(1,0,1) | false | false | 5 | nan nan nan |
| 3 near-singular PD | diag(1,1e-14,1) | true | true | 0 | 1 1 1e14 |
| 4 tiny negative | diag(1,−1e-10,1) | **true** | false | 5 | 1 1 1e20 |
| 5 rank-deficient | [[1,1,0],[1,1,0],[0,0,1]] | **true** | false | 5 | 0.382 1 2.618 |
| 6 same + 1e-16 noise | | **true** | false | 5 | 0.382 1 2.618 |
| 7 same + 1e-10 | | true | true | 0 | 0.5 1 2e10 |

### 2.2 Assessment

The mechanism the review describes is exactly right: `GeneralizedSelfAdjointEigenSolver.h:175` computes `LLT<MatrixType> cholB(matB)` and never reads `cholB.info()`; `info()` is inherited from the tridiagonal QR step. Cases 1, 4, 5, 6 are indefinite or singular, `dsbgv` and `LLT` both say so, and the PR's new code does not. **The regression relative to `main` is real.**

Two corrections to the table as published:

* **Rows 3 and 7 are not failures.** $B = \mathrm{diag}(1, 10^{-14}, 1)$ *is* positive definite, and $\lambda = 1/10^{-14} = 10^{14}$ is the exact solution of $x = \lambda Bx$. Every solver, including `dsbgv` (`info = 0`), returns it. Calling these "spurious" conflates ill-conditioning with indefiniteness — and that distinction turns out to be the crux of finding 3.2.
* **Row 1** is the cleanest demonstration: an indefinite $B$ yields a finite, innocuous-looking spectrum (1 1 1) with `info == Success`. LAPACK's answer to the same input is equally meaningless — the difference is that `dsbgv` *says so* (`info = 5`), which is the contract the PR lost.

## 3. Finding 3.2 — colliding delta positions

### 3.1 Exact duplicate (the review's reproduction)

`repro/dup.yaml` (two identical `{x: 0.0, strength: -1.0}` entries): `tise_solver` exits 1 with the "overlap matrix is most likely not positive definite" message. Confirmed. The mechanism is as described: `buildStrategicRadialGrid` inserts $\text{order}-2 = 6$ extra copies per knot entry, so the knot at 0 reaches multiplicity $1 + 6 + 6 = 13 > \text{order}$; `BSpline::init` then produces basis functions with empty support, $S$ acquires exact zero rows, and this happens to land in case 2 above (caught).

### 3.2 Near-duplicate (the review's prediction, now measured)

The review predicted that "a less exact collision … could instead silently return finite-looking wrong numbers." `repro/near1e-*.yaml` place the second delta at $x = \Delta x$:

| $\Delta x$ | exit | $E_0$ (solver) | exact ($\to -g_{\rm tot}^2/2 = -2$ as $\Delta x\to0$) | min eig $S$ | cond $S$ | `LLT(S)` |
|---|---|---|---|---|---|---|
| $10^{-3}$ | 0 | −1.9960099734 | −1.9960099734 (even state of two unit deltas $\Delta x$ apart: $\kappa = 1 + e^{-\kappa\Delta x}$) ✓ | 3.3e-08 | 1.5e07 | Success |
| $10^{-6}$ | 0 | −1.9999959969 | −1.9999960000 ✓ (3e-9; the basis is starting to strain) | 3.3e-11 | 1.5e10 | Success |
| $10^{-9}$ | 0 | **−2.07 × 10⁴** | −2.0 | 3.3e-14 | 1.5e13 | Success |
| $10^{-11}$ | 0 | **−2.40 × 10⁸** | −2.0 | 3.3e-16 | 1.5e15 | Success |
| $10^{-13}$ | 0 | **−2.12 × 10¹²** | −2.0 | −7.7e-18 | 6e16 | Success |

(`figures/near_duplicate_deltas.png`.) Three things the review did not have:

1. **The corruption is not "plausible-looking".** It is four to twelve orders of magnitude off, with four negative eigenvalues where there should be one, and the solver's only output is a routine wall-collision warning. Exit code 0.
2. **It is not caught by Cholesky.** `repro/llt.cpp` reads the solver's own `overlap.dat` back and runs `Eigen::LLT` on it: `info() == Success` for every $\Delta x$ down to $10^{-13}$, even where the smallest eigenvalue has gone (barely) negative. The review's recommended fix for 3.1 — "explicitly check `Eigen::LLT<MatrixXd>(B).info()`" — would therefore change nothing here. What goes wrong is that the two knot clusters at 0 and $\Delta x$ generate B-splines of support width $\Delta x$, whose overlap entries are $O(\Delta x)$; $S$ stays positive definite but its condition number scales as $1/\Delta x$ and crosses $1/\varepsilon_{\rm machine}$ near $\Delta x \approx 10^{-11}$. Already at $10^{-9}$ (cond $10^{13}$) the Cholesky-reduced problem has lost the digits that carry the physics.
3. **The reference the review used is wrong by a factor of four.** It compares the duplicated config to the single-delta control ($E_0 = -0.5$). Two coincident unit deltas are one delta of strength 2, $E_0 = -2$; the $\Delta x = 10^{-3}$ and $10^{-6}$ rows above reproduce exactly that, which is the real control for this experiment.

### 3.3 Delta on an existing join

Not tested by the review and not tested here: the same knot-stacking logic (`tise.cpp:1938` for Singular joins, `:1952` for deltas) would sum multiplicities for a delta placed at a `Step`/`StitchedKink`/`Singular` join. The mechanism is identical to §3.1, so the conclusion carries over; a regression config for it is listed in §8.

## 4. Finding 3.3 — `n_pts_eigenstate ≤ 1`

| config | `tise_solver` exit | `eigenstates.dat` | per-state files | via `controller.py` |
|---|---|---|---|---|
| `n_pts_eigenstate: 1` | 0 | header + one row of 92 `-nan` | each one `-nan` row | exit 1: `analysis.py: malformed eigenstate_001.dat … non-finite field` |
| `n_pts_eigenstate: 0` | 0 | header only | zero rows | exit 1: raw `ValueError` traceback (= finding 3.4) |

Confirmed. The formula `(rMax − rMin)·(ix−1)/(npts−1)` at five sites (§4.4) divides by zero for `npts = 1`; for `npts = 0` the loop body never runs. Neither `controller.py` (no `n_pts` validation at all — `grep -n n_pts controller.py` is empty) nor `tise_solver_main.cpp:272` (`.as<int>()` with no range check) guards it. One nuance the review's "silently reported as successful" omits: in the *shipped pipeline* the run is not silent — `analysis.py`'s reader (`_read_data_rows`) rejects the NaN and the controller exits 1 with a clear message. The silence is confined to anyone using `tise_solver` directly or with `run_analysis: false`, which the PR's own README recommends for custom analysis.

## 5. Findings 3.4 – 3.6 — Python crashes and the validation gap

All three reproduced both as isolated expressions (the review's method) and end-to-end through `controller.py` (`repro/npts0_full.yaml`, `repro/visnull.yaml`, `repro/nulldelta_full.yaml`):

* **3.4** `analysis.py:600` — `max(abs(p.phi) for p in st)` with no `default` on an empty `st`. Raw `ValueError` traceback; not an `AnalysisError`. Reachable from `n_pts_eigenstate: 0`.
* **3.5** `controller.validate_config()` passes a config whose `potential_deltas` is `None` (both `potential_deltas:` and `potential_deltas: null` parse to `None`); `tise_solver` then exits 1. yaml-cpp's null node is `IsDefined()` but not `IsSequence()`, exactly as described. The controller's docstring promise ("gives the error before a subprocess is launched") is broken for this input.
* **3.6** `analysis.py:676` — `cfg.get("visualization", {})` returns `None` for a present-but-null key; `.get` on it raises. Confirmed pre-existing at lines 678/681/683 in the base commit. The new line runs first and defaults to `True`, so it is the first to trip.

## 6. Code-quality findings 4.1 – 4.6

Each was checked against the source rather than the diff:

* **4.1** Confirmed: Eigen 3.4.0's `GeneralizedSelfAdjointEigenSolver.h` uses `selfadjointView<Lower>()` in all three branches (lines 184, 197, 210). The `A(i,j)`/`B(i,j)` upper writes at `tise.cpp:1348–1349` are never read. Harmless; worth one line of comment if kept.
* **4.2** Confirmed verbatim (`tise.hpp:540`, `tise_solver_main.cpp:469–471`). `dsbgv_` appears nowhere in code; only in comments.
* **4.3** Confirmed: identical fetch block, option spelled `TISE_FETCH_DEPENDENCIES` in `TISE/CMakeLists.txt:57` and `FETCH_DEPENDENCIES` in `student-base-code/CMakeLists.txt:19`.
* **4.4** Confirmed: the grid formula appears at `tise.cpp` 1245, 1701, 1751, 1808, 1829; `setprecision(16)` nine times. A `linspace(npts ≥ 2)` helper fixes 3.3 at all sites at once.
* **4.5** Confirmed by inspection (`tise.cpp:746–747`); not timed. For $n \approx 100$ and a handful of deltas the cost is microseconds — correctness of the hoist is the only reason to do it.
* **4.6** Partly overstated. The C++ delta feature **is** tested: `TISETests` contains seven `DeltaPotentialTest` cases, including the double-delta doublet to $10^{-9}$ and rejection of an out-of-domain delta. The accurate statement is that no test covers any of the six *failure modes* in §3 — duplicate or near-duplicate deltas, an indefinite or ill-conditioned $S$, `n_pts ≤ 1`, an empty state in `plot_spectrum_overview`, `potential_deltas: null`, `visualization: null` — and that the Python side has zero coverage of `potential_deltas` at all.

## 7. Assessment of the review's fix plan

| Review item | Keep? | Comment |
|---|---|---|
| 1. `LLT` PD check after the solve | keep, **demote** | Restores `main`'s behaviour for truly indefinite $S$ (cases 1, 4, 5, 6) — correct and cheap. But it does not touch the §3.2 corruption, so it is not "the root dependency for 3.2's severity" as the review states. |
| 2. Reject colliding deltas / delta-on-join in `validateDeltaTerms` + `validate_potential_deltas` | **promote to first** | The only change that prevents the garbage runs in §3.2. Tolerance: the data say $\Delta x < 10^{-6}$ is already physically indistinguishable from a merged delta and $\Delta x \le 10^{-9}$ is corrupt; a tolerance of order the uniform node spacing $\times10^{-3}$, or simply "merge deltas closer than `1e-6` into one of summed strength", is defensible. |
| (new) Conditioning guard on $S$ | **add** | After filling $S$, compute its extreme eigenvalues (or `LDLT` with `rcond`) and refuse when $\lambda_{\min} \le 0$ or $\lambda_{\max}/\lambda_{\min} > 10^{12}$. This catches every row of §3.2 and any future basis pathology, independent of its source. An `isfinite` check, which the review also suggests, catches **none** of them. |
| 3. `n_pts ≥ 2` validation both sides | keep | Add to `validate_config` and to `tise_solver_main.cpp:272/313`; or centralise in the `linspace` helper of 4.4. |
| 4. One-line fixes for 3.4–3.6 | keep | `default=0.0`; reject or normalise `None` consistently; `(cfg.get("visualization") or {})` at four sites. |
| 5. Regression tests | keep, with the list in §8 | |

## 8. Regression configs produced by this verification

All in `repro/`, runnable as `tise_solver --config <file> --output-dir <dir>` (or via `controller.py` for the `*_full.yaml` ones). Expected behaviour *after* the fixes:

| file | expected |
|---|---|
| `control.yaml` | $E_0 = -0.5 \pm 10^{-12}$ |
| `two-sep.yaml` | $E_0 = -0.6147825363$, $E_1 = -0.3174547853$ ($\pm10^{-10}$) |
| `dup.yaml`, `near1e-13.yaml`, `near1e-9.yaml` | rejected at validation with a message naming both deltas |
| `near1e-6.yaml` | either rejected (if tolerance ≥ 1e-6) or $E_0 = -1.999996 \pm 10^{-6}$ |
| `npts1.yaml`, `npts0.yaml` | rejected at validation (`n_pts_eigenstate >= 2`) |
| `nulldelta.yaml`, `nulldelta2.yaml` | same outcome from `validate_config` and from `tise_solver` |
| `visnull.yaml` | runs; `visualization: null` treated as `{}` |
| `eig.cpp` | reference for a `SolveEigenTest` asserting that cases 1, 4, 5, 6 throw |
| `llt.cpp` | demonstrates that `LLT` alone is insufficient on the real `overlap.dat` files |

## 9. Recommendations beyond the review

1. **Merge coincident deltas instead of rejecting them** if the use case is a student sweeping delta positions: two deltas within tolerance become one of summed strength, with a warning. Physically exact; avoids a hard failure on a benign input.
2. **Make `validate_config` the single source of truth** by having `tise_solver_main.cpp` call the same checks (or a shared schema): 3.3 and 3.5 are both instances of the two validators drifting.
3. **Keep `dsbgv` reachable behind a CMake option** on platforms that have LAPACK. It is the only solver here whose `info` reports both convergence and definiteness, and it is banded; the review's note that the dense Eigen path is "O(n³) rather than exploiting the band" will matter the first time someone sets `n_nodes: 2001`.
