# TDSE Solver — Task Breakdown for Two Engineers

**Date:** 2026-08-21
**Status:** Draft — ready for two-person parallel execution
**Scope:** SDD §5.3 ("TDSE Solver") internal numerics only — i.e. SDD §10.2 **Phase 8** ("TDSE Solver implementation"). This document does **not** cover Phases 5–7 (Controller↔TDSE, TDSE↔Analysis, Controller↔Analysis-extended interface contracts — the `tdse_solver` CLI binary, `yaml-cpp` config ingestion, the formal `data/tdse/` output-directory contract including `warnings.json`) or REQ-F-060 Analysis-module work (`analysis.py`). Mirrors exactly how `docs/planning/tise-task-breakdown.md` scoped itself to TISE's Phase 4.

---

## 1. Context & Scope

### Why this document exists

`docs/planning/engineer-a-plan-cleanup.md` (item 6) names this document directly: *"Create plan for TDSE for engineers A and B."* The SDD (`docs/SDD.md` §5.3, Figure 8, Figure 9) and `docs/planning/bsplines.md` ("The Round-Trip: Enabling TDSE") together describe the TDSE solver's architecture, but — unlike TISE's Figure 6/7, which fully specified both branches of that solver's remaining numerics — **Figure 8 leaves the driven-case propagator method explicitly unresolved**: *"numerically integrate forward in time (propagator scheme is an implementation detail, not user-configurable)."* That gap is now closed by **`docs/adr/0006-tdse-propagator-choice.md`**, which compares four candidate schemes and adopts one (§"The propagator decision" below summarizes it). This document does two things TISE's didn't have to: divide the remaining work into two parallelizable workstreams, **and** carry the task-level consequences of ADR-0006's propagator decision through A1–A6/B1–B6.

### What's already done (do not re-implement)

Confirmed by direct code read of `TISE/time_evolution.hpp`/`.cpp` and exercised by `TISE/tests/test_time_evolution.cpp` (~15 passing `TEST`/`TEST_F` cases):

- `tevol::gaussianWavepacket`, `computeGaussianOverlaps`, `projectToEigenBasis`, `timeEvolveState`, `transformToSpaceBasis`, `writeTimestep`, `runTimeEvolution` (`time_evolution.cpp:22–154`) implement the **field-free** round trip of Figure 9 exactly: Gaussian initial state → B-spline overlaps → eigenbasis projection (`Phi_G = Cinv * B_G`) → per-eigenstate phase `e^{-iE_nt/ħ}` → back-transform to B-spline coefficients → real-space evaluation.
- This is genuinely exact numerics for the field-free case (Figure 8's `STATIC` node) — no time-stepping, no truncation error, just a per-snapshot closed-form phase evaluation.
- The change-of-basis machinery it depends on (`EigenResult{values, vectors, dim, ldz}`, `Eigen::Map<const Eigen::MatrixXd> C(er.vectors.data(), nEn, nEn)`) is TISE output, already built and tested.

### What's genuinely missing

Confirmed by exhaustive grep across the whole repo for `TDSE|TimeDependent|Propagator|CrankNicolson|SplitOperator|Magnus|Chebyshev|Lanczos|van Dijk`: **zero hits outside the archived, superseded `docs/TDSE-original-design/`.** Figure 8's `DRIVEN` branch ($H(t)=H_0+H_\text{int}(t)$, genuine population transfer) has no code at all — not a stub, not a placeholder. Neither does any dipole/momentum operator assembly, field-expression evaluation, generalized (momentum-carrying) initial-state construction, or observable computation (norm/energy/population) — `tevol::` only ever writes raw $\psi(x,t)$ snapshots, nothing derived from them.

### The propagator decision

Four candidate schemes were compared for the driven case ($H(t)=H_0+H_\text{int}(t)$, $H_\text{int}(t)=-\hat O\,\mathcal E(t)$); full derivations, two independent equivalence proofs, and a concrete efficiency table are in **`docs/adr/0006-tdse-propagator-choice.md`**. Summary:

- **Dense Crank-Nicolson (Cayley), eigenbasis.** Exactly unitary, standard — but each step needs a fresh dense $n_\text{En}\times n_\text{En}$ complex factorization, since the field changes every step and nothing can be cached between steps.
- **Generalized (banded) Crank-Nicolson, B-spline basis.** Proven *exactly* equivalent to the eigenbasis version above (two independent proofs: a Löwdin/$\mathbf S^{1/2}$ substitution, and a direct substitution of TISE's own $C$ matrix showing the two produce identical trajectories) — keeps the B-spline basis's banded structure (`TISE/BSpline.hpp`'s compact support), at the cost of a new general-banded LAPACK binding.
- **Symmetric split-operator (Strang splitting), eigenbasis — adopted.**
  $$|\psi(t+\Delta t)\rangle = e^{-i\hat H_0\Delta t/(2\hbar)}\;e^{+i\,\Delta t\,\mathcal E(t+\Delta t/2)\,\hat O/\hbar}\;e^{-i\hat H_0\Delta t/(2\hbar)}\;|\psi(t)\rangle$$
  The two free half-steps are exactly the existing, tested `tevol::timeEvolveState` called with $t=\Delta t/2$ — zero new numerics. The interaction step is exact once $\hat O_\text{eig}$'s eigendecomposition ($\hat O_\text{eig}=W\Lambda W^\dagger$) is precomputed **once**, up front, via `Eigen::SelfAdjointEigenSolver` — no per-step factorization at all.
- **Split-operator + Crank-Nicolson hybrid.** Considered and rejected — even a "smart" implementation collapses to the split-operator method's own cost profile while giving up its exactness, for no offsetting benefit.

**Why split-operator wins, not just "is also valid":** unlike the Crank-Nicolson family, this method reduces to the already-built field-free `STATIC` path **exactly** (bit-for-bit, by the exponential group law) when the field is zero — not merely to $O(\Delta t^2)$, a strictly stronger correctness property. It is also the cheapest of the four options at this project's actual scale ($n_\text{En}\approx59$, ~1000 steps per default `config.yaml`) — roughly 59× fewer flops than dense CN — with **no accuracy tradeoff**, and needs zero new dependencies (unlike the banded-CN alternative's new LAPACK binding). See ADR-0006 for the full comparison and the deferred alternatives' revisit triggers.

### How the split works

A literal fork at Figure 8's `FIELDQ{field.enabled?}` node — the way TISE forked at its `CLASSIFY` node — would be **unbalanced**: `STATIC` is fully implemented (above) and `DRIVEN` is 100% of the missing work. Splitting there would hand one engineer polish and the other everything.

Instead, the split cuts **across** Figure 8: one engineer owns the spine every state vector flows through (regardless of which branch); the other owns the physical inputs that spine consumes but never touches the spine itself.

```
IN ──▶ INIT{alpha0} ──▶ FIELDQ{field.enabled?}
        ▲                   │no                 │yes
        │                   ▼                   ▼
   built by B5          STATIC (A2)          DRIVEN sandwich driver (A3)
   (Engineer B)                            half-free ─ kick ─ half-free
                                              ▲                ▲
                                              │ diagonalized   │ evaluated
                                              │ once (A1)      │ every step
                                          Oeig ┘            fieldEval(t) ┘
                                              ▲                ▲
        │                          produced by B4          produced by B3
        └──────────────────┬───────────────────────────────────┘
                            ▼
                  RECORD: norm, population, <H0> (A4)
                            ▼
                  WRITE: snapshot_NNNNN.dat, observables.dat (A5)
                  norm-drift check (A6)
```

**The two handoffs**, agreed up front and never re-litigated (per `docs/adr/0006-tdse-propagator-choice.md`'s Consequences — the adopted split-operator method wants two plain-function handoffs, not one fused closure):

```cpp
Eigen::MatrixXcd buildCouplingOperatorEigenbasis(const bspline::BSpline&, const tise::EigenResult&,
                                                  const std::string &gauge, double hbar);  // B4 -> A1
double evaluateField(const std::string &expression, double t);                            // B3 -> A3
```

Engineer B's code *produces* both. Engineer A's code *consumes* both — A1 diagonalizes the first once, outside the time loop; A3 calls the second fresh every step — and never constructs either in production code, only synthetic values in A's own tests. **Every task on both sides is independently buildable and testable without the other engineer's real implementation** — Engineer A validates against hand-built toy Hermitian matrices and an analytically-solvable two-level Rabi system; Engineer B validates against brute-force `bs.integral` calls and a real, already-working `solveTISE` run. This is a stronger independence guarantee than TISE's own split had (which still needed one soft "grid stability" coordination point) — here there's no runtime coupling at all, only an up-front type/sign-convention agreement (§5).

---

## 2. Engineer A — Propagation Numerics, Observables & Output

New files only: `TISE/propagator.hpp`/`.cpp`, namespace `tevol::`. **`time_evolution.hpp`/`.cpp` are not modified** — every existing tested function stays exactly as-is; A2/A3 call into them.

### A1. Coupling-operator eigendecomposition + exponential-apply primitive

**What:** Given a fixed Hermitian $\hat O_\text{eig}$ (from Engineer B's B4), diagonalize it **once** via `Eigen::SelfAdjointEigenSolver` to get $\hat O_\text{eig}=W\Lambda W^\dagger$, then provide a cheap per-call primitive applying $e^{i\theta\hat O_\text{eig}}$ for any scalar $\theta$:

```cpp
struct DiagonalizedOperator { Eigen::MatrixXcd W; Eigen::VectorXd lambda; };

DiagonalizedOperator diagonalizeCouplingOperator(const Eigen::MatrixXcd &Oeig);

Eigen::VectorXcd applyExponential(const DiagonalizedOperator &diag, double theta,
                                   const Eigen::VectorXcd &v);
// diag.W * ( (i*theta*diag.lambda).array().exp().matrix().asDiagonal() * (diag.W.adjoint() * v) )
// i.e. W diag(e^{i*theta*lambda_k}) W^dagger v — project onto W's basis, phase, rotate back
```

Per `docs/adr/0006-tdse-propagator-choice.md`'s sign derivation, callers pass $\theta=\Delta t\cdot\text{fieldValue}/\hbar$ and get $e^{+i\theta\hat O_\text{eig}}v$ — the `+i` (not `-i`) sign is load-bearing and is this primitive's responsibility, not Engineer B's.

**Source:** `docs/adr/0006-tdse-propagator-choice.md` (Decision).

**Touches:** new file only; consumes B4's `buildCouplingOperatorEigenbasis` output as `Oeig`.

**Done when:** (a) for a random synthetic Hermitian matrix, `applyExponential` at several $\theta$ values matches Eigen's generic `(std::complex<double>(0,1)*theta*Oeig).exp()` (via `unsupported/Eigen/MatrixFunctions`, used only in the test, never in production code) to machine precision; (b) $\lVert\text{applyExponential}(\cdot)\rVert_2$ is exactly norm-preserving for any $\theta$ (unitarity — a product of unit-modulus phases can't change a vector's norm); (c) for the exactly-solvable two-level toy system $\hat O_\text{eig}=\begin{pmatrix}0&V\\V&0\end{pmatrix}$, a resonant $\theta$ chosen to give a known rotation angle reproduces the closed-form rotated state directly. **Zero dependency on Engineer B's code** — synthetic matrices suffice for full coverage.

### A2. STATIC-path snapshot dispatch

**What:** Figure 8's `no` branch — the exact-phase formula itself does not change (Figure 8 is explicit: "no time-stepping needed"). What's missing is a snapshot-time dispatcher around it: given `alpha0` (built by B5) and a snapshot cadence, evaluate `tevol::timeEvolveState` directly at each requested time (an $O(1)$ exact evaluation per snapshot, no incremental stepping). A small shared helper, reused by A3:

```cpp
std::vector<double> snapshotTimes(double dt, double tFinal, int snapshotInterval);
```

**Source:** SDD Figure 8 (`STATIC` node); SDD §6.1 `tdse.snapshot_interval` (currently unused by `runTimeEvolution`, which writes every step — `time_evolution.cpp:133`).

**Touches:** calls (does not modify) `tevol::timeEvolveState`, `tevol::transformToSpaceBasis`.

**Kept as an explicit fast path, with a note on why that's now optional rather than required:** under the adopted split-operator method (ADR-0006), A3's driven-case output with a zero field is now **exactly** — not approximately — equal to this task's output, since the two half-free steps compose by the exponential group law to the same closed-form phase `timeEvolveState` already computes. That makes A2 no longer strictly load-bearing for correctness (A3 alone could serve both cases); keeping it separate is a legitimate choice (already-written, already-tested code, and it doubles as an exact-match oracle for A3's own "Done when (a)" below) but merging it away later is now a legitimate simplification too — which was **not** true under the Crank-Nicolson approach originally considered (there the match was only $O(\Delta t^2)$, making A2 load-bearing, not optional). Whoever refactors this later should know which case applies.

**Done when:** for `initial_state.type: eigenstate, index: k`, population stays exactly 1.0 in state $k$ and 0.0 elsewhere at every snapshot time; for a Gaussian init, sparse-snapshot output matches evaluating `timeEvolveState` at every single step and discarding non-snapshot steps — confirming the snapshot cadence introduces no discretization error at all (it can't, since each evaluation is independent and exact).

### A3. Driven-case sandwich trajectory driver

**What:** The `DRIVEN` loop, implementing ADR-0006's adopted split-operator method — per step: half-free evolution (`timeEvolveState` at $t+\Delta t/2$), exact interaction kick (A1's `applyExponential` with $\theta=\Delta t\cdot\text{fieldEval}(t+\Delta t/2)/\hbar$), then another half-free evolution; record every `snapshotInterval`-th step:

```cpp
struct TdseSnapshot { double t; Eigen::VectorXcd alpha; };

std::vector<TdseSnapshot> propagateDriven(const Eigen::VectorXcd &alpha0,
                                           const std::vector<bspline::Real> &diagE,
                                           const DiagonalizedOperator &OeigDiag,
                                           const std::function<double(double)> &fieldEval,
                                           double dt, double tFinal,
                                           int snapshotInterval,
                                           double hbar = 1.0);
```

**Source:** `docs/adr/0006-tdse-propagator-choice.md` (Decision); SDD Figure 8 `DRIVEN`/`POP` nodes.

**Note for Engineer B:** this consumes B4's `Oeig` (diagonalized once via A1, *outside* the step loop, not per step) and B3's `evaluateField` (called fresh *every* step) as two plain-function handoffs, not a fused closure. Confirm this shape with Engineer A before either side is deep into implementation.

**Touches:** new file; consumes A1's `diagonalizeCouplingOperator`/`applyExponential` and `tevol::timeEvolveState` (called twice per step, for the two half-steps).

**Done when:** (a) with `fieldEval` returning zero for all $t$, output at every snapshot matches A2's `STATIC` output to near machine precision — an **exact** cross-check, not merely a converging one, per A2's note above; (b) for a resonant, constant-amplitude `fieldEval` on the two-level toy system, reproduces the Rabi rotation already validated in A1; (c) snapshot count/timing matches the requested cadence exactly. **Zero dependency on Engineer B's real dipole/field code** — synthetic `fieldEval`/`Oeig` inputs suffice for full coverage.

### A4. Observables from a trajectory

**What:**
```cpp
double computeNorm2(const Eigen::VectorXcd &alpha);
double computeFieldFreeEnergyExpectation(const Eigen::VectorXcd &alpha,
                                          const std::vector<bspline::Real> &diagE);
double computeBoundPopulation(const Eigen::VectorXcd &alpha,
                               const std::vector<int> &boundIndices);
```

**Interpretive note (state this explicitly in code comments, don't silently guess):** SDD §6.3's `observables.dat` 4th column, $P(t)$, is not fully specified anywhere, and neither is which "energy" the 3rd column means. This task uses $\langle H_0\rangle(t)=\sum_n|\alpha_n(t)|^2E_n$ (population-weighted field-free energy) rather than the full instantaneous $\langle\psi(t)|H(t)|\psi(t)\rangle$ — deliberately, so A4 needs only $\alpha(t)$ and the already-available $E_n$, with **no dependency on Engineer B's operators**. `boundIndices` is a plain parameter (no persisted bound/continuum classification exists in `data/tise/` output to read it from — §4 cleanup item), mirroring TISE's own precedent of plain-parameter, no-config-wiring tasks (A5/B2 in `tise-task-breakdown.md`). Note the notational collision: this $P(t)$ is unrelated to the interval-probability $P_{[x_a,x_b]}(t)$ used elsewhere in the SDD's Analysis section — same letter, different quantity.

**Source:** SDD §6.3 storage table; `docs/planning/architecture-06-26.md` line 103 ("compute escape probability P(t) and survival probability S_n(t)").

**Touches:** new file; consumes `EigenResult` only.

**Done when:** for field-free evolution, `computeFieldFreeEnergyExpectation` is *constant* across all snapshots (energy conservation — a strong invariant, not just a tolerance check); `computeNorm2` matches an independently-computed brute-force $\int|\psi(x,t)|^2dx$ (via `bs.integral` on the back-transformed B-spline coefficients) to quadrature precision — the single most important check in this document, since it directly validates the "S is identity in the eigenbasis" claim end-to-end rather than by assertion.

### A5. Output formatting

**What:**
```cpp
std::string snapshotFilename(int index); // "snapshot_" + 5-digit zero-pad + ".dat"
void writeSnapshot(std::ostream &out, const bspline::BSpline &bs,
                    const Eigen::VectorXcd &psiBSplineCoeffsPadded,
                    int nBSplines, int npts, bspline::Real xMin, bspline::Real xMax);
void writeObservablesRow(std::ostream &out, double t, double norm2,
                          double energy, double population);
```
matching SDD §6.3's `snapshot_NNNNN.dat` (x, Re ψ, Im ψ) / `observables.dat` (t, norm, ⟨E⟩, P(t)) column contracts, same `std::scientific`/`setprecision(16)` convention as `tise::writeEigenstate` (`tise.cpp:719–736`) and `tevol::writeTimestep` (`time_evolution.cpp:82–109`).

**Source:** SDD §6.3 (verbatim file-format table).

**Explicitly not in scope:** wiring through `--output-dir`/`data/tdse/` — that's Phase 5 (excluded, §6). Write plain functions callable from tests, mirroring TISE B4's identical exclusion.

**Done when:** filenames match `snapshot_\d{5}\.dat`; `observables.dat` rows have exactly 4 columns in the documented order.

### A6. Norm-drift warning (detection only)

**What:**
```cpp
struct NormDriftCheck { bool driftedBeyondTolerance; double deviation; };
NormDriftCheck checkNormDrift(double norm2, double referenceNorm2 = 1.0,
                               double tol = 1e-6, std::ostream &warnOut = std::cerr);
```
mirroring `tise::classifyAsymptote`'s `warnOut` parameter convention (`tise.hpp:120–123`).

**Source:** SDD §5.3.4 ("norm drift beyond tolerance during propagation... a correctness warning"); §8's warning taxonomy.

**Explicitly not in scope:** writing `data/tdse/warnings.json` (Phase 5, mirroring how the equivalent TISE sidecar was established as part of TISE's Phase-1 contract, not its Phase-4 numerics) and **any corrective action** (renormalization, adaptive step-size) — detection/warning only, per SDD §5.3.4's own framing.

**Done when:** an artificially corrupted trajectory (constructed only inside the test, e.g. via a deliberately non-unitary synthetic step) triggers the warning at the documented tolerance boundary; a cleanly-propagated trajectory (via A2 or A3) never does.

---

## 3. Engineer B — Physical Operator & Initial-State Construction

New files only: `TISE/coupling_operator.hpp`/`.cpp`, `TISE/field.hpp`/`.cpp`, `TISE/initial_state.hpp`/`.cpp`, namespace `tevol::`.

**Scope assumption (see §6):** all tasks below assume the *default* B-spline drop-set (`{1, nBSplines}`, i.e. `nEn = nBSplines - 2`) — the same assumption `tevol::computeGaussianOverlaps` already makes (`time_evolution.cpp:42`: loops physical indices `iBs2 = 2..nEn+1`). This is deliberately narrower than `tise::fillBandedMatrices`'s generalized `dropSet`-aware loop (`tise.cpp:345–398`, using `resolveDropSet`/`colOf`) — supporting TISE's generalized drop-set here is an explicit, flagged follow-up, not attempted in this task.

### B1. Dipole operator (length gauge) — B-spline assembly + eigenbasis rotation

**What:** $x_{ij} = \int B_i(x)\,x\,B_j(x)\,dx$ via `bs.integral(fX, iBs1, iBs2)` with `fX = [](double x, const double*){ return x; }`, looped exactly like `computeGaussianOverlaps` (`time_evolution.cpp:42–50`: physical indices `iBs2 = 2..nEn+1`, `iBs1 = max(2, iBs2-order+1)..iBs2`) into a dense `Eigen::MatrixXd` (set both `(i,j)` and `(j,i)` explicitly — this is plain dense storage, not LAPACK band storage, so symmetry isn't implicit). Rotate into the eigenbasis using the same `C` construction `runTimeEvolution` already uses (`time_evolution.cpp:128`):

$$O_\text{eig} = C^T\,O_\text{BS}\,C \qquad \text{(per \texttt{docs/planning/bsplines.md}, "Eigenstates as Linear Transformations")}$$

```cpp
Eigen::MatrixXd buildDipoleOperatorEigenbasis(const bspline::BSpline &bs,
                                               const tise::EigenResult &er);
```

**Source:** REQ-F-010 (operator scope); SDD Figure 8 ($\hat O=\hat x$, length gauge); `docs/planning/bsplines.md` line 393's rotation formula.

**Touches:** `bs.integral`, `EigenResult::vectors` — read-only use of existing TISE output.

**Done when:** for a harmonic-oscillator potential run through a real `solveTISE` (e.g. `{"domain":"(-inf,inf)","function":"0.5*x^2"}`), $|O_{\text{eig}}[n,n+1]|$ matches the closed form $\left|\langle n|\hat x|n{+}1\rangle\right|=\sqrt{(n+1)\hbar/(2m\omega)}$ within `test_tise.cpp`'s existing tolerance, and off-tridiagonal elements are ≈0 (selection rule). **Check magnitude, not sign** — LAPACK `DSBGV` fixes each eigenvector only up to an arbitrary overall sign, so an exact-signed comparison is fragile by construction; the sign-robust cross-check is B6. Structural check: result is real-symmetric to machine precision.

### B2. Momentum operator (velocity gauge) — B-spline assembly + eigenbasis rotation

**What:** $p_{ij} = -i\hbar\int B_i(x)\,B_j'(x)\,dx$. Compute $K_{ij}=\int B_i B_j'\,dx$ via `bs.integral(fUnity, iBs1, iBs2, 0, 1)` (same physical-index loop as B1), **explicitly antisymmetrize** $K_{ij}\leftarrow(K_{ij}-K_{ji})/2$ (robust against quadrature-order edge cases rather than relying on integration-by-parts holding to full floating-point precision), then $p_\text{BS}=-i\hbar K$, rotate as in B1.

```cpp
Eigen::MatrixXcd buildMomentumOperatorEigenbasis(const bspline::BSpline &bs,
                                                  const tise::EigenResult &er,
                                                  double hbar = 1.0);
```

**Source:** REQ-F-010; SDD Figure 8 ($\hat O=\hat p_x$, velocity gauge).

**Touches:** same pattern as B1.

**Done when:** $|O_\text{eig}[n,n+1]|$ matches $\sqrt{(n+1)m\hbar\omega/2}$ (magnitude check, same LAPACK sign-convention caveat as B1); structural checks provable directly from the antisymmetrized construction: real part exactly zero, diagonal exactly zero; Hermiticity $O_\text{eig}^\dagger=O_\text{eig}$ to machine precision.

### B3. Field-expression evaluation

**What:** reuses `evaluateFunction`'s `mu::Parser` mechanism (`tise.cpp`), simplified — `tdse.field.expression` is a single global expression in `t`, not a piecewise domain list:

```cpp
double evaluateField(const std::string &expression, double t);
```

**Source:** SDD §6.1 (`tdse.field.expression`, "same rule as `potential.function`"); config-schema-design doc's `tdse.field` block.

**Touches:** `mu::Parser`, already linked via `muparser::muparser` (`TISE/CMakeLists.txt:22–25`) — zero new dependency.

**Done when:** reproduces the config-schema-design doc's own worked example, `"0.05 * exp(-((t-50)^2)/(2*10^2)) * cos(0.2*t)"`, at hand-checked points (peak $\mathcal E(50)=0.05$; decays far from $t=50$); malformed expressions throw a clear error, mirroring `evaluateFunction`'s existing throw-on-uncovered-domain behavior.

### B4. Gauge dispatch → coupling operator

**What:** dispatch by gauge string to B1 (length) or B2 (velocity), returning the dense Hermitian coupling-operator matrix Engineer A's A1 diagonalizes:

```cpp
Eigen::MatrixXcd buildCouplingOperatorEigenbasis(const bspline::BSpline &bs,
                                                  const tise::EigenResult &er,
                                                  const std::string &gauge, // "length"|"velocity"
                                                  double hbar = 1.0);
```

**Source:** SDD Figure 8 ($\hat O=\hat x$ length / $\hat p_x$ velocity gauge); SDD §6.1 `tdse.operator.gauge`; `docs/adr/0006-tdse-propagator-choice.md` (Consequences — this task builds only the fixed-matrix half of what a Crank-Nicolson-family implementation would have needed; the closure-producing half is deferred with the CN family, not built).

**Touches:** B1/B2's outputs.

**Done when:** dispatch delegates correctly to B1/B2 by gauge string; **and**, as the actual integration checkpoint, Engineer A's A1 (`diagonalizeCouplingOperator`) consumes this output with zero code changes on A1's side.

### B5. Generalized initial-state construction

**What:** extend the Gaussian prototype to the full `{position, momentum, width}` config triple — `tevol::gaussianWavepacket(x, r0, mOmega, hbar)` currently has **no momentum parameter at all** (confirmed directly against `time_evolution.hpp:20`) despite `config.yaml`'s `tdse.initial_state` schema requiring one for `type: gaussian` — plus the trivial `eigenstate` case:

$$\psi_G(x,0) = (2\pi\sigma^2)^{-1/4}\exp\!\left(-\frac{(x-x_0)^2}{4\sigma^2}\right)e^{ik_0x}$$

integrated as separate real/imaginary parts via `bs.integral` (real-valued quadrature) into a genuinely complex `Eigen::VectorXcd` overlap, then projected to the eigenbasis (reuse `tevol::projectToEigenBasis`, generalized to accept a complex overlap vector — currently `Eigen::VectorXd`-only).

```cpp
Eigen::VectorXcd buildInitialState(const bspline::BSpline &bs,
                                    const tise::EigenResult &er,
                                    const std::string &type,       // "eigenstate"|"gaussian"
                                    int eigenstateIndex,
                                    double position, double momentum, double width,
                                    double hbar = 1.0);
```

**Source:** SDD §6.1 `tdse.initial_state`; SDD Figure 8 `INIT` node.

**Note for Engineer A:** this is the only place `alpha0` is ever built in production code — Engineer A never constructs one, only synthetic ones in A's own tests.

**Done when:** `eigenstate` case produces an exact one-hot vector; `gaussian` case with `momentum=0` reduces to the *existing, already-tested* `computeGaussianOverlaps`/`gaussianWavepacket` output after projection (a regression-consistency check against the untouched prototype); norm² is within a documented tolerance of 1.0 for parameter choices well inside the box.

### B6. Gauge-equivalence validation

**What:** the commutator-derived identity relating length- and velocity-gauge matrix elements between exact $H_0$ eigenstates (from $[\hat x,\hat H_0]=i\hbar\hat p_x/m$):

$$\langle m|\hat p_x|n\rangle = \frac{im(E_m-E_n)}{\hbar}\langle m|\hat x|n\rangle$$

checked numerically between B1's and B2's real, rotated operators for the same solved system. **This check is sign-convention-robust where B1/B2's individual textbook comparisons are not** — both operators are built from the *same* eigenvectors, so any arbitrary overall sign LAPACK assigned to a given eigenvector cancels consistently on both sides of the identity.

**Source:** standard QM commutator identity; SDD §9.3 ("agreement with known physics, not merely 'tests pass'"); Bachau et al. 2001 as the project's authoritative reference (cite per SDD §9.3's precedent — the in-repo PDF has not been directly re-verified against this identity's exact notation, so treat the citation as directional, not a confirmed page reference).

**Done when:** the identity holds to numerical tolerance across several $(m,n)$ pairs for the same real solved harmonic-oscillator system used in B1/B2.

---

## 4. Cleanup Tasks

Small, real gaps surfaced during research. Pick up whichever fits naturally alongside the adjacent workstream — none require interface/CLI work.

1. **Update `1D-QM-Playground/README.md`.** It doesn't mention `controller.py`, `config.yaml`, or `analysis.py` at all, despite all three already existing at the repo root — worse than merely stale, since it doesn't describe a pipeline that already exists on disk. (Distinct from `TISE/README.md`'s staleness, already tracked in `tise-task-breakdown.md`'s own cleanup list.)
2. **Flag (do not fix) `time_evolution.cpp`'s output naming.** `runTimeEvolution` writes `outputDir/Timestep_NNN` (`time_evolution.cpp:145–146`), not the documented `snapshot_NNNNN.dat` (SDD §6.3). A5 implements the correct convention as new functions; the old demo path is left alone — mirroring TISE cleanup #4's "flag don't fix" precedent for the same reason (fixing it means touching the CLI/output-dir contract, which is Phase 5, out of scope here).
3. **`TISE/CMakeLists.txt` hardcodes `/usr/include/eigen3` twice** (lines 68 and 87) instead of `find_package(Eigen3 REQUIRED)`. Not fixed here, but the new `tdse_lib` target (§5) should not repeat the pattern.
4. **Cross-reference (don't re-litigate) the unconditional time-evolution call.** `main.cpp:121–133` still unconditionally calls `tevol::runTimeEvolution` after every TISE solve, no `run_tdse` gating — already flagged in TISE's own cleanup list (#4), still open a month later. Independently corroborated: `controller.py` explicitly rejects `run.run_tdse: true` today with "not yet supported (Phase 5)."
5. **`BSpline.hpp`'s `integral()` docstring claims an `r² dr` weight** the actual `driverIntegral` implementation doesn't apply — a stale radial-Fortran-heritage comment. Worth a one-line fix so B1/B2 don't second-guess a nonexistent factor.

---

## 5. Shared Integration Checklist

- **File layout.** New files only, listed in §2/§3; `time_evolution.hpp`/`.cpp` stays untouched by both engineers — a stronger conflict-avoidance guarantee than TISE had for the shared `tise.hpp`/`.cpp` (applying that lesson learned up front rather than discovering it mid-implementation).
- **The two handoffs.** `buildCouplingOperatorEigenbasis` (B4 → A1, a fixed matrix diagonalized once) and `evaluateField` (B3 → A3, called fresh every step) are two plain functions, not a fused closure — see `docs/adr/0006-tdse-propagator-choice.md`'s Consequences for why. **Agree explicitly:** the `+i` sign in $e^{+i\theta\hat O_\text{eig}}$ is Engineer A's responsibility, applied in A1's `applyExponential` — Engineer B's `Oeig` and `evaluateField` are handed over sign-neutral (no minus sign baked in on B's side, unlike the earlier closure-based design this document originally considered). Sign-convention mismatches are the classic silent propagator bug; this is called out in both A1 and A3 above.
- **`CMakeLists.txt` (shared touch point).** One new target, `tdse_lib`, sourcing all four new `.cpp` files (`propagator.cpp` from Engineer A; `coupling_operator.cpp`, `field.cpp`, `initial_state.cpp` from Engineer B — each engineer adds their own lines, low conflict risk), linked against `bspline_lib`, `tise_lib`, Eigen, `muparser::muparser`. Use `find_package(Eigen3 REQUIRED)`/`Eigen3::Eigen` properly rather than repeating the hardcoded-path pattern (§4 cleanup #3), even though the *existing* targets aren't touched this phase. Corresponding `tests/CMakeLists.txt` entries follow the existing `add_executable`/`target_link_libraries`/`add_test` pattern.
- **`hbar`/mass convention.** All new functions accept `hbar`/mass as explicit plain parameters defaulting to `1.0`, matching the existing `tevol::` functions' own convention — never hardcode atomic units inside a new function body.
- **Shared small-basis test fixture.** Both engineers should develop and test against the same modest basis size already used in `test_time_evolution.cpp`, so independently-written tests don't silently diverge in numerical assumptions. Not a blocking dependency — a convention, mirroring TISE §5's "grid stability" soft-coordination bullet.
- **Do the norm cross-check (A4) early.** It's the strongest single validation of the whole architecture (confirms "S is identity in the eigenbasis" end-to-end) and only needs A4 plus either engineer's initial-state construction — worth running as soon as both exist, before either side builds much further on top of the assumption.
- **Late-stage joint physics benchmark — needs both engineers' real code, belongs to neither alone.** Two-level Rabi flopping using *real* components end to end: a configured potential with well-separated lowest two bound states (e.g. a finite square well), the real dipole element $d_{12}=\langle\phi_1|\hat x|\phi_2\rangle$ from B1, real eigenvalues from an actual `solveTISE` run, a near-resonant `field.expression` at $\omega_{12}=(E_2-E_1)/\hbar$, propagated through A3's real `propagateDriven`. Compare the observed flopping period against the RWA prediction $T_\text{Rabi}=2\pi\hbar/(|d_{12}|\mathcal E_0)$. Distinct from A1's synthetic-2-level unit test (made-up numbers, constant coupling, validates the exponential-apply primitive in isolation) — this is the real end-to-end integration checkpoint.

---

## 6. Explicitly Out of Scope

- **SDD §10.2 Phases 5–7** — Controller↔TDSE interface (the `tdse_solver` CLI, `--tise-dir`/`--output-dir` wiring), TDSE↔Analysis interface, Controller↔Analysis extended. Confirmed still genuinely unbuilt: no `tdse_solver_main.cpp` exists anywhere (unlike TISE's Phase-1 `tise_solver_main.cpp` stub), and `controller.py` explicitly rejects `run_tdse: true` with "not yet supported (Phase 5)."
- **`data/tdse/warnings.json` sidecar persistence** — A6 implements detection/warning logic only; the JSON-file mechanism is Phase-5 interface work, mirroring how `data/tise/warnings.json` was established as part of TISE's Phase-1 contract, not its Phase-4 numerics.
- **REQ-F-060 Analysis-module work** (`analysis.py`) — bound-state populations, asymptotic distribution, expectation values, interval probability. Confirmed `analysis.py` today only checks `--tdse-dir` existence without reading its contents.
- **Continuum-state initial conditions** — `tdse.initial_state.type` is `eigenstate | gaussian` only; no `continuum` option exists in the schema, Figure 8, or Figure 9. (Also currently blocked transitively: TISE's own continuum-construction work, `tise-task-breakdown.md`'s B1–B5, is unmerged.)
- **FEDVR / Krylov / Chebyshev propagators** — deferred per **ADR-0001**; their attractiveness is explicitly tied to FEDVR's diagonal-$V$/identity-$S$ properties, which the current (even eigenbasis-rotated) representation doesn't have.
- **Crank-Nicolson-family propagators (dense and banded) and the split-operator/CN hybrid** — deferred/rejected per **ADR-0006** (`docs/adr/0006-tdse-propagator-choice.md`), which adopts the symmetric split-operator method instead; see that ADR for the full comparison, the two CN-equivalence proofs, and the banded-CN revisit trigger.
- **The archived van Dijk/Chebyshev-expanded propagator design** (`docs/TDSE-original-design/`) — explicitly superseded; kept only for archival traceability, not a candidate to revive without a deliberate new decision superseding ADR-0006.
- **Norm-drift corrective action** (renormalization, adaptive step-size control) — SDD §5.3.4 requires monitoring/warning only; A6 is detection-only.
- **Generalized (non-default) B-spline drop-set support in TDSE** — every task in §3 assumes the default `{1, nBSplines}` convention (stated explicitly at the top of §3); supporting TISE's generalized `dropSet` (`tise.hpp`'s A4b work) here is a genuine open follow-up, flagged not attempted, mirroring how TISE's own doc flagged Coulomb-tail matching as a follow-up rather than attempting it.
- **Higher-order/adaptive-step propagators** — the adopted split-operator method's fixed-step, second-order scheme is the default; a natural future upgrade the two-handoff design doesn't foreclose, not attempted now.
