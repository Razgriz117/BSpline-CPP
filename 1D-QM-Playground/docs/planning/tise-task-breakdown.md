# TISE Solver — Task Breakdown for Two Engineers

**Date:** 2026-07-26
**Status:** Draft — ready for two-person parallel execution
**Scope:** SDD §5.2 ("TISE Solver") internal numerics only — i.e. SDD §10.2 **Phase 4** ("TISE Solver implementation"). This document does **not** cover Phases 1–3 (Controller↔TISE, TISE↔Analysis, Controller↔Analysis interface contracts — CLI flags, `yaml-cpp` config.yaml ingestion, the formal `data/tise/` output-directory contract) or Phases 5–8 (TDSE). Those remain separately-scoped future work.

---

## 1. Context & Scope

### Why this document exists

The SDD (`docs/SDD.md` §5.2.3, Figure 6, Figure 7) and the source worksheet `PHY5606_F25_ContinuumEigenstates.pdf` together fully specify the TISE solver's remaining numerics, but as prose/diagrams they don't divide into two independent, parallelizable workstreams. This document does that division so two engineers can work simultaneously with minimal blocking dependencies.

### What's already done (do not re-implement)

Confirmed by direct code read of `TISE/tise.cpp`/`tise.hpp` and exercised by `TISE/tests/test_tise.cpp`:

- `fillBandedMatrices` (`tise.cpp:113-155`) already assembles the banded $\mathbf{H}$/$\mathbf{S}$ matrices for a **general, config-shaped piecewise potential** — `evaluateFunction`/`inInterval` (`tise.cpp:45-111`) evaluate a `map<domain-string, muparser-expression-string>`, matching the `potential` DSL described in the config schema spec and SDD §6.1 (Figure 5). This is decoupled from the old hydrogenic-only `radialPotential` — it is **not** hardcoded to Coulomb.
- `solveGeneralizedEigenproblem` (`tise.cpp:157-193`) solves $\mathbf{Hc}=E\mathbf{Sc}$ via LAPACK `DSBGV` for any potential passed to `fillBandedMatrices`.
- `bspline::BSpline` (`BSpline.hpp`/`.cpp`) provides `eval` and `integral` (Gauss-Legendre), which both workstreams below reuse directly — no new quadrature or evaluation machinery is needed.

Note: `docs/SDD.md:587` and `TISE/README.md` currently describe `tise.cpp` as still hydrogenic-only. That description is **stale** relative to the actual code — see the cleanup tasks in §4.

### What's genuinely missing

Confirmed by exhaustive grep (zero hits for "continuum", "phase", "B_N", "boundary coupling" anywhere under `TISE/`): everything in Figure 6's right-hand branch (`CLASSIFY` → `CONTEN` → `PRECOMP` → `LOOP` → `BUILDPSI` → `MATCH`) and all of Figure 7 (the boundary-condition asymptote classifier). This is the actual content of the two workloads below.

### How the split works

The split follows Figure 6's fork out of `CLASSIFY` — one engineer per branch:

```
SOLVE → CLASSIFY{E_n below threshold?}
          ├─ yes → BOUND ─────────────────────────┐   Engineer A
          └─ no  → CONTEN{continuum.enabled?}      │
                     └─ yes → PRECOMP → LOOP        │   Engineer B
                              → BUILDPSI → MATCH ───┘
```

**These two workloads are parallelizable, not sequential.** The continuum-construction recipe (Engineer B) sums over the *entire* confined eigenbasis $\{\phi_n, E_n\}_{n=1}^{N-2}$ — per the PDF, with no threshold filter — so Engineer B's work only depends on the existing `solveTISE`/`fillBandedMatrices` output, not on Engineer A's classification or node-placement work. See §5 for the one soft coordination point (grid stability).

---

## 2. Engineer A — Bound-State Generalization & Boundary Conditions

Figure 6's `BOUND` path plus all of Figure 7.

### A1. Boundary-condition asymptote classifier (REQ-F-030)

**What:** For each side of the spatial domain (`bspline.domain`), determine bounded vs. unbounded. For unbounded sides, evaluate the assembled potential's asymptote (as $x \to$ that boundary) and dispatch:
- **Case 1** — no finite asymptote (diverges/grows, e.g. $x^2$, $x$): hard Dirichlet wall at $R$; all states are discrete pseudostates.
- **Case 2** — known analytic asymptote (flat, or Coulomb $\sim 1/r$): bound states diagonalized with Dirichlet at $R$ as usual; continuum matched to the analytic asymptotic form (this is the hand-off point to Engineer B's B2/B3 for the flat-asymptote sub-case — see note below).
- **Case 3** — unknown/irregular asymptote (e.g. $1/r^{1.5}$): approximate $V'(x)$ as flat beyond $R$; match continuum to a shifted sine as in Case 2's flat branch; **emit a user-facing warning** that this introduces a discontinuity at $x=R$ and that continuum normalization is approximate.

**Source:** `docs/planning/architecture-06-20.md`, "Extra boundary conditions" (2026-07-03 stakeholder feedback — the three-domain-configuration table and Case 1/2/3 dispatch logic, quoted in full); SDD §5.2.3 Figure 7.

**Touches:** new function, e.g. `tise::classifyAsymptote(potential, domain, side) -> AsymptoteCase`, called before/alongside `fillBandedMatrices`/`solveTISE`.

**Note for Engineer B:** Cases 2 (flat sub-branch) and 3 both reduce to the *same* flat-asymptote matching formula that Engineer B implements (B2/B3). Case 2's Coulomb-tail sub-branch (matching to Coulomb functions instead of $\sin(kx+\delta)$) is **not yet derived by any source document** (SDD §5.2.3 says so explicitly) — do not attempt it as part of this task; flag it as an open follow-up (§6).

**Done when:** unit tests cover at least one potential per case (e.g. harmonic-oscillator-like for Case 1, a flat/step potential for Case 2's flat sub-branch, a $1/r^{1.5}$-tailed potential for Case 3) and assert the correct case is selected, plus the Case 3 warning is emitted (to `stderr`, per SDD §8's warning taxonomy — non-fatal, does not affect exit code).

### A2. Ionization-threshold determination + bound/continuum classification (REQ-F-020)

**What:** Implement the `CLASSIFY{E_n below threshold?}` node. The ionization threshold is the potential's asymptotic value at the relevant (unbounded) boundary, from A1's classification. Label each diagonalized $E_n$ as bound (below threshold) or above-threshold pseudostate; report only bound states as "bound states" in output.

**Source:** SDD §5.2.3/Figure 6; `docs/planning/architecture-06-20.md`, "Number of bound states" (2026-07-03 feedback — "bound-state count is output, not input").

**Touches:** a small classification pass over `EigenResult::values`, likely returning a `std::vector<bool>` or count of bound states — implement as a standalone function (e.g. `tise::classifyBoundStates(EigenResult, threshold)`) rather than mutating `EigenResult` itself, to avoid the merge-conflict risk noted in §5.

**Done when:** for a known finite-bound-state potential (e.g. a finite square well) and a known infinite-bound-state potential (e.g. Coulomb), the classifier correctly separates bound vs. above-threshold count; regression-tested against the existing hydrogenic analytic energies already used in `test_tise.cpp`.

### A3. Well-containment diagnostic

**What:** For each reported bound state, compute $\psi'(x_\text{max})$. Nonzero (beyond numerical tolerance) flags the state as "not well-contained" — its energy/wavefunction may be inaccurate.

**Source:** `docs/planning/architecture-06-20.md`, "Number of bound states" — "check the derivative of the eigenstate at the box boundary... a red flag indicating the state is not well-contained"; SDD §5.2.3 Figure 6 (`BOUND` node), §6.4, §8 (physics-warning taxonomy).

**Touches:** reuses `bs.eval` with a derivative order argument (already supported by `BSpline::eval`/`integral`'s derivative-order parameters, per `docs/planning/bsplines.md`); attach the flag alongside each bound state's output.

**Done when:** a deliberately too-small box (states "colliding" with the wall) produces flagged states; a comfortably large box does not.

### A4. Strategic node placement (REQ-F-050)

**What:** Detect potential structure from the parsed `potential` pieces and insert/remove knots accordingly, per the required-treatment table:

| Potential type | Required treatment |
|---|---|
| Delta potential $\delta(x-x_0)$ | Pile up degenerate knots at $x_0$ (discontinuity in $\psi'$) |
| Potential step | Pile up degenerate knots at the step location (discontinuity in $\psi''$) |
| Stitched potentials, continuous derivative | Knot degeneracy at the join point (discontinuity in $\psi'''$) |
| Singular potentials (e.g. $1/r$) | Remove B-splines at the singular point (enforce regularity of $\psi$) |

**Explicitly out of scope:** WKB-proportional node *density* is deferred per **ADR-0002** — only *strategic* placement (the table above) is in scope. Do not implement a reference-energy-driven density scheme.

**Source:** `docs/planning/architecture-06-20.md`, "Collocation scheme" (2026-07-03 feedback); `docs/adr/0002-defer-wkb-collocation.md`.

**Touches:** the grid-construction step upstream of `bspline::BSpline::init` (currently `tise::buildUniformRadialGrid`) — extend or wrap it to accept strategic knots derived from the `potential` map's domain boundaries/expression forms.

**Done when:** a potential with a domain-boundary discontinuity (e.g. the existing "particle in a box with a rectangular barrier" example from the config schema spec) produces a grid with degenerate knots at the boundary, and eigenvalue accuracy improves (smaller `eigenvalueError`-style residual) relative to a purely uniform grid at the same node count.

### A5. `E_acc` continuum-accuracy warning (REQ-F-040, warning half)

**What:** Compute the basis accuracy ceiling from node spacing, $E_\text{acc} \gtrsim \dfrac{\pi^2}{2m\,\Delta x_\text{node}^2}$, and warn if the requested continuum `E_max` (accept as a plain function parameter for now — no config schema wiring, that's out of scope) exceeds it.

**Source:** `docs/planning/architecture-06-20.md`, "Continuum range" (2026-07-03 feedback — full derivation quoted); SDD §5.2.3/§6.4/§8.

**Touches:** a small standalone function, e.g. `tise::computeEAcc(nodeSpacing, mass)`, plus a warning emission at the call site once Engineer B's energy-grid loop (B2) exists — coordinate on the call site, but the function itself has no dependency on Engineer B's code.

**Done when:** unit test asserts $E_\text{acc}$ matches the closed-form formula for a known node spacing, and that the warning fires (to `stderr`) exactly when `E_max > E_acc`.

---

## 3. Engineer B — Continuum-State Construction

Figure 6's `CONTEN` → `PRECOMP` → `LOOP` → `BUILDPSI` → `MATCH` path, per SDD §5.2.3 and `PHY5606_F25_ContinuumEigenstates.pdf` (the two are a verbatim match — SDD §5.2.3 is a direct transcription of the PDF's algorithm, generalized only in the energy-grid endpoint, see B2 below).

### B1. Precompute boundary-coupling elements

**What:** Compute $\langle\phi_n|H|B_N\rangle$ and $\langle\phi_n|B_N\rangle$ for every confined eigenstate $n=1,\dots,N-2$, where $B_N$ is the *last* B-spline (normally dropped to enforce $\psi(R)=0$ for bound states, but not dropped for genuine continuum solutions, which need not vanish at $R$). Compute once, reusing the same $\mathbf{H}$ and basis already assembled for bound states.

**Source:** SDD §5.2.3 (verbatim): *"The only energy-independent quantities this requires precomputing... are $\langle\phi_n|H|B_N\rangle$ and $\langle\phi_n|B_N\rangle$ for each $n$."* PDF, "Now, we need to figure out how to find the generalized eigenfunctions for arbitrary values of the energy $E$."

**Touches:** new function, e.g. `tise::precomputeBoundaryCoupling(bs, H, EigenResult) -> BoundaryCoupling{ std::vector<Real> phiHB, phiB }`. **Implement as a standalone companion struct, not as new fields on `EigenResult`** — both engineers touch `tise.hpp`/`tise.cpp`, and keeping this as an additive, separate struct avoids merge conflicts with Engineer A's classification additions (SDD §6.2 flags this struct extension as needed but doesn't mandate where it lives).

**Done when:** unit test validates $\langle\phi_n|B_N\rangle$ and $\langle\phi_n|H|B_N\rangle$ against a directly-computed (brute-force `bs.integral` call, no shortcuts) reference for a small basis size.

### B2. Energy-grid loop and closed-form coefficients

**What:** Loop over a uniform energy grid $\varepsilon_i$ and, for each energy $E$, compute:

$$c_n = \frac{\langle\phi_n|H|B_N\rangle - E\langle\phi_n|B_N\rangle}{E - E_n}, \qquad |\bar\psi_E\rangle = \sum_{n=1}^{N-2}|\phi_n\rangle c_n + |B_N\rangle$$

This is an algebraic evaluation per energy (no re-diagonalization).

**Grid endpoint note:** the PDF's grid is $\varepsilon_i = \dfrac{E_\text{max}}{N_E}i$ (implicitly $[0, E_\text{max}]$, sampling starts at $i=1$ to avoid the $k=0$ singularity in the matching formulas of B3). SDD §5.2.3 generalizes this to a user-supplied `[E_threshold, E_max]` window — implement the generalized form (e.g. $\varepsilon_i = E_\text{threshold} + \dfrac{E_\text{max}-E_\text{threshold}}{N_E}i$), accepting `E_threshold`, `E_max`, `N_E` as plain function parameters (again, no config schema wiring — out of scope).

**Source:** SDD §5.2.3; PDF page 1, the $c_n$ derivation from $\langle\phi_n|(E-H)|\psi_E\rangle=0$.

**Touches:** new function, e.g. `tise::buildContinuumState(EigenResult, BoundaryCoupling, bs, E) -> std::vector<Real>` (B-spline coefficient vector for $\bar\psi_E$, including the $B_N$ term), plus the grid-loop driver.

**Done when:** unit test confirms $\langle\phi_n|(E-H)|\bar\psi_E\rangle \approx 0$ for every confined $n$ at a sampled energy (the defining property used to derive $c_n$ in the first place) — a strong, direct correctness check independent of the matching step in B3.

### B3. Asymptotic matching at $x=R$

**What:** Evaluate $\bar\psi_E(R)$ and $\bar\psi_E'(R)$ (via `bs.eval` with derivative order 0 and 1), then extract:

$$A_E = \sqrt{\frac{2/\pi}{k\,\bar\psi_E(R)^2 + \bar\psi_E'(R)^2/k}}, \qquad \delta(E) = \arctan\!\left[\frac{k\,\bar\psi_E(R)}{\bar\psi_E'(R)}\right] - kR, \qquad k=\sqrt{2E}$$

and the numerically stable phase-shift derivative (avoids differentiating $\arctan$'s branch cuts directly):

$$\frac{d\delta}{dE} = \frac{1}{2\cos(2\delta)}\frac{d\sin(2\delta)}{dE}$$

**Source:** SDD §5.2.3 (verbatim, matches PDF exactly, including the derivative formula).

**Scope boundary:** this is the **flat-asymptote case only** — the general recipe behind Figure 7's Case 2 (flat sub-branch) and Case 3. The Coulomb-tail sub-branch of Case 2 is explicitly out of scope (§6).

**Touches:** new function, e.g. `tise::matchAsymptotic(psiR, psiPrimeR, E, R) -> {A_E, delta, dDeltaDE}`.

**Done when:** for a simple analytically-known flat-asymptote test case (e.g. a finite square well beyond $R$, or free-particle-like potential), the extracted $\delta(E)$ matches the known analytic phase shift within tolerance.

### B4. Output formatting

**What:** Produce the in-memory content (or a simple write-to-stream function, following the existing `tise::writeEigenstate` pattern at `tise.cpp:218-235`) matching the documented column formats:

| Content | Format |
|---|---|
| `phase_shifts.dat` | 3-col: $\varepsilon_i$, $\delta(\varepsilon_i)$, $d\delta/dE$ |
| `continuum_state_NNN.dat` | 2-col: $x$, $\psi_{\varepsilon_i}(x)$, tabulated on $x_i = \dfrac{R}{N_x-1}(i-1)$ |

**Explicitly not in scope:** wiring this through `--output-dir`/`data/tise/` — that's the (excluded) Phase 1/2 interface contract. Write plain functions callable from tests, following the same shape as the existing `writeEigenstate(std::ostream&, ...)` signature.

**Source:** SDD §6.3 (verbatim file-format table); §5.2.3.

**Done when:** output matches the documented column count/order/precision conventions already used by `writeEigenstate` (`std::scientific`, `setprecision(16)`).

### B5. Validation against the Bachau reference

**What:** Benchmark computed phase shifts $\delta(E)$ against H. Bachau et al., *Rep. Prog. Phys.* **64**, 1815 (2001) (`H_Bachau_2001_Rep._Prog._Phys._64_1815.pdf`), per SDD §9.3's stated validation strategy ("agreement with known physics, not merely 'tests pass'").

**Source:** SDD §9.3.

**Done when:** at least one worked example from Bachau (or a simple square-well/step potential with an independently-derivable analytic phase shift) is reproduced within a stated numerical tolerance, documented alongside the test.

**Explicitly flagged as a follow-up, not part of this task:** general Coulomb-tail matching (Case 2's non-flat sub-branch) remains undocumented by any source material — SDD §5.2.3 says so outright. Do not attempt it here; record it as an open item (§6) for whoever picks up boundary-condition work beyond this phase.

---

## 4. Cleanup Tasks

Small, real gaps surfaced during research. Pick up whichever fits naturally alongside the adjacent workstream — none require interface/CLI work.

1. **Fix `TISE/make_and_run.sh`.** It currently invokes `./build/H-BoundStates` with zero arguments; `main.cpp:60` unconditionally dereferences `argv[1]` via `parsePiecewise(argv[1])`, so this will crash. Update the script to pass a valid JSON potential array argument matching the actual CLI contract (e.g. `'[{"domain": "(0, 100]", "function": "-1/x + 1/x^2"}]'`).
2. **Update `TISE/README.md`.** It still describes the project as hydrogenic-only with angular momentum set via a `constexpr int L` — it doesn't mention the muparser-based general piecewise-potential mechanism that's actually implemented (see §1 above). Update it to describe the current JSON-argv/muparser mechanism.
3. **Reconcile the expression-parser decision record.** `docs/SDD.md` §11.1 and `docs/planning/resources.md` both describe the expression-parser choice (FunctionParser vs. NFParam) as "not yet chosen definitively" — but `TISE/CMakeLists.txt:16-18` has already committed to **muparser** (via `pkg_check_modules(MUPARSER REQUIRED IMPORTED_TARGET muparser)`), which isn't even one of the two shortlisted candidates. Update SDD §11.1 / `resources.md` to record muparser as the adopted choice, or open a short ADR if there's genuine appetite to revisit.
4. **Flag (do not fix) the unconditional time-evolution call.** `main.cpp:121-133` always calls `tevol::runTimeEvolution(...)` after solving the TISE problem, regardless of whether a bound/continuum-only run was intended — there is no `run_tise`/`run_tdse` gating, since no config-flag plumbing exists in `main.cpp` at all yet. A proper fix requires the (out-of-scope) config-driven run-gating from the interface phases. Document this as a known limitation; do not attempt to fix it in this phase.

---

## 5. Shared Integration Checklist

- **`tise.hpp`/`tise.cpp` is a shared file.** Both engineers add new functions here. To minimize merge conflicts: Engineer A's classification/diagnostic work (A1–A3) should be additive standalone functions taking `EigenResult` as input, not modifications to the `EigenResult` struct itself; Engineer B's boundary-coupling data (B1) should live in its own companion struct (e.g. `BoundaryCoupling`), not as new `EigenResult` fields. Neither engineer should need to touch the other's new functions directly.
- **Grid stability.** Engineer A's strategic node placement (A4) changes the grid `solveTISE` builds. Engineer B's continuum construction (B1–B3) is grid-agnostic in principle (it just consumes whatever `EigenResult`/basis it's given), but should be re-validated against the final grid once A4 lands, since node placement affects numerical accuracy of the boundary-coupling integrals. This is not a blocking dependency — Engineer B should develop and validate against the current uniform grid from day one.
- **Validate via the existing GTest pattern, not via `main.cpp`/CLI.** Follow `tests/test_tise.cpp`'s existing style (hand-built `map<string,string>` potentials, direct calls into `tise::` namespace functions) for all new functions in both workloads. Neither workload requires touching `main.cpp` except optionally, for manual/demo purposes.
- **Case 2/Case 3 hand-off (A1 ↔ B3).** Engineer A's asymptote classifier determines *which* matching formula applies; Engineer B's B3 implements the flat-asymptote formula that both Case 2 (flat sub-branch) and Case 3 use. Confirm at integration time that A1's classifier calls into B3's matching function for these two cases, and that Case 1 (hard wall, no continuum) correctly skips the continuum branch entirely.

---

## 6. Explicitly Out of Scope

- **SDD §10.2 Phases 1–3** — Controller↔TISE interface (CLI `--config`/`--output-dir`, `yaml-cpp` config.yaml ingestion), TISE↔Analysis interface, Controller↔Analysis interface. The full `data/tise/` output-directory contract (as opposed to the in-memory/stream-based output functions in A/B above) belongs to this excluded phase.
- **SDD §10.2 Phases 5–8** — TDSE solver work of any kind.
- **FEDVR** as an alternative basis — deferred, **ADR-0001**.
- **WKB-proportional node density** — deferred, **ADR-0002** (strategic placement per A4 is in scope; the density scheme is not).
- **Coulomb-tail continuum matching** (Figure 7 Case 2's non-flat sub-branch) — not yet derived by any source document per SDD §5.2.3; flagged in both A1 and B5 above as a follow-up for whoever picks up boundary-condition work next.
- **CAP / outgoing-wave (Siegert) boundary conditions / exterior complex scaling (ECS)** — still a genuinely open question per SDD Appendix B ("no REQ, no ADR — the next time this is discussed with the stakeholder, it should be promoted to either a new REQ or a new ADR"). Not part of this task breakdown either way.
