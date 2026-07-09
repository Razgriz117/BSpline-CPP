# Diagram Reference — 1D-QM-Playground

*Generated from documentation only, not from the C++/Python source, since the
implementation is expected to change substantially while the documented
design is comparatively stable. Every diagram below cites the doc(s) it was
drawn from. Diagrams marked "reused verbatim" are copied unmodified from an
existing doc; everything else is newly authored here to visualize prose,
tables, or pseudocode that had no diagram of its own.*

*The docs disagree on the TDSE propagation method (see Section D). Only the
current design is diagrammed here; the superseded Chebyshev/van Dijk propagator
diagram lives at `docs/TDSE-original-design/tdse-propagator-diagram.md`, alongside
the archived docs describing it.*

## Contents

- **A. System Architecture & Data Flow** — 1–3
- **B. Configuration** — 4–5
- **C. TISE Solver** — 6–7
- **D. TDSE Solver** — 8–9
- **E. B-Spline Numerics** — 10–11
- **F. Process & History** — 12–14

---

## Section A — System Architecture & Data Flow

## 1. System Component Architecture (`flowchart`)

*Source: `docs/planning/architecture-06-26.md`, §1 (current). Reused verbatim.*

```mermaid
flowchart TD
    CTRL["**Controller**\n(Python · controller.py)\nParse config.yaml → orchestrate pipeline"]

    TISE["**TISE Solver**\n(C++ · tise_solver)"]
    TDSE["**TDSE Solver**\n(C++ · tdse_solver)"]
    ANALYSIS["**Analysis**\n(Python · analysis.py)"]

    TISE_STORE[("**data/tise/**\neigenvalues · eigenvectors\nH/S matrices\ncontinuum states · phase shifts")]
    TDSE_STORE[("**data/tdse/**\nwavefunction snapshots\nnorm · energy · P(t)")]

    CTRL -->|"subprocess call"| TISE
    CTRL -->|"subprocess call"| TDSE
    CTRL -->|"subprocess call"| ANALYSIS

    TISE -->|"writes"| TISE_STORE
    TISE_STORE -->|"reads"| TDSE
    TISE_STORE -->|"reads"| ANALYSIS

    TDSE -->|"writes"| TDSE_STORE
    TDSE_STORE -->|"reads"| ANALYSIS
```

The system is four loosely-coupled components orchestrated by a central Python **Controller**. The C++ solvers are independent executables; **Analysis** is a Python script; all persistent data moves through the `data/` directory.

---

## 2. Orchestration Sequence (`sequenceDiagram`)

*Source: controller pseudocode in `docs/planning/architecture-06-26.md` §8, and the `run.*` flags in `config.yaml`.*

```mermaid
sequenceDiagram
    actor User
    participant Controller as Controller (controller.py)
    participant TISE as TISE Solver (tise_solver)
    participant TDSE as TDSE Solver (tdse_solver)
    participant Analysis as Analysis (analysis.py)
    participant Storage as data/ (tise/, tdse/)

    User->>Controller: run(config.yaml)
    Controller->>Controller: parse and validate config.yaml

    alt run.run_tise == true
        Controller->>TISE: subprocess: tise_solver --config --output-dir data/tise
        TISE->>Storage: write eigenvalues, eigenvectors, H, S, continuum states, phase shifts
        TISE-->>Controller: exit code
    end

    alt run.run_tdse == true
        Controller->>TDSE: subprocess: tdse_solver --config --tise-dir data/tise --output-dir data/tdse
        TDSE->>Storage: read H, S, eigenpairs
        TDSE->>Storage: write wavefunction snapshots, observables
        TDSE-->>Controller: exit code
    end

    alt run.run_analysis == true
        Controller->>Analysis: subprocess: analysis.py --config --tise-dir --tdse-dir
        Analysis->>Storage: read data/tise/, data/tdse/
        Analysis-->>Controller: exit code
    end

    Controller-->>User: run complete
```

Each stage is gated by a boolean in `config.yaml`'s `run` block; the Controller "does not implement any physics — it is purely an orchestration layer."

---

## 3. Data Artifact Map (`flowchart`)

*Source: Data Flow Summary (§3) and Local Storage Format (§6) tables in `docs/planning/architecture-06-26.md`.*

```mermaid
%%{init: {'flowchart': {'nodeSpacing': 45, 'rankSpacing': 90, 'curve': 'basis'}}}%%
flowchart LR
    subgraph TISE_OUT["Produced by TISE Solver"]
        direction TB
        F1["eigenvalues.dat\nindex, E_n"]
        F2["eigenvectors.dat\ncolumns = c_n coefficient vectors"]
        F3["hamiltonian.dat\nH matrix, banded"]
        F4["overlap.dat\nS matrix, banded"]
        F5["phase_shifts.dat\neps_i, delta(eps_i), d(delta)/dE"]
        F6["continuum_state_NNN.dat\nx, psi_eps(x) per energy"]
    end

    TDSE["TDSE Solver"]

    subgraph TDSE_OUT["Produced by TDSE Solver"]
        direction TB
        F7["snapshot_NNNNN.dat\nx, Re(psi), Im(psi) per time step"]
        F8["observables.dat\nt, norm, E_avg, P(t)"]
    end

    ANALYSIS["Analysis"]

    F1 --> TDSE
    F2 --> TDSE
    F3 --> TDSE
    F4 --> TDSE

    TDSE --> F7
    TDSE --> F8

    F1 --> ANALYSIS
    F2 --> ANALYSIS
    F3 --> ANALYSIS
    F4 --> ANALYSIS
    F5 --> ANALYSIS
    F6 --> ANALYSIS
    F7 --> ANALYSIS
    F8 --> ANALYSIS
```

The layout now follows the actual pipeline order left to right (TISE files → TDSE Solver → TDSE files → Analysis), with the TDSE Solver node explicitly writing `F7`/`F8` so that column has something to anchor to; wider rank/node spacing gives the long TISE-to-Analysis edges (for `eigenvalues`, `eigenvectors`, `H`, `S`, which per the Data Flow Summary table feed both TDSE and Analysis) room to route around the middle columns instead of crossing through them. Plain text is preferred initially "for transparency and ease of inspection with standard tools"; migration to HDF5 is a possible later change that would not alter this producer/consumer map.

---

## Section B — Configuration

## 4. config.yaml Schema Structure (`classDiagram`)

*Source: `docs/superpowers/specs/2026-06-28-config-yaml-schema-design.md` + inline comments in `config.yaml`.*

```mermaid
classDiagram
    class Config {
        +run: Run
        +physics: Physics
        +bspline: BSpline
        +tise: TISE
        +tdse: TDSE
        +analysis: Analysis
        +visualization: Visualization
    }
    class Run {
        +run_tise: bool
        +run_tdse: bool
        +run_analysis: bool
        +output_dir: string
    }
    class Physics {
        +mass: float
        +hbar: float
    }
    class BSpline {
        +n_nodes: int
        +order: int
        +domain: float
    }
    class PotentialPiece {
        +domain: string
        +function: string
    }
    class TISE {
        +n_pts_eigenstate: int
        +error_threshold: float
    }
    class Continuum {
        +enabled: bool
        +E_threshold: float
        +E_max: float
        +n_energies: int
        +n_pts: int
    }
    class TDSE {
        +dt: float
        +t_final: float
        +snapshot_interval: int
    }
    class InitialState {
        +type: string
        +index: int
        +position: float
        +momentum: float
        +width: float
    }
    class Operator {
        +gauge: string
    }
    class Field {
        +enabled: bool
        +expression: string
    }
    class Analysis {
        +bound_state_populations: bool
        +asymptotic_populations: bool
        +asymptotic_distribution: bool
    }
    class ExpectationValues {
        +x: bool
        +p: bool
        +T: bool
        +V: bool
        +H: bool
    }
    class IntervalProbability {
        +enabled: bool
        +intervals: float
    }
    class Visualization {
        +eigenstates: bool
        +phase_shifts: bool
        +time_evolution: bool
        +bound_state_populations: bool
        +asymptotic_populations: bool
        +asymptotic_distribution: bool
        +expectation_values: bool
    }

    Config "1" *-- "1" Run
    Config "1" *-- "1" Physics
    Config "1" *-- "1" BSpline
    Config "1" *-- "1..*" PotentialPiece : potential
    Config "1" *-- "1" TISE
    Config "1" *-- "1" TDSE
    Config "1" *-- "1" Analysis
    Config "1" *-- "1" Visualization
    TISE "1" *-- "1" Continuum : continuum
    TDSE "1" *-- "1" InitialState : initial_state
    TDSE "1" *-- "1" Operator : operator
    TDSE "1" *-- "0..1" Field : field
    Analysis "1" *-- "1" ExpectationValues : expectation_values
    Analysis "1" *-- "0..1" IntervalProbability : interval_probability
```

`bspline` is top-level (shared infrastructure for both solvers); `continuum` nests under `tise` because "continuum states are computed by the TISE solver — nesting expresses that dependency." `analysis` and `visualization` are deliberately separate top-level blocks: "conflating them mislabeled TISE/TDSE raw outputs as 'analysis'."

---

## 5. Potential Definition DSL (`flowchart`)

*Source: the `potential` block comment in `config.yaml` and its schema entry in the config-schema spec.*

```mermaid
flowchart TD
    CFG["config.yaml: potential (YAML list)"]
    ENTRY["List entry: string-encoded dict\n{'domain': ..., 'function': ...}"]
    PARSE["Parse string as dict\n(ast.literal_eval in Python,\nequivalent parser in C++)"]
    DOM["domain: interval notation\n[a,b]  (a,b)  [a,b)  (a,b]\nsupports inf, e.g. [0, inf)"]
    FN["function: expression in x\nliteral numeric constants only\nno references to other config fields"]
    TILE["Pieces tile the spatial domain\nwithout gaps or overlaps\n(behavior on gaps/overlaps undefined)"]
    VX["Assembled piecewise V(x)"]

    CFG --> ENTRY --> PARSE
    PARSE --> DOM
    PARSE --> FN
    DOM --> TILE
    FN --> TILE
    TILE --> VX
```

All constants, "including centrifugal terms, must be baked in" — there is no cross-referencing between config fields inside a potential expression.

---

## Section C — TISE Solver

## 6. TISE Solver Internal Flow (`flowchart`)

*Source: `docs/planning/architecture-06-26.md` §2.2, and Layer 2 + the bound-state stakeholder decision in `docs/planning/architecture-06-20.md`.*

```mermaid
flowchart TD
    START(["config.yaml: bspline, potential, tise.*"])
    BUILD["Build B-spline basis\non the spatial grid\n(bspline.n_nodes, order, domain)"]
    ASSEMBLE["Assemble banded H and S\n(Gauss-Legendre quadrature)"]
    SOLVE["Solve generalized eigenproblem\nH c = E S c   (LAPACK DSBGV)"]
    CLASSIFY{"E_n below\nionization threshold?"}
    BOUND["Bound state\nreport E_n, phi_n(x)\ncheck psi'(x_max): if != 0,\nflag as not well-contained"]
    CONTEN{"tise.continuum.enabled?"}
    LOOP["Loop energy grid\neps_i = E_max / N_E * i"]
    MATCH["Compute psi_bar_E,\nextract A_E, delta(E), d(delta)/dE"]
    WRITE["Write to data/tise/:\neigenvalues, eigenvectors, H, S,\ncontinuum states, phase_shifts"]

    START --> BUILD --> ASSEMBLE --> SOLVE --> CLASSIFY
    CLASSIFY -->|yes| BOUND --> WRITE
    CLASSIFY -->|no, above threshold| CONTEN
    CONTEN -->|yes| LOOP --> MATCH --> WRITE
    CONTEN -->|no| WRITE
```

Bound-state count is an **output**, not an input: "the solver always computes all sub-threshold states." The well-containment check (`psi'(x_max) != 0`) is the documented diagnostic for flagging states whose energy/wavefunction may be inaccurate.

---

## 7. Boundary Condition Decision Tree (`flowchart`)

*Source: "Extra boundary conditions" stakeholder-feedback subsection (2026-07-03) in `docs/planning/architecture-06-20.md`.*

```mermaid
flowchart TD
    SIDE{"For each domain side:\nbounded or unbounded?"}
    BOUNDED["Bounded side\nApply Dirichlet: psi = 0 at wall\n(only user-facing BC for bounded sides)"]
    UNBOUNDED["Unbounded side\nProgram evaluates the asymptote\nof the assembled potential"]
    CASE1{"What is the asymptote?"}
    C1["Case 1: No finite asymptote\n(e.g. x^2, x -- diverges)\nHard wall (Dirichlet) at box edge R\nAll states = discrete pseudostates"]
    C2["Case 2: Known analytic asymptote\n(flat, or Coulomb ~ 1/r)\nBound states: diagonalize w/ Dirichlet at R\nContinuum: match to Coulomb functions or\na shifted sine sin(kx+delta) at R"]
    C3["Case 3: Unknown/irregular asymptote\n(e.g. 1/r^1.5)\nApproximate V'(x) as flat beyond R\nMatch continuum to a shifted sine\nWarn user: discontinuity at x=R,\nnormalization is approximate"]

    SIDE -->|bounded| BOUNDED
    SIDE -->|unbounded| UNBOUNDED
    UNBOUNDED --> CASE1
    CASE1 -->|diverges or grows without bound| C1
    CASE1 -->|flat or Coulomb 1/r| C2
    CASE1 -->|any other asymptote| C3
```

Three domain configurations are supported end-to-end: finite box (both sides bounded), half-line (left bounded, right unbounded — the radial equation), and full line (both sides unbounded).

---

## Section D — TDSE Solver

The docs disagree on the TDSE propagation method. `docs/planning/architecture-06-26.md` (current, matches the live `config.yaml` schema) describes a general driven-Hamiltonian design using eigenstate expansion; an older copy of the same-named document (archived under `docs/TDSE-original-design/` and duplicated at the top-level `Quantum Mechanics BSplines C++ Project/architecture-06-26.md`) describes an explicit Chebyshev/van Dijk propagator that the config-schema spec's Design Decisions table explicitly says was removed (`tdse.chebyshev_order` — "propagator method is an implementation detail, not a user-facing config field"). This doc only diagrams the current design below; the superseded Chebyshev/van Dijk propagator diagram has been moved to `docs/TDSE-original-design/tdse-propagator-diagram.md`, alongside the archived docs it illustrates, until the team formally retires it for good.

## 8. TDSE Solver Internal Flow — Eigenstate Expansion `[CURRENT]` (`flowchart`)

*Source: `docs/planning/architecture-06-26.md` §2.3 (current) + `docs/superpowers/specs/2026-06-28-config-yaml-schema-design.md` §`tdse`.*

```mermaid
flowchart TD
    IN["Inputs: data/tise/ (E_n, c_n, S)\ntdse.initial_state, tdse.operator.gauge,\ntdse.field, dt, t_final, snapshot_interval"]
    INIT["Expand initial state in H0 eigenbasis:\nalpha_n(0) = (phi_n | psi(0))"]
    FIELDQ{"tdse.field.enabled?"}
    STATIC["Field-free special case:\nalpha_n(t) = alpha_n(0) * e^(-i E_n t / hbar)\n|alpha_n(t)|^2 conserved --\nno population transfer between H0 eigenstates\n(e.g. Rydberg wave-packet dephasing/revivals)"]
    DRIVEN["General driven case:\nH(t) = H0 + H_int(t)\nH_int(t) = -O_hat * field(t)\nO_hat = x_hat (length gauge) or\np_x_hat (velocity gauge)\nnumerically integrate forward in time\n(propagator scheme is an implementation\ndetail, not user-configurable)"]
    POP["Genuine population transfer between\nbound levels and into the continuum"]
    RECORD["At each snapshot: record\npsi(x,t), norm, E_avg(t)"]
    WRITE["Write to data/tdse/:\nsnapshot_NNNNN.dat, observables.dat"]

    IN --> INIT --> FIELDQ
    FIELDQ -->|false| STATIC --> RECORD
    FIELDQ -->|true| DRIVEN --> POP --> RECORD
    RECORD --> WRITE
```

This is the mode that makes bound-state population vs. time a *meaningful, non-constant* analysis output only when a field is applied — the field-free case conserves populations by construction.

---

## 9. TDSE Eigenstate Round-Trip — math notation (`flowchart`)

*Source: `docs/planning/bsplines.md`, "The Round-Trip: Enabling TDSE". Reused verbatim.*

```mermaid
flowchart LR
    d0["Initial state\nin B-spline coords\n$$\mathbf{d}_0$$"]
    a0["Eigenstate\ncoefficients\n$$\boldsymbol{\alpha}(0)$$"]
    at["Time-evolved\neigenstate coeffs\n$$\boldsymbol{\alpha}(t)$$"]
    dt["B-spline coeffs\nat time t\n$$\mathbf{d}(t)$$"]
    psi["Real-space\nwavefunction\n$$\Psi(r,t)$$"]

    d0 -->|"$$C^{-1}$$"| a0
    a0 -->|"$$e^{-iE_n t}$$\ncomponent-wise"| at
    at -->|"$$C$$"| dt
    dt -->|"evaluate\nB-splines"| psi
```

This is the concrete algorithm behind Diagram 8's field-free special case, as implemented (prototype) in `TISE/main.cpp`: `Phi_G = C_dagger * B_G` → `V_t = evolution.cwiseProduct(Phi_G)` → `CV_t = C * V_t` → `bspline.eval`. It closes only because the eigenstates live in the B-spline subspace by construction, making `C` a well-defined, invertible change-of-basis matrix.

---

## Section E — B-Spline Numerics

## 10. B-Spline Construction via de Boor Recursion (`flowchart`)

*Source: `docs/planning/bsplines.md`, "Recursion Tree for $B_i^3$". Reused verbatim.*

```mermaid
flowchart BT
    B1i["$$B_i^1$$"]
    B1i1["$$B_{i+1}^1$$"]
    B1i2["$$B_{i+2}^1$$"]
    B1i3["$$B_{i+3}^1$$"]

    B2i["$$B_i^2$$"]
    B2i1["$$B_{i+1}^2$$"]
    B2i2["$$B_{i+2}^2$$"]

    B3i["$$B_i^3$$"]

    B1i  --> B2i
    B1i1 --> B2i
    B1i1 --> B2i1
    B1i2 --> B2i1
    B1i2 --> B2i2
    B1i3 --> B2i2
    B2i  --> B3i
    B2i1 --> B3i
```

Each node is a weighted blend of the two nodes below it. At any evaluation point `x`, only `k` B-splines are nonzero simultaneously, so exactly `k` leaves of this tree contribute. The project currently uses `k = 12`.

---

## 11. FEDVR vs. B-Spline Operator Structure `[future / deferred]` (`flowchart`)

*Source: comparison table (§4.4) in `docs/planning/fedvr-exploration.md`. Status: "deferred. To be revisited once the B-spline code is working."*

```mermaid
flowchart LR
    subgraph BSP["B-spline basis (current)"]
        direction TB
        BS_S["S: banded, bandwidth k"]
        BS_V["V(x): banded, bandwidth k"]
        BS_T["T: banded, bandwidth k"]
        BS_H["H = T + V: banded"]
        BS_EIG["H c = E S c\ngeneralized eigenproblem"]
    end

    subgraph FEDVR["FEDVR basis (deferred future option)"]
        direction TB
        FD_S["S: identity I\northonormal DVR functions"]
        FD_V["V(x): diagonal\nV(x_alpha) on the diagonal"]
        FD_T["T: banded\nblock-tridiagonal"]
        FD_H["H = T + V: banded"]
        FD_EIG["H c = E c\nstandard eigenproblem"]
    end

    NOTE["FEDVR = B-spline basis in the limit where\nevery interior knot has multiplicity k-1\n(DVR transform of x_hat within each element)"]

    BSP -.->|"deferred migration path\nsee fedvr-exploration.md"| FEDVR
    NOTE -.-> FEDVR
```

Not part of the current architecture — included here as a documented future option, not a present component.

---

## Section F — Process & History

## 12. Simulation Run State Diagram (`stateDiagram-v2`)

*Source: `run.run_tise` / `run.run_tdse` / `run.run_analysis` flags in `config.yaml`, and the controller pseudocode in `docs/planning/architecture-06-26.md` §8.*

```mermaid
stateDiagram-v2
    [*] --> Idle
    Idle --> ParseConfig: run(config.yaml)
    ParseConfig --> TISECheck

    state TISECheck <<choice>>
    TISECheck --> TISERunning: run.run_tise == true
    TISECheck --> TDSECheck: run.run_tise == false
    TISERunning --> TDSECheck: writes data/tise/

    state TDSECheck <<choice>>
    TDSECheck --> TDSERunning: run.run_tdse == true
    TDSECheck --> AnalysisCheck: run.run_tdse == false
    TDSERunning --> AnalysisCheck: writes data/tdse/

    state AnalysisCheck <<choice>>
    AnalysisCheck --> AnalysisRunning: run.run_analysis == true
    AnalysisCheck --> Done: run.run_analysis == false
    AnalysisRunning --> Done: reads data/tise/, data/tdse/

    Done --> [*]
```

---

## 13. Original 4-Layer Architecture Sketch `[historical, 2026-06-20]` (`flowchart`)

*Source: `docs/planning/architecture-06-20.md`, "Architecture Diagram" (transcribed from handwritten notes `06-20-Notes-Part-1.jpg`). Reused verbatim.*

```mermaid
flowchart LR
    subgraph CFG["Input / Configuration"]
        direction TB
        POT["Define potential\nand domain"]

        subgraph STRUCT_IN["STRUCTURAL"]
            direction TB
            SI["How many bound states\nExtra boundary conditions\nWhich range of continuum\nHow to collocate points"]
        end

        subgraph PROP_IN["PROP."]
            direction TB
            PI["TDSE input\n→ Define operator\n→ Define TD field\n→ Time stepping"]
        end

        ANA["Analysis:\nwhat to compute"]
    end

    BS["B-spline nodes\nchosen accordingly"]

    subgraph PROG["Programs  (C++)"]
        direction TB
        SP["Structured part\nC++"]
        TD["TDSE\nC++"]
    end

    subgraph OUT["Output"]
        direction TB
        OUT_SP["Energies · Wavefunctions\nMatrix elements · Phase shifts"]
        OUT_TD["Wave packets\nTD observables"]
    end

    PY["Python\nAnalysis"]
    PLT["Spectra &\nfancy plots"]

    POT  -->|node placement| BS
    BS   --> SP
    POT  --> SP
    SI   --> SP
    PI   --> TD
    ANA  --> PY

    SP     --> OUT_SP
    OUT_SP -->|eigenbasis| TD
    TD     --> OUT_TD

    OUT_SP --> PY
    OUT_TD --> PY
    PY     --> PLT
```

This is the earliest architecture sketch (predates the Controller and the local-storage schema, both introduced in the 06-26 doc). Kept here as the historical starting point of the design.

---

## 14. Design Evolution Timeline (`timeline`)

*Source: dates and content of all planning docs cross-referenced above.*

```mermaid
timeline
    title Design Evolution -- 1D-QM-Playground
    2026-06-20 : Initial 4-layer sketch from handwritten notes
               : Input -> Structured part (TISE) -> TDSE -> Python Analysis
    2026-06-26 : Architecture Design Document formalized
               : Controller introduced as 5th component (Python orchestrator)
               : subprocess + YAML chosen for the Python-C++ interface
               : Local storage schema (data/tise/, data/tdse/) defined
               : TDSE originally specified as a Chebyshev / van Dijk propagator
    2026-06-28 : config.yaml schema locked
               : tise.n_states removed -- bound count is output, not input
               : tdse.chebyshev_order removed -- propagator is an implementation detail
               : analysis block split from visualization block
    2026-07-02 : TISE / TDSE / Analysis I-O notation reference written
    2026-07-03 : Stakeholder decisions incorporated into architecture-06-20.md
               : dynamics operators restricted to H, x, p_x -- stay strictly 1D
               : boundary conditions -- 3 domain configs, Case 1-2-3 asymptote logic
               : bound-state count confirmed as solver output, not input
               : node placement -- uniform default plus strategic knots, WKB deferred
               : FEDVR explored as a deferred future basis option
               : SDD.md opened as the authoritative document -- skeleton, sections pending
    Post-06-28 : architecture-06-26.md revised in place to the general
               : driven-TDSE design (gauge choice, optional field), superseding
               : its own original Chebyshev description, which is preserved
               : separately under docs/TDSE-original-design/
```

Two documents share the name `architecture-06-26.md` with different content: the one under `docs/planning/` was updated in place to match the current, driven-Hamiltonian TDSE design; the one under `docs/TDSE-original-design/` (and its duplicate at the top-level `Quantum Mechanics BSplines C++ Project/architecture-06-26.md`) preserves the original Chebyshev/van Dijk design as an archival snapshot.
