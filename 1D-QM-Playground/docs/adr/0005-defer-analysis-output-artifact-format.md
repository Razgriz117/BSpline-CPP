# ADR-0005: Defer Analysis Output-Artifact Format

- **Status:** Accepted (deferred)
- **Date:** 2026-07-28 (identified and decided during Phase 3 — Controller↔Analysis interface — implementation planning, SDD §7.2.3/§10.2)

## Context

SDD §7.2.3 (Controller to Analysis) describes Analysis's outputs as "plots/derived data (format not yet fixed — matures alongside ADR-0004's revisit trigger)." On close reading, that cross-reference is imprecise: ADR-0004's own Context/Decision/Consequences (and every other place ADR-0004 is cited — §1.2, §3.3, §5.4.2, §6.1) are scoped narrowly to the `visualization` config block's plot *parameters* — axis ranges, which state indices or energy windows to include, log vs. linear scale, etc. — for a plot that is assumed to already exist. ADR-0004 says nothing about, and was never intended to cover, the broader and logically prior question of whether Analysis produces an output artifact at all, and if so, in what form: a file (and if so, what format, and where — e.g. a `data/analysis/` directory analogous to `data/tise/`/`data/tdse/`, never documented in SDD §6.3 or §11.2), stdout data, or nothing.

This gap was surfaced during Phase 3 implementation planning (SDD §10.2: "Controller ↔ Analysis interface," SDD §7.2.3): Phase 3's job is to implement the subprocess/exit-code contract and `--tise-dir` consumption described in §7.2.3, which is fully specifiable and testable without resolving this question — REQ-F-060's actual computable quantities (bound-state populations, expectation values, spectral distributions) all require TDSE output per §6.1's `analysis` block ("all fields require `run_tdse: true`"), and TDSE doesn't exist until Phase 8 (§9.4 ties REQ-F-060's real test creation to Phase 8, not Phase 3). So Phase 3 has nothing yet to plot or persist, and no document anywhere (SDD.md, any ADR, `docs/planning/`, or the archived `docs/TDSE-original-design/`) actually specifies what Analysis should produce once it does.

## Decision

Defer specification of Analysis's output artifact — its format, its location (file vs. stdout; if file-based, its directory/naming, and whether a `--output-dir`-style CLI flag is ever added, unlike `tise_solver`/`tdse_solver`'s existing one per §7.1) — until REQ-F-060's quantities are actually computable, i.e. until TDSE output exists (SDD §10.2 Phase 8) and/or §7.2.3 is revisited and extended for TDSE outputs (SDD §10.2 Phase 7).

Phase 3 (this phase) implements only the subprocess/exit-code/`--tise-dir`-consumption contract of §7.2.3 and produces **no output artifact of any kind** — no plot files, no derived-data files, no `data/analysis/` directory, no placeholder stdout payload. `analysis.py`'s CLI exit code (`0` success / non-zero failure, per §8) is Phase 3's entire observable success signal, exactly as `tise_solver`'s was before Phase 4 filled in real numerics.

## Consequences

- SDD §7.2.3's "Outputs" bullet is corrected, in the same documentation pass that introduces this ADR, to say Phase 3 produces no output artifact and to cite ADR-0005 (not ADR-0004) for the still-open format/location question.
- ADR-0004 remains scoped exactly as it already is everywhere it's cited (narrowly, to `visualization` plot *parameters*) — this ADR does not change or reinterpret ADR-0004; it documents a previously-uncaptured, broader question ADR-0004 was never actually answering.
- A future phase must decide, before or as part of implementing it: whether Analysis writes files (and if so, under what directory/naming — e.g. `data/analysis/`, added to §6.3's persistent-storage table and §11.2's directory layout, neither of which mentions it today), prints structured data to stdout, or some combination; and whether `analysis.py` needs its own `--output-dir` flag added to §7.1's External Interfaces table (which today, unlike `tise_solver`/`tdse_solver`, deliberately has none for `analysis.py`).
- **Revisit trigger:** SDD §10.2 Phase 7 (Controller↔Analysis interface, revisited/extended for TDSE outputs) and/or Phase 8 (TDSE Solver implementation) — whichever first actually requires Analysis to produce a real, consumable output from computed REQ-F-060 quantities. Until then, `analysis.py`'s CLI is expected to keep producing no output artifact, by design, not by omission.

## Source

SDD §7.2.3 (Controller to Analysis) and §10.2 (Phase 3 / Phase 7 rows); identified during Phase 3 (Controller↔Analysis interface) implementation planning, 2026-07-28. Cross-reference: [ADR-0004](0004-defer-visualization-plot-parameters.md) for the distinct, narrower plot-*parameter* schema question this ADR does not reopen.
