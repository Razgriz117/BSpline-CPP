# ADR-0016: Analysis Plot-Artifact Format (Supersedes ADR-0005)

- **Status:** Accepted (implemented)
- **Date:** 2026-09-11 (formalized during the TISE-Generalization final pre-merge
  pass, after `analysis.py`'s own module docstring was found contradicting its
  code — see Context)

## Context

ADR-0005 deferred Analysis's output-artifact format entirely, reasoning that
REQ-F-060's computable quantities all require TDSE output (unavailable until
Phase 8), so Phase 3 had "nothing yet to plot or persist." It decided Phase 3
would produce "no output artifact of any kind — no plot files, no
derived-data files, no `data/analysis/` directory, no placeholder stdout
payload."

That is no longer true, and has not been true since commits `8236239`/
`5cd28e4` added `plot_tise`, `plot_eigenstates`, and `plot_phase_shifts` to
`analysis.py` — none of which require TDSE output at all, since they plot
*TISE* output (continuum states, bound eigenstates, phase shifts), not the
REQ-F-060 quantities ADR-0005 was actually reasoning about. Despite this,
`analysis.py`'s module docstring and `run()`'s own docstring still asserted
"No plotting" / "Computes and writes nothing (ADR-0005)" directly above the
three `plt.savefig`/`fig.savefig` calls that do exactly that — a
self-contradiction discovered during a branch review, not a subtle drift.
`docs/SDD.md` (§1.2, its ADR-0005 traceability-table row, and §7.2.3) made
the identical stale claim.

## Decision

Supersede ADR-0005. Document the plotting behavior that already exists and
is already tested as the real, accepted output-artifact format for the TISE
half of Analysis's output (the REQ-F-060/TDSE-dependent half ADR-0005 was
actually about remains open — see Consequences):

1. **Location:** plots are written directly into the `--tise-dir` passed to
   `analysis.py` (i.e. `<output_dir>/tise/`) — not a separate
   `data/analysis/` directory. `analysis.py` gained no `--output-dir` flag of
   its own; it reuses the directory the Controller already told it about.
2. **`continuum_{idx}.png`** (`plot_tise`) — one per continuum energy,
   written unconditionally whenever continuum states are present (tolerates
   absence the same way the rest of Analysis tolerates `tise.continuum.
   enabled: false`, per §5.4.4). Not gated by any `visualization.*` flag.
3. **`eigenstate_{idx}.png`** (`plot_eigenstates`) — gated on
   `visualization.eigenstates`. Raw $\phi_n(x)$ by default; a state
   confirmed bound (see `docs/SDD.md` §6.1's `visualization.
   bound_states_squared` field for the exact boundness rule) is plotted as
   $|\phi_n(x)|^2$ instead when `visualization.bound_states_squared` is also
   set.
4. **`phase_shifts.png`** (`plot_phase_shifts`) — gated on
   `visualization.phase_shifts`; two stacked subplots, $\delta(E)$ and
   $d\delta/dE$ vs. $E$. Not gated on any TDSE-related flag.
5. Analysis's process exit-code contract (`0`/non-zero, §8) is unchanged;
   plotting failures are not distinguished from other Analysis failures.

## Consequences

- `analysis.py`'s module/`run()` docstrings and `docs/SDD.md` (§1.2's open-
  questions list, its ADR-0005 traceability row, §7.2.3's Controller-to-
  Analysis outputs bullet) are corrected, in the same pass that introduces
  this ADR, to describe this real behavior instead of asserting no artifact
  exists.
- ADR-0005's actual, narrower unresolved question — what Analysis should
  produce for REQ-F-060's TDSE-dependent quantities (bound-state
  populations, expectation values, spectral distributions) — remains
  genuinely open. This ADR does not resolve it; it only retires the parts of
  ADR-0005 that this branch's TISE-side plotting already answered by
  implementation. A future ADR should address the TDSE-dependent half if/
  when Phase 8 makes it concrete, rather than reopening this one.
- `visualization`'s other fields (`time_evolution`, `bound_state_
  populations`, `asymptotic_populations`, `asymptotic_distribution`,
  `expectation_values`, `interval_probability`) remain schema-only, awaiting
  TDSE — unaffected by this ADR.
- No `data/analysis/` directory was ever added, and none is planned by this
  decision; if a future TDSE-dependent artifact needs one, that's a decision
  for whichever future ADR addresses the remaining open question above.

## Source

`analysis.py` (`plot_tise`, `plot_eigenstates`, `plot_phase_shifts`, `run()` —
implemented in `8236239`/`5cd28e4`); `docs/SDD.md` §6.1 (`visualization.*`
field definitions, already accurate); [ADR-0005](0005-defer-analysis-output-artifact-format.md)
(the deferral this supersedes); `docs/planning/TISE-Generalization-Review.md`
(the review that surfaced the docstring/ADR contradiction).
