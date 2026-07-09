# ADR-0004: Defer Visualization Plot-Parameter Schema

- **Status:** Accepted (deferred)
- **Date:** 2026-06-28 (config schema design); reaffirmed 2026-07-09 (SDD population pass)

## Context

`config.yaml`'s `visualization` block controls which plots the pipeline produces, spanning both raw solver output (TISE eigenstates, phase shifts; TDSE time evolution) and `analysis`-computed quantities (bound-state populations, asymptotic distribution, expectation values). As specified today, every field in `visualization` is a boolean toggle — it says *whether* to produce a given plot, not *how* (axis ranges, which state indices or energy windows to include, log vs. linear scale, etc.).

Those finer plot parameters were raised during config schema design but not specified, because they depend on the concrete plotting logic of the Analysis module (SDD [§5.4](../SDD.md#54-analysis-module)), which had not yet been designed in detail at that time. Fixing a plot-parameter schema before that logic exists risks guessing at a shape the eventual implementation won't need, or missing one it will.

## Decision

Defer a detailed plot-parameter schema. `visualization` remains boolean-toggle-only (SDD [§6.1](../SDD.md#61-configuration-schema)) until the Analysis module's detailed design (SDD [§5.4](../SDD.md#54-analysis-module)) is complete enough to know what parameters its plotting logic actually needs.

## Consequences

- Analysis/plotting code in the near term either hard-codes reasonable plot defaults (ranges, subsets) or accepts them as ad hoc script arguments outside `config.yaml`, not as part of the versioned schema.
- **Revisit trigger: completion of SDD [§5.4](../SDD.md#54-analysis-module) (Analysis Module detailed design).** This is an explicit instruction from this decision's stakeholder, not just a general "later" — once [§5.4](../SDD.md#54-analysis-module) is fleshed out, this ADR should be revisited immediately as part of that same work, and `visualization` extended with whatever parameters the finished Analysis design requires.

## Source

`docs/superpowers/specs/2026-06-28-config-yaml-schema-design.md`, `visualization` section and Design Decisions table ("`visualization` scope"); referenced from SDD [§6.1](../SDD.md#61-configuration-schema) and [§5.4.2](../SDD.md#542-interfaces).
