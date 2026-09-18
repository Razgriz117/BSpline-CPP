# ADR-0017: `physics.mass`/`physics.hbar` Generalization (Completes G12)

- **Status:** Accepted (implemented)
- **Date:** 2026-09-18

## Context

`physics.mass`/`physics.hbar` were documented config fields from early on (`docs/SDD.md`'s schema table), but every physics formula in the TISE solver — the kinetic-energy matrix element, `computeEAcc`, and `matchAsymptotic`'s continuum wavenumber/Sommerfeld parameter/energy-normalization amplitude — hardcoded atomic units (`mass=1`, `hbar=1`) internally, silently ignoring the config values.

`docs/planning/tise-release-readiness-plan.md`'s gap **G12** flagged this and Part C closed it as a **guard-rail, not full generalization**: `tise_solver_main.cpp` was made to reject any `physics.mass`/`physics.hbar` value other than exactly `1.0`, turning a previously silent wrong-physics bug into an honest config error. Full generalization was explicitly deferred.

## Decision

Complete G12's deferred half: `mass`/`hbar` become real, honored parameters throughout the TISE solver.

1. **Kinetic-energy matrix element** (`fillBandedMatrices`, `TISE/tise.cpp`): `hbar^2/(2*mass)` replaces the fixed `1/2` factor.
2. **`computeEAcc`**: gains an `hbar` parameter (it already took `mass`); formula becomes `hbar^2*pi^2/(2*mass*nodeSpacing^2)`.
3. **`matchAsymptotic`**: the continuum wavenumber becomes `k=sqrt(2*mass*E)/hbar`; the Sommerfeld parameter becomes `eta=mass*C/(hbar^2*k)` (derived from the radial Schrödinger equation with a `C/r` tail — the reduced Coulomb ODE itself needs no change, only correct `eta`/`rho` inputs); the energy-normalization amplitude `A_E` (both the flat-tail and Coulomb branches) picks up the general `sqrt(mass)/hbar`-flavored factor derived from the `dk/dE=mass/(hbar^2*k)` Jacobian that converts `delta(k-k')`- to `delta(E-E')`-normalization — this is the one part of the change that required an actual re-derivation, not a mechanical substitution, verified against `docs/SDD.md:555`'s documented (atomic-units-only) target formula, which the general form collapses to exactly at `mass=hbar=1`.
4. **`solveTISE`** and **`analyticHydrogenEnergy`/`eigenvalueError`** gain the same `mass`/`hbar` parameters, so the all-in-one entry point (used by unit tests and the legacy `TISE/main.cpp` driver) and the hydrogen analytic-comparison helper stay consistent with the rest of the solver.
5. **Guard-rail replaced**: `tise_solver_main.cpp` now defaults `mass`/`hbar` to `1.0` when absent (unchanged behavior for every config that never set `physics:`) and rejects only non-positive values, not non-unity ones.
6. All new parameters are **trailing, defaulted** (`=1.0`), matching this codebase's existing extension convention for `fillBandedMatrices`/`matchAsymptotic` — no existing call site needed updating.

`radialPotential` (`tise.cpp:55-58`, dead code "kept for reference," not on the production path) was deliberately left untouched — generalizing it would add risk for zero production benefit.

## Consequences

- A config can now set `physics.mass`/`physics.hbar` to any positive value and get physically correct results, not a config error.
- `TISE/main.cpp`'s `PARTICLE_MASS=3.75` constant (used for its Gaussian wavepacket shape) is now actually threaded into `solveTISE`/`eigenvalueError`, closing a pre-existing latent inconsistency where the eigenbasis being propagated silently assumed `mass=1` regardless.
- A future TDSE propagator built on top of TISE's `H`/`S` matrices and eigenbasis inherits correctly-generalized eigenstates for any configured mass/ħ, rather than eigenstates silently locked to atomic units — this was true whether or not this ADR landed before TDSE work began, but landing it first avoids TDSE ever having to reconcile a general-mass propagator against a mass-locked eigenbasis.
- `docs/guides/input-file-authoring-guide.md` §3.2/§7/§8, `docs/SDD.md`'s field table and glossary, and two near-duplicate schema-table copies (`docs/TDSE-original-design/`, `docs/superpowers/specs/`) are updated in the same pass to describe this real behavior instead of the old "always 1.0, anything else errors" framing.
- Covered by new unit tests (`TISE/tests/test_tise.cpp`: `ComputeEAccTest.ScalesWithHbarSquared`, `FillBandedMatricesTest.KineticTermScalesWithMassAndHbarOverlapUnchanged`, `FreeParticleContinuumGeneralUnitsTest` (3 cases), `AnalyticHydrogenEnergyTest.ScalesWithMassAndInverseHbarSquared`, `EigenvalueErrorTest.RespectsMassAndHbar`, `SolveTISETest.GroundStateMatchesAnalyticHydrogenAtNonUnitMassAndHbar`) and new real-subprocess integration tests (`tests/free_particle_general_units.yaml` — the first `tests/*.yaml` reference config to set a non-default `physics:` block — plus `TestFreeParticleGeneralUnitsRealSubprocess` and the replacement guard-rail tests in `test_controller_integration.py`).

## Source

`TISE/tise.cpp`/`tise.hpp` (`fillBandedMatrices`, `computeEAcc`, `matchAsymptotic`, `solveTISE`, `analyticHydrogenEnergy`/`eigenvalueError`); `TISE/tise_solver_main.cpp` (guard-rail); `TISE/main.cpp`; `docs/planning/tise-release-readiness-plan.md` (G12/Part C, the deferral this completes); `docs/planning/tise-mass-hbar-generalization-plan.md` (this work's own detailed plan + implementation record); `docs/SDD.md:555` (the atomic-units-only target formula the general `A_E` derivation was checked against).
