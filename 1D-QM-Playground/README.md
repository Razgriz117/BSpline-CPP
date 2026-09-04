# 1D-QM-Playground

A B-spline based quantum mechanics solver for arbitrary 1D potentials, driven end-to-end by a single `config.yaml`. Covers:

- **Bound and continuum eigenstates** of the TISE (Time-Independent Schrödinger Equation) — implemented and tested against known closed-form solutions (hydrogen, free particle, harmonic oscillator, finite square well).
- **Time-domain dynamics** via the TDSE (Time-Dependent Schrödinger Equation) under external pulses — not yet implemented; see `docs/planning/tdse-task-breakdown.md`.

The long-term goal is a publication-quality solver suitable for submission to the *American Journal of Physics*.

---

## Quick start (the real pipeline)

The actual, current pipeline is `config.yaml` → `controller.py` → the compiled `tise_solver` binary → `analysis.py`. Build the C++ solver once, then drive everything through `controller.py`:

```bash
# 1. Build the C++ TISE solver (produces TISE/build/tise_solver, requires
#    yaml-cpp in addition to H-BoundStates' own dependencies -- see
#    TISE/README.md for the full dependency list and manual g++ fallback).
cmake -S TISE -B TISE/build -DBUILD_TESTING=ON
cmake --build TISE/build

# 2. Install the Python dependencies (pyyaml, matplotlib, pytest).
pip install -r requirements.txt

# 3. Edit config.yaml (bspline basis, potential, continuum settings,
#    visualization toggles -- see the file's own inline comments for the
#    full schema), then run the whole pipeline:
python3 controller.py --config config.yaml
```

`controller.py` runs `tise_solver` (writing `eigenvalues.dat`, `eigenvectors.dat`, `eigenstate_NNN.dat`, and — if `tise.continuum.enabled` — `phase_shifts.dat`/`continuum_state_NNN.dat`, plus a `warnings.json` sidecar for anything physics-relevant it noticed) into `<run.output_dir>/tise`, then — if `run.run_analysis: true` — runs `analysis.py` against that output to produce plots (`continuum_*.png`, and, per the `visualization.eigenstates`/`visualization.phase_shifts` toggles, `eigenstate_*.png`/`phase_shifts.png`).

Each stage is also independently runnable — see `docs/SDD.md` §7 for the full inter-component contract, and `tests/*.yaml` (below) for ready-to-run example configs.

### Known-solution reference configs

`tests/{free_particle,finite_square_well,harmonic_oscillator,hydrogen,interior_singularity,right_edge_singularity,case3_irregular_tail}.yaml` are runnable configs for potentials with a closed-form analytic answer, meant for visually/numerically spot-checking solver output against known physics (not just passing automated tests):

```bash
python3 controller.py --config tests/free_particle.yaml
```

produces `./data/free_particle/tise/continuum_*.png` etc. — see each file's own settings for what it's checking (e.g. `finite_square_well.yaml`'s genuinely energy-dependent phase shift vs. `free_particle.yaml`'s flat `delta=0`). These same configs are also exercised by the automated suite in `tests/test_analysis_integration.py`.

---

## Directory Structure

```
1D-QM-Playground/
├── config.yaml            # Single source of truth for a real pipeline run
├── controller.py           # Orchestrates tise_solver (+ analysis.py) per config.yaml
├── analysis.py              # Reads data/tise/ output, produces plots
├── requirements.txt         # Python dependencies (pyyaml, matplotlib, pytest)
├── TISE/                  # C++ solver
│   ├── tise_solver_main.cpp #   CLI for the real, config-driven `tise_solver` binary
│   ├── tise.cpp/.hpp        #   Core TISE numerics (bound + continuum states)
│   ├── BSpline.cpp/.hpp    #   B-spline basis implementation
│   ├── main.cpp             #   Standalone `H-BoundStates` hydrogen demo (see TISE/README.md)
│   ├── utils/                #   Utility functions
│   ├── tests/                #   GoogleTest suite (BSpline, utils, TISE)
│   └── CMakeLists.txt
├── tests/                  # Python test suite (pytest) + known-solution reference configs
├── docs/                   # SDD, ADRs, planning docs
├── moduleBspline.f90         # Original Fortran B-spline module (Argenti) -- historical
├── Template.f90              # Fortran template for the eigenvalue problem -- historical
├── plot.py / heatmap.py       # Older standalone plotting scripts (pre-analysis.py)
└── output/                  # Precomputed eigenstate data from the original Fortran runs -- historical
```

For C++ build instructions and dependencies (including the standalone `H-BoundStates` demo), see [TISE/README.md](TISE/README.md).

---

## Scope and Roadmap

Implemented: general 1D piecewise potentials (not just hydrogenic), bound-state eigenvalues/eigenvectors, continuum eigenstates (phase shifts + scattering wavefunctions) with strategic B-spline node placement for potentials with internal structure, and a Python analysis/plotting layer — all driven by `config.yaml` through `controller.py`.

Not yet implemented: TDSE propagation under external time-dependent fields (`run.run_tdse` is explicitly rejected by `controller.py` today) and the Coulomb-tail continuum-matching formula for potentials with a genuine `1/r`-type tail beyond the box (the implemented continuum matching assumes a flat/free-particle-like asymptote; see `docs/SDD.md` §5.2.3 and `docs/adr/`).
