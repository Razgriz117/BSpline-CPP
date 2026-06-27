# 1D-QM-Playground

A full B-spline based quantum mechanics solver for hydrogen-like potentials in one dimension. This directory serves as the starting point for ongoing development toward a complete 1D QM toolkit covering:

- **Bound and continuum eigenstates** of the TISE (Time-Independent Schrödinger Equation)
- **Time-domain dynamics** via the TDSE (Time-Dependent Schrödinger Equation) under external pulses

The long-term goal is a publication-quality solver suitable for submission to the *American Journal of Physics*.

---

## Directory Structure

```
1D-QM-Playground/
├── TISE/                  # C++ solver for the time-independent Schrödinger equation
│   ├── BSpline.cpp/.hpp   #   B-spline basis implementation
│   ├── main.cpp           #   Hydrogenic eigenvalue problem + time evolution prototype
│   ├── utils/             #   Utility functions
│   ├── CMakeLists.txt
│   └── make_and_run.sh
├── moduleBspline.f90      # Original Fortran B-spline module (Argenti)
├── Template.f90           # Fortran template for the eigenvalue problem
├── plot.py                # Python script for plotting eigenstates
├── heatmap.py             # Python script for heatmap visualization
└── output/                # Precomputed eigenstate data
    ├── L0/                #   ℓ = 0 bound states
    ├── L1/                #   ℓ = 1 bound states
    └── original/          #   Reference output from initial runs
```

For build instructions, dependencies, and code details, see [TISE/README.md](TISE/README.md).

---

## Scope and Roadmap

The `TISE/` C++ solver currently computes hydrogenic bound-state eigenvalues and eigenfunctions and includes a prototype time-evolution of a Gaussian initial state. Planned extensions include:

- Continuum eigenstates (scattering states above the ionization threshold)
- Full TDSE propagation under external pulses (e.g., laser fields)
- Python-based analysis pipeline for spectra and Fourier plots
