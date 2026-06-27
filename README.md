# BSpline-CPP

A C++ B-spline toolkit for solving quantum mechanics problems on a radial grid. This repository serves two purposes:

1. **Student base code** for the PHY5606 Hydrogen Atom project
2. **A full working solver** that is being extended into a publication-quality 1D QM toolkit

---

## Repository Structure

```
BSpline-CPP/
├── student-base-code/     # Base code distributed to PHY5606 students
└── 1D-QM-Playground/     # Full solver — starting point for ongoing development
```

### `student-base-code/`

A trimmed version of the B-spline codebase provided to students as the foundation for the [PHY5606 Hydrogen Atom project](../PHY5606_F25_Projects.pdf). It includes the B-spline basis implementation, a skeleton `main.cpp`, and a unit test suite. Students build on this code to compute hydrogenic bound-state eigenvalues and eigenfunctions.

See [student-base-code/README.md](student-base-code/README.md) for build instructions and code details.

### `1D-QM-Playground/`

The complete solver, expanding from bound states to include continuum eigenstates (TISE) and time-domain dynamics under external pulses (TDSE). This directory also contains the reference Fortran B-spline module and Python analysis scripts.

See [1D-QM-Playground/README.md](1D-QM-Playground/README.md) for an overview of the structure and roadmap, and [1D-QM-Playground/TISE/README.md](1D-QM-Playground/TISE/README.md) for C++ build and usage details.

---

## Background

The B-spline method represents the radial wavefunction in a basis of piecewise polynomials defined on a uniform radial grid. The TISE reduces to a generalized eigenvalue problem

$$H \mathbf{c} = E\, S \mathbf{c}$$

solved with LAPACK's `DSBGV` routine. The basis and method follow the approach described in:

> H. Bachau *et al.*, *Rep. Prog. Phys.* **64**, 1815 (2001).
