# H-BoundStates: Hydrogenic Bound States with B-Splines

## Table of Contents
- [Code Structure](#code-structure)
- [Changing the Angular Momentum](#changing-the-angular-momentum)
- [Dependencies](#dependencies)
- [Building with CMake](#building-with-cmake)
- [Building with `g++` (manual build)](#building-with-g-manual-build)
- [Output Files](#output-files)

This project computes bound-state eigenvalues and eigenfunctions for a hydrogen-like atom using a **B-spline basis** and a **generalized eigenvalue problem**.

The radial Schrödinger equation is solved for the potential

\[
V_\ell(r) = \frac{\ell(\ell+1)}{2r^2} - \frac{1}{r},
\]

where \(\ell\) is the orbital angular-momentum quantum number.  
The code constructs a B-spline basis on a radial grid, builds the Hamiltonian and overlap matrices, and solves

\[
H \mathbf{c} = E\, S \mathbf{c}
\]

using LAPACK’s `DSBGV` (symmetric banded generalized eigenproblem). The resulting eigenvalues are compared to the analytic hydrogenic energies

\[
E_n^{(\ell)} = -\frac{1}{2(n_{\text{eff}})^2}, \quad n_{\text{eff}} = n + \ell,
\]

and the corresponding eigenfunctions are written to files `EigenState_XXX` in the build directory.

## Code Structure

- `BSpline.hpp` / `BSpline.cpp`  
  C++ implementation of a B-spline basis (closely modeled after Luca Argenti’s Fortran `ModuleBSpline`).  
  Handles:
  - B-spline grid construction with knot multiplicities at the endpoints
  - Polynomial coefficients for each B-spline segment
  - Evaluation of B-splines and their derivatives at arbitrary points
  - Numerical integration via Gauss–Legendre quadrature

- `main.cpp`  
  Sets up the hydrogenic radial problem:
  - Defines the angular-momentum parameter `L` (an integer \(\ell \ge 0\))
  - Builds a uniform radial grid \([r_{\min}, r_{\max}]\)
  - Constructs Hamiltonian and overlap (mass) matrices in banded storage
  - Calls LAPACK’s `DSBGV` to solve for eigenvalues/eigenvectors
  - Prints accurate eigenvalues and errors
  - Writes eigenfunctions to `EigenState_XXX` files sampled on a radial grid

## Changing the Angular Momentum

In `main.cpp`, there is a compile-time constant:

```cpp
// Select orbital angular momentum l (integer >= 0).
constexpr int L = 0;  // set to 0, 1, 2, ...
````

Set `L` to the desired (\ell) value and rebuild. Both the potential and the analytic comparison used for the eigenvalue error will update accordingly.

---

## Dependencies

* A C++ compiler with C++17 support (e.g. `g++ >= 7`)
* LAPACK and BLAS libraries (for `dsbgv_`)
* Eigen C++ library.

On many Linux systems these can be installed with your package manager, e.g.:

```bash
sudo apt-get install liblapack-dev libblas-dev
```

---

## Building with CMake

If you have CMake installed, you may run the following commands from the project root:

```bash
cmake -S . -B build
cmake --build build
```

This will create the `H-BoundStates` executable inside the `build` directory.

### Running the executable

Change into the build directory and run from there:

```bash
cd build
./H-BoundStates
```

The program will:

* Print the number of accurate eigenvalues and their energies/errors.
* Write eigenfunctions to files `EigenState_XXX` in the **build directory**.

### Make and Run Script

On Linux, you may also combine the above steps by invoking the `make_and_run.sh` script as follows:
```
# You may need to add execution permission with chmond +x prior to first usage,
# but this will not be necessary on subsequent executions of the script.
# Example:
#   chmod +x make_and_run.sh

./make_and_run.sh
```

---

## Building with `g++` (manual build)

From the project root directory:

1. Compile the source files:

   ```bash
   g++ -O2 -std=cpp17 -c BSpline.cpp
   g++ -O2 -std=cpp17 -c main.cpp
   ```

2. Link them with LAPACK and BLAS:

   ```bash
   g++ -O2 -std=cpp17 -o H-BoundStates main.o BSpline.o -llapack -lblas
   ```

3. Run the executable **from the project directory (or wherever `H-BoundStates` resides)**:

   ```bash
   ./H-BoundStates
   ```

This will print eigenvalues and their errors to the terminal and create files like `EigenState_001`, `EigenState_002`, … in the same directory.

---

## Output Files

* **Terminal output**
  For each accurate eigenvalue:

  ```text
   n   E_n (numeric)           E_n - E_n^analytic (or similar error)
  Number of Accurate Eigenvalues : N
  ```

* **Eigenstate files** (`EigenState_001`, `EigenState_002`, …)
  Each file contains two columns:

  ```text
  r        ψ_n(r)
  ```

  sampled on a uniform grid between `BS_GRMIN` and `BS_GRMAX` (as set in `main.cpp`).

You can plot these eigenfunctions using your favorite tool (Python/Matplotlib, gnuplot, etc.) to visualize the radial bound states.

