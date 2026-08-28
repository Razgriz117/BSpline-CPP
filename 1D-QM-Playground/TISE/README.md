# H-BoundStates: Hydrogenic Bound States with B-Splines

## Table of Contents
- [Code Structure](#code-structure)
- [Specifying the Potential](#specifying-the-potential)
- [Dependencies](#dependencies)
- [Building with CMake](#building-with-cmake)
- [Building and Running Tests](#building-and-running-tests)
- [Building with `g++` (manual build)](#building-with-g-manual-build)
- [Coverage](#coverage)
- [Output Files](#output-files)
- [Known Limitations](#known-limitations)

This project computes bound-state eigenvalues and eigenfunctions for a 1D/radial quantum system using a **B-spline basis** and a **generalized eigenvalue problem**. The potential is supplied at runtime as a piecewise expression (see [Specifying the Potential](#specifying-the-potential)) rather than hardcoded, so the same binary solves any potential expressible as a sum of closed-form pieces — not just the hydrogen atom.

The code constructs a B-spline basis on a radial grid, builds the Hamiltonian and overlap matrices from the supplied potential, and solves

\[
H \mathbf{c} = E\, S \mathbf{c}
\]

using LAPACK’s `DSBGV` (symmetric banded generalized eigenproblem). If the potential you supply matches the hydrogenic form \(V_\ell(r) = \ell(\ell+1)/2r^2 - 1/r\) for the angular momentum \(\ell\) set at compile time (see [Specifying the Potential](#specifying-the-potential)), the resulting eigenvalues are also compared to the analytic hydrogenic energies

\[
E_n^{(\ell)} = -\frac{1}{2(n_{\text{eff}})^2}, \quad n_{\text{eff}} = n + \ell,
\]

printed alongside each eigenvalue as an error column. This comparison is only meaningful when the supplied potential actually is hydrogenic — for any other potential the printed "error" isn't validated against what you passed in and can be ignored. The corresponding eigenfunctions are written to files `EigenState_XXX` in the build directory.

## Code Structure

- `BSpline.hpp` / `BSpline.cpp`  
  C++ implementation of a B-spline basis (closely modeled after Luca Argenti’s Fortran `ModuleBSpline`).  
  Handles:
  - B-spline grid construction with knot multiplicities at the endpoints
  - Polynomial coefficients for each B-spline segment
  - Evaluation of B-splines and their derivatives at arbitrary points
  - Numerical integration via Gauss–Legendre quadrature

- `main.cpp`  
  Drives the solve, given a potential passed on the command line:
  - Parses the potential (see [Specifying the Potential](#specifying-the-potential)) from `argv[1]`
  - Builds a uniform radial grid \([r_{\min}, r_{\max}]\)
  - Constructs Hamiltonian and overlap (mass) matrices in banded storage from the parsed potential
  - Calls LAPACK’s `DSBGV` to solve for eigenvalues/eigenvectors
  - Prints accurate eigenvalues and (hydrogenic-comparison) errors
  - Writes eigenfunctions to `EigenState_XXX` files sampled on a radial grid
  - Propagates a Gaussian wavepacket in time after the bound-state solve (see [Known Limitations](#known-limitations))

## Specifying the Potential

The potential is not hardcoded — it is passed as `argv[1]`, a JSON array of `{"domain": "<interval>", "function": "<expression in x>"}` pieces, one per region of a piecewise-defined potential. `function` is parsed and evaluated by [muparser](https://beltoforion.de/en/muparser/) at each quadrature point; `domain` is a standard interval notation (`[`/`(` and `]`/`)` for closed/open ends, `inf`/`-inf` allowed).

For example, `main.cpp`'s own doc comment gives:

```bash
./H-BoundStates '[{"domain": "[0, 20)", "function": "x"}, {"domain": "[20, 40]", "function": "x^2"}]'
```

which solves a potential that is linear on \([0, 20)\) and quadratic on \([20, 40]\). To reproduce the hydrogenic potential \(V_\ell(r) = \ell(\ell+1)/2r^2 - 1/r\) for a given \(\ell\), supply that full expression yourself — e.g. for \(\ell=1\):

```bash
./H-BoundStates '[{"domain": "(0, 100]", "function": "-1/x + 1/x^2"}]'
```

`main.cpp` also has a compile-time constant:

```cpp
// Select orbital angular momentum l (integer >= 0).
constexpr int L = 0;  // set to 0, 1, 2, ...
```

Unlike in earlier versions of this project, `L` no longer shapes the potential itself — the potential is entirely whatever `argv[1]` describes. `L` now only feeds `eigenvalueError`'s analytic-hydrogen comparison column (\(E_n^{(\ell)} = -1/(2(n+\ell)^2)\)); it is not validated against the potential you actually passed in. If you set `L` to a value inconsistent with your supplied potential (e.g. leaving `L=0` while passing the \(\ell=1\) expression above), the printed "error" column no longer means anything. Set `L` to match whatever hydrogenic \(\ell\) your supplied potential corresponds to (or ignore the error column entirely for non-hydrogenic potentials) and rebuild.

---

## Dependencies

* A C++ compiler with C++17 support (e.g. `g++ >= 7`)
* [nlohmann's C++ JSON library](https://github.com/nlohmann/json#quick-reference)
* [muparser](https://beltoforion.de/en/muparser/) — expression parser used to evaluate each potential piece's `function` string
* LAPACK and BLAS libraries (for `dsbgv_`)
* Eigen C++ library.

On many Linux systems these can be installed with your package manager, e.g.:

```bash
sudo apt-get install liblapack-dev libblas-dev libeigen3-dev nlohmann-json3-dev libmuparser-dev
```

---

## Building with CMake

From the project root:

```bash
cmake -S . -B build
cmake --build build
```

This will create the `H-BoundStates` executable inside the `build` directory.

### Running the executable

Change into the build directory and run from there, passing the potential as a JSON argument (see [Specifying the Potential](#specifying-the-potential)). The example below is the \(\ell=0\) hydrogenic potential, matching `main.cpp`'s default `constexpr int L = 0`, so the printed error column is meaningful out of the box:

```bash
cd build
./H-BoundStates '[{"domain": "(0, 100]", "function": "-1/x"}]'
```

The program will:

* Print the number of accurate eigenvalues and their energies/errors.
* Write eigenfunctions to files `EigenState_XXX` in the **build directory**.

---

## Building and Running Tests

The test suite uses [GoogleTest](https://github.com/google/googletest) and is built separately from the main executable via the `BUILD_TESTING` CMake option. GTest must be installed on your system first:

```bash
sudo apt-get install libgtest-dev
```

To configure and build with tests enabled:

```bash
cmake -S . -B build -DBUILD_TESTING=ON
cmake --build build
```

To run all tests:

```bash
cd build
ctest --output-on-failure
```

This runs two test suites:

| Suite | Binary | What it covers |
|---|---|---|
| `BSplineTests` | `test_bspline` | `BSpline` class — init error codes, `free()`, getters, partition of unity, derivative of PoU, linear combination eval, `skipFirst`, `Bsmin`/`Bsmax`, integrals (symmetry, bounds, derivative orders, parvec), multiple orders, non-uniform grids |
| `UtilsTests` | `test_utils` | `buildAdjointMatrix` correctness, `printBandMatrixAsDense` smoke tests |

The test sources live in `tests/`.

---

## Building with `g++` (manual build)

From the project root directory (`TISE/`). This mirrors what `CMakeLists.txt` does for the `H-BoundStates` target:

1. Compile each translation unit `H-BoundStates` actually links (note `-std=c++17`, not `-std=cpp17`):

   ```bash
   g++ -O2 -std=c++17 -c BSpline.cpp
   g++ -O2 -std=c++17 -c tise.cpp $(pkg-config --cflags muparser)
   g++ -O2 -std=c++17 -I/usr/include/eigen3 -c time_evolution.cpp
   g++ -O2 -std=c++17 -I/usr/include/eigen3 -c main.cpp
   ```

2. Link them with LAPACK, BLAS, and muparser:

   ```bash
   g++ -O2 -std=c++17 -o H-BoundStates main.o tise.o time_evolution.o BSpline.o \
       -llapack -lblas $(pkg-config --libs muparser)
   ```

3. Create the `timesteps` output directory. `main.cpp` runs time evolution after the bound-state solve and writes per-step output under a directory named by the `TIMESTEPS_DIR` macro; the CMake build defines this to an absolute path and creates it via a custom target automatically, but this manual recipe doesn't define the macro, so `main.cpp` falls back to the relative path `./timesteps`, which must exist beforehand:

   ```bash
   mkdir -p timesteps
   ```

4. Run the executable **from the project directory (or wherever `H-BoundStates` resides)**, passing the potential as a JSON argument (see [Specifying the Potential](#specifying-the-potential)):

   ```bash
   ./H-BoundStates '[{"domain": "(0, 100]", "function": "-1/x"}]'
   ```

This will print eigenvalues and their errors to the terminal and create files like `EigenState_001`, `EigenState_002`, … in the same directory, plus per-timestep output under `timesteps/`.

---

## Coverage

Coverage instrumentation is opt-in via the `ENABLE_COVERAGE` CMake option (adds `--coverage` compile/link flags for `gcov`/`gcovr`):

```bash
cmake -S . -B build-coverage -DBUILD_TESTING=ON -DENABLE_COVERAGE=ON
cmake --build build-coverage -j
find build-coverage -name '*.gcda' -delete
ctest --test-dir build-coverage --output-on-failure
mkdir -p build-coverage/coverage
gcovr --root . --filter '.*tise\.cpp' --exclude '.*/tests/.*' --exclude '.*main\.cpp' \
      --print-summary --html --html-details -o build-coverage/coverage/index.html build-coverage
```

This project's coverage gate (≥80%) is tracked against `tise.cpp`'s own reported line coverage specifically, not a blended total across every source file — `gcovr --print-summary`'s combined number includes the always-100%-covered test file itself, which isn't the signal that matters. `build-coverage/` is excluded from version control by the existing `build*` pattern in `.gitignore`.

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

---

## Known Limitations

* **Time evolution always runs.** After solving the bound-state problem, `main.cpp` unconditionally calls `tevol::runTimeEvolution(...)` — there is currently no way to request a TISE-only run. `config.yaml` already defines a `run.run_tdse` flag for exactly this purpose (see `docs/SDD.md` §6.1), but `main.cpp` does not parse `config.yaml` (it has no YAML dependency at all) and so never reads it. Wiring this up is deferred to the Controller↔TISE configuration plumbing (`docs/planning/tise-task-breakdown.md` §4 item 4).

* **Continuum phase-shift matching assumes a regular right boundary.** If the supplied potential is singular at the domain's right edge (`x = rMax`), `solveTISE` still attempts continuum construction and prints a warning to stderr, but `matchAsymptotic`'s flat-asymptote matching formula is not valid there — treat any `phase_shifts.dat`/`continuum_state_NNN.dat` output from such a run with suspicion. See `docs/planning/engineer-a-plan-A4-wiring-design.md`.

