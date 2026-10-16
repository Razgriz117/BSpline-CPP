# TISE: B-Spline Time-Independent Schrödinger Equation Solver

> **This file documents two separate binaries.** `tise_solver` (below) is the real, config-driven binary the project's actual pipeline (`config.yaml` → `controller.py` → `tise_solver` → `analysis.py`, see the top-level [README.md](../README.md)) uses — build and run that first if you just want to solve a potential. `H-BoundStates` (everything after it in this file) is an earlier, standalone hydrogen demo taking its potential as a raw JSON `argv` string instead of `config.yaml` — still maintained and tested, useful for quick one-off checks, but not part of the config-driven pipeline.

## Table of Contents
- [Building and Running `tise_solver`](#building-and-running-tise_solver)
- [H-BoundStates: Hydrogenic Bound States with B-Splines](#h-boundstates-hydrogenic-bound-states-with-b-splines)
- [Code Structure](#code-structure)
- [Specifying the Potential](#specifying-the-potential)
- [Dependencies](#dependencies)
- [Why there is no LAPACK here](#why-there-is-no-lapack-here)
- [Building with CMake](#building-with-cmake)
- [Building and Running Tests](#building-and-running-tests)
- [Building with `g++` (manual build)](#building-with-g-manual-build)
- [Coverage](#coverage)
- [Output Files](#output-files)
- [Known Limitations](#known-limitations)

## Building and Running `tise_solver`

```bash
# From the project root (one level above this TISE/ directory):
cmake -S TISE -B TISE/build -DBUILD_TESTING=ON -DCMAKE_BUILD_TYPE=Release
cmake --build TISE/build --config Release
```

This builds both `H-BoundStates` and `tise_solver`. You do not need to install
any C++ libraries first: anything missing is downloaded and built by CMake (see
[Dependencies](#dependencies)). If you want `H-BoundStates` alone, pass
`-DBUILD_TISE_SOLVER=OFF`.

Run it directly against any `config.yaml` (see the top-level README's "Quick start" and `config.yaml`'s own inline comments for the schema):

```bash
./TISE/build/tise_solver --config config.yaml --output-dir data/tise
```

This writes `eigenvalues.dat`, `eigenvectors.dat`, `eigenstates.dat`, `eigenstate_NNN.dat` (one per bound state, unconditionally), `hamiltonian.dat`, `overlap.dat`, and `warnings.json` always; `phase_shifts.dat`, `continuum_states.dat` and `continuum_state_NNN.dat` additionally when `tise.continuum.enabled: true`.

**`eigenstates.dat` and `continuum_states.dat` are the ones to reach for if you are writing your own analysis.** Each is a single table — column 1 is `x`, column `n+1` is the n-th wavefunction — with every eigenvalue (and, for the continuum, every phase shift) repeated in the `#` header, so the file needs nothing else opened alongside it:

```python
d = np.loadtxt('eigenstates.dat')     # (n_points, n_states + 1)
plt.plot(d[:, 0], d[:, 3])            # psi_3
```

```gnuplot
plot 'eigenstates.dat' using 1:4 with lines   # psi_3
```

The per-state `eigenstate_NNN.dat` / `continuum_state_NNN.dat` files are still written and unchanged; the tables are additive. Note that those filenames are 1-based while `eigenvalues.dat`'s index column is 0-based — an off-by-one the consolidated tables avoid entirely. Non-zero exit means no partial output is left behind (see `docs/SDD.md` §7.2.1). Ordinarily you'd drive this via `controller.py` rather than invoking it directly — see the top-level README.

---

# H-BoundStates: Hydrogenic Bound States with B-Splines

This project computes bound-state eigenvalues and eigenfunctions for a 1D/radial quantum system using a **B-spline basis** and a **generalized eigenvalue problem**. The potential is supplied at runtime as a piecewise expression (see [Specifying the Potential](#specifying-the-potential)) rather than hardcoded, so the same binary solves any potential expressible as a sum of closed-form pieces — not just the hydrogen atom.

The code constructs a B-spline basis on a radial grid, builds the Hamiltonian and overlap matrices from the supplied potential, and solves

\[
H \mathbf{c} = E\, S \mathbf{c}
\]

using Eigen's `GeneralizedSelfAdjointEigenSolver` (see [Why there is no LAPACK here](#why-there-is-no-lapack-here)). If the potential you supply matches the hydrogenic form \(V_\ell(r) = \ell(\ell+1)/2r^2 - 1/r\) for the angular momentum \(\ell\) set at compile time (see [Specifying the Potential](#specifying-the-potential)), the resulting eigenvalues are also compared to the analytic hydrogenic energies

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
  - Calls `tise::solveTISE(...)` (see below) to do the actual work
  - Prints accurate eigenvalues and (hydrogenic-comparison) errors
  - Writes eigenfunctions to `EigenState_XXX` files sampled on a radial grid
  - Propagates a Gaussian wavepacket in time after the bound-state solve (see [Known Limitations](#known-limitations))

- `tise.hpp` / `tise.cpp`  
  The core solver library. Both `H-BoundStates` (via `main.cpp`, above) and `tise_solver` (via `tise_solver_main.cpp`, see [Building and Running `tise_solver`](#building-and-running-tise_solver)) are built on top of this same code, so a fix or feature added here benefits both binaries. A few pieces worth knowing about if you're reading or extending this file:

  - **`solveTISE(...)` and `SolveTISEResult`.** `solveTISE` is the top-level driver `main.cpp` calls: it builds a B-spline basis, fills and diagonalizes the Hamiltonian/overlap matrices, then constructs and writes the continuum states and phase shifts. It returns a `SolveTISEResult` bundling everything a caller needs to correctly interpret that solve afterward — the diagonalized eigenvalues/eigenvectors (`eigen`), the exact `BSpline` basis that was diagonalized against (`bs`), the exact physical grid passed to that basis (`grid`), the basis size (`nBSplines`), and the full set of physical B-spline indices excluded from the basis (`dropSet`). `main.cpp` now consumes these fields directly (for computing eigenstate coefficients, writing eigenstate files, running time evolution, etc.) instead of independently rebuilding its own grid/basis/drop-set — necessary because, per the next point, that basis isn't always the old fixed "uniform grid, drop B-spline 1" shape.

  - **Strategic node placement and singular-join B-spline removal**, done automatically and unconditionally by `buildStrategicGridAndDropSet(...)` on every solve (not something a caller opts into): it scans the supplied piecewise potential for structure at each piece-to-piece join. A join where `V` itself jumps ("Step") or where `V` is continuous but its slope `V'` jumps ("StitchedKink") gets extra degenerate B-spline knots inserted there, giving the basis just enough reduced smoothness to represent that discontinuity. A join where `V` actually *diverges* (e.g. an interior `1/x`-like term) is handled differently: instead of adding knots, the B-splines whose support touches that point are removed from the basis entirely. This removal only fires for **interior** joins — a singularity that coincides with the domain's own edge (e.g. hydrogen's Coulomb singularity sitting at `x = rMin`) is left to the classic single-B-spline wall exclusion instead (dropping B-spline 1, which already enforces `psi(rMin) = 0`); removing the whole edge-touching cluster there as well was tried and found to gut the basis exactly where a near-origin-peaked wavefunction has most of its amplitude (a 5x hydrogen ground-state error was observed empirically at B-spline order 8). A potential with no detectable structure at all gets back an unchanged uniform grid and the classic `{1}` drop-set, byte-for-byte identical to this project's original behavior.

  - **Shared between both binaries.** `buildStrategicGridAndDropSet` is called identically by `solveTISE` (the `H-BoundStates` path) and by `tise_solver_main.cpp` (the `tise_solver` path) — before this was unified, `tise_solver_main.cpp` built a plain uniform grid and hardcoded the classic `{1}` drop-set, structurally unable to reach either strategic node placement or singular-join removal. Now both entry points construct their grid and drop-set the same way and therefore behave consistently.

  - **Generalized drop-set.** The matrix fill and every continuum-construction function accept an explicit set of physical B-spline indices to exclude, not just the classic `{1}` (or `{1, N}`, once the last B-spline is folded in for boundary coupling) — needed now that strategic placement can add interior indices to that set. Omitting it reproduces the original hardcoded behavior bit-for-bit, so existing callers that don't pass one are unaffected.

  - **`minInterNodeGap` / `computeEAcc` / `warnIfContinuumExceedsEAcc`: a basis-accuracy ceiling check.** A bound or continuum state whose de Broglie wavelength is too short for the B-spline node spacing to resolve can't be represented accurately; `computeEAcc` estimates the energy above which that starts to matter, from the grid's node spacing. Because the grid can now be non-uniform (see strategic placement above), `minInterNodeGap` reduces it to the minimum spacing between *distinct* physical nodes first (collapsing the degenerate/repeated knots that strategic placement inserts, so the raw minimum doesn't come out as zero). **Current wiring:** this check is wired into `tise_solver`'s continuum path (`tise_solver_main.cpp`, gated on `tise.continuum.enabled`) — it calls `minInterNodeGap` on the grid actually used for that solve and warns into `warnings.json` if the configured `E_max` exceeds the resulting ceiling. It is *not* called anywhere inside `solveTISE` itself, so the `H-BoundStates` binary does not currently perform or report this check.

  - **Newly-configurable classifier thresholds.** A number of the numeric thresholds behind `classifyAsymptote`, `detectPotentialStructure`, and `classifySequenceConvergence` (sample counts, jump tolerances, geometric-probe ratios, and similar) are exposed as optional C++ parameters with documented defaults rather than hardcoded magic numbers. See `docs/planning/tise-classifier-configurability-exposure.md` for an analysis of which of these are actually worth threading further up into a `config.yaml`/CLI knob, versus which should stay reachable only from C++ code that calls these functions directly.

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

**What you need to install yourself:**

* A C++17 compiler — `g++ >= 7`, Clang, or MSVC (Visual Studio 2019 or newer)
* CMake >= 3.16
* Git
* Python 3.9+ (only for `controller.py`/`analysis.py`, not for the C++ build)

That is the whole list, on macOS, Linux and Windows alike.

**What CMake handles for you:** every C++ library below is looked for on your
system first and, if absent, downloaded and built automatically as part of the
configure step. No package manager, no admin rights, no platform-specific
instructions.

| Library | Used for |
|---|---|
| [Eigen](https://eigen.tuxfamily.org) | the generalized eigensolver, and time evolution |
| [nlohmann/json](https://github.com/nlohmann/json) | the JSON-encoded piecewise potential, and the `warnings.json` sidecar |
| [muparser](https://beltoforion.de/en/muparser/) | evaluating each potential piece's `function` expression at runtime |
| [yaml-cpp](https://github.com/jbeder/yaml-cpp) | parsing `config.yaml` (`tise_solver` only) |
| [GoogleTest](https://github.com/google/googletest) | the test suite (`-DBUILD_TESTING=ON` only) |

The first configure therefore needs a network connection and takes a couple of
minutes; later ones reuse what was fetched into the build directory.

If you would rather use system packages — they build faster, and this is what
CI's `playground-system-deps` job does — install them and CMake will prefer them:

```bash
# Debian/Ubuntu
sudo apt-get install libeigen3-dev nlohmann-json3-dev libmuparser-dev \
     libyaml-cpp-dev libgtest-dev

# macOS (Homebrew)
brew install eigen nlohmann-json muparser yaml-cpp googletest
```

Add `-DTISE_FETCH_DEPENDENCIES=OFF` to make a missing system package a hard
error instead of silently downloading it.

### Why there is no LAPACK here

The generalized eigenproblem `H c = E S c` is solved by Eigen's
`GeneralizedSelfAdjointEigenSolver`, in `tise.cpp`'s
`solveGeneralizedEigenproblem`. It used to be a call to LAPACK's banded
`dsbgv`, and the matrices are still assembled in LAPACK's upper banded layout.

LAPACK was dropped because it is Fortran: every route to it from source needs a
Fortran toolchain, which is exactly the barrier that made this project
unbuildable on Windows. Since that is the only place LAPACK was used, removing
it makes every remaining dependency header-only or a small CMake project.

The obvious objection is that Eigen has no *banded* generalized eigensolver, so
this throws away the band structure. Measured, that costs nothing here, because
the solver asks for **all** eigenvectors:

| N (basis size) | LAPACK `dsbgv`, banded | Eigen, dense |
|---:|---:|---:|
| 61 (default `config.yaml`) | 0.79 ms | 0.35 ms |
| 250 | 12.7 ms | 10.6 ms |
| 1000 | 549 ms | 577 ms |
| 2000 | 4986 ms | 4790 ms |

LAPACK does exploit the band, but only in the reduction to tridiagonal form.
Accumulating the full N x N eigenvector matrix afterwards is O(N^3) regardless
of how the input was stored, and that dominates. The same routine asked for
eigenvalues *only* (`jobz='N'`) takes 184 ms at N=2000 rather than 4986 ms --
i.e. 96% of its time is the band-agnostic part.

Where the band would genuinely win, and where this decision would need
revisiting: **memory** (dense is O(N^2), so N=20000 would need 3.2 GB per
matrix against 3 MB banded), and any future need for eigenvalues only or for a
subset of the spectrum (`dsbgvx`). Neither applies at the sizes this code is
used at, and `solveGeneralizedEigenproblem` is the single function that would
have to change.

Note that eigenvectors are sign-normalized by `solveGeneralizedEigenproblem`
itself rather than left as whichever sign the solver happened to produce (each
state is scaled so its first significant B-spline coefficient is positive, i.e.
psi_n > 0 just inside the left boundary). An eigenvector is only defined up to
a sign, and without pinning it the same `config.yaml` can plot some states
upside down on one machine versus another.

---

## Building with CMake

From this `TISE/` directory (there is no CMakeLists.txt one level up, so
`-S .` only works from here -- the top-level README's `-S TISE -B TISE/build`
is the equivalent invocation from the project root):

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

This creates the `H-BoundStates` and `tise_solver` executables in `build/`
(`build/Release/` if your generator is a multi-config one, which is the default
for Visual Studio and Xcode).

`--config Release` and `-DCMAKE_BUILD_TYPE=Release` are both given because
different generators honour different ones: single-config generators (Unix
Makefiles, Ninja) read `CMAKE_BUILD_TYPE` at configure time and ignore
`--config`, multi-config generators do the reverse. Passing both means the same
two commands work everywhere. Without them, a Visual Studio build silently
produces an unoptimized Debug binary.

### Windows

Nothing platform-specific is needed. The two commands above work in a
"Developer Command Prompt for VS" or any shell where MSVC is on `PATH`, and
CMake fetches and builds the C++ dependencies itself. Building inside WSL works
too, and is then just the Linux instructions.

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
   # $EIGEN is Eigen's include root: /usr/include/eigen3 on Debian/Ubuntu,
   # /opt/homebrew/include/eigen3 with Homebrew, /opt/local/include/eigen3 with
   # MacPorts. The CMake build finds this for you via find_package(Eigen3), and
   # downloads Eigen if it is not installed at all -- this manual recipe
   # assumes you have it.
   EIGEN=/usr/include/eigen3

   g++ -O2 -std=c++17 -c BSpline.cpp
   g++ -O2 -std=c++17 -I$EIGEN -c tise.cpp $(pkg-config --cflags muparser)
   g++ -O2 -std=c++17 -I$EIGEN -c time_evolution.cpp
   g++ -O2 -std=c++17 -I$EIGEN -c main.cpp
   ```

2. Link them with muparser (Eigen is header-only, and there is no LAPACK to
   link -- see [Why there is no LAPACK here](#why-there-is-no-lapack-here)):

   ```bash
   g++ -O2 -std=c++17 -o H-BoundStates main.o tise.o time_evolution.o BSpline.o \
       $(pkg-config --libs muparser)
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

* **Time evolution always runs.** After solving the bound-state problem, `main.cpp` unconditionally calls `tevol::runTimeEvolution(...)` — there is currently no way to request a TISE-only run from this binary. `config.yaml` already defines a `run.run_tdse` flag for exactly this purpose (see `docs/SDD.md` §6.1), but `main.cpp` has no YAML dependency at all and never reads `config.yaml`, so wiring that flag up here would mean giving this small demo driver a config-parsing capability it doesn't otherwise need — out of scope for `H-BoundStates` specifically, and left that way. `tise_solver` doesn't share this limitation: it never calls `runTimeEvolution` at all (time-dependent evolution stays entirely out of `tise_solver`), so `run.run_tdse` only affects whether `controller.py` invokes a separate TDSE step, not anything in this binary.

* **Continuum phase-shift matching assumes a regular right boundary.** If the supplied potential is singular at the domain's right edge (`x = rMax`), the two binaries diverge on purpose: `H-BoundStates`/`solveTISE` still attempts continuum construction and only prints a warning to stderr (an accepted limitation of this standalone demo binary), while `tise_solver` refuses continuum construction outright — no `phase_shifts.dat`/`continuum_state_NNN.dat` are written — and surfaces the diagnostic through `warnings.json` instead. Both binaries learn this from the same `StrategicGridResult::rightEdgeSingular` flag, set by the shared `buildStrategicGridAndDropSet` described in [Code Structure](#code-structure); `matchAsymptotic`'s flat-asymptote matching formula assumes the potential is negligible at `x = R`, which isn't true when the potential itself diverges there. `tise_solver` applies the identical refusal when a potential has a *genuine interior* singular join instead (e.g. `tests/interior_singularity.yaml`'s box-union-repulsive-Coulomb split) — `StrategicGridResult::interiorSingularSplit` — since matching at `x = rMax` has no meaning for a domain that's really two physically decoupled regions. Coulomb-tail continuum matching (for a potential whose tail is genuinely `1/r`-like beyond the box, as opposed to flat) **is implemented** — see [ADR-0013](../docs/adr/0013-coulomb-tail-continuum-matching.md) and [§3.6.1 of the input-file authoring guide](../docs/guides/input-file-authoring-guide.md#361-the-l-field-and-coulomb-tail-matching) for how to enable it via `tise.continuum.l`.

