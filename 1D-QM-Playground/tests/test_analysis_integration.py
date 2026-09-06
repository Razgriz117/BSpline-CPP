"""Real-subprocess integration tests for two related Analysis-side contracts
(docs/SDD.md Sec 9.2): the TISE Solver<->Analysis file-shape contract (Sec
7.2.2) and analysis.py's OWN Controller<->Analysis CLI/subprocess contract
(Sec 7.2.3).

The first three classes below (TestNoContinuumRealRoundTrip,
TestContinuumEnabledRealRoundTrip, TestRepeatedRunOverwritesConsistently)
cover Sec 7.2.2: unlike test_analysis_unit.py -- where every fixture is a
hand-written .dat file written directly into tmp_path, shaped by hand to
match the documented file formats, and the real tise_solver binary is never
built or run -- each drives controller.run_tise_solver against the REAL,
compiled tise_solver binary, then reads its REAL on-disk output with
analysis.read_tise_output (called in-process). Sec 7.2.2 itself has no
subprocess relationship between TISE and Analysis (only the Controller
invokes TISE as a subprocess; Analysis just reads files off disk), so what
these three classes prove is that the real tise_solver binary's actual
output and analysis.py's readers genuinely agree end to end on file shape --
NOT any unit-level parsing/error-path/tolerance logic, which
test_analysis_unit.py already covers exhaustively with hand-written
fixtures.

TestAnalysisCliRealSubprocess, below, instead covers Sec 7.2.3: it invokes
`python3.10 analysis.py --config ... --tise-dir ... --tdse-dir ...` as a raw
subprocess (subprocess.run([sys.executable, <path-to-analysis.py>, ...])) --
never calling read_tise_output() or anything else from analysis.py
in-process, and never importing or invoking controller.py except for
run_tise_solver, used purely to produce real data/tise/ fixture input
exactly as the three classes above already do. This proves analysis.py's
own CLI is "independently runnable and testable without the Controller in
the loop" (SDD Sec 2.4) -- distinct from test_controller_integration.py's
TestRunAnalysisStageRealSubprocess, which proves controller.py's correct
*usage* of this same Sec 7.2.3 contract (via controller.run_analysis_stage/
run(), which resolve controller.DEFAULT_ANALYSIS_SCRIPT internally), not
analysis.py's own CLI behavior under malformed/incomplete argv -- which only
TestAnalysisCliRealSubprocess covers.

CRITICAL test-hygiene note: exactly like test_controller_integration.py, the
real config.yaml has run.output_dir: "./data". NOTHING in this module ever
runs the solver against that raw config.yaml path -- every test goes
through the `tmp_config` fixture (tests/conftest.py), which copies
config.yaml and redirects run.output_dir into pytest's tmp_path, so a real
solver invocation can never write into the actual project's data/
directory.

Marked `integration` throughout (see pytest.ini) so this slower,
real-subprocess suite can be selected/excluded independently of the fast
mocked/fixture-based unit suites, e.g. `pytest -m "not integration"` to skip
it.
"""

from __future__ import annotations

import json
import math
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from analysis import ContinuumPoint, EigenstatePoint, EigenvalueRow, PhaseShiftRow, read_tise_output
from controller import run_tise_solver

# The tests/ dir itself -- where the known-solution reference configs
# (free_particle.yaml, finite_square_well.yaml, harmonic_oscillator.yaml,
# hydrogen.yaml, interior_singularity.yaml) live, next to this file.
_TESTS_DIR = Path(__file__).resolve().parent

pytestmark = pytest.mark.integration

# analysis.py's own script path, computed locally -- deliberately NOT
# imported from controller.DEFAULT_ANALYSIS_SCRIPT (see
# TestAnalysisCliRealSubprocess below), so this file's independence from
# controller.py stays literal for the one thing under test in that class:
# analysis.py's CLI itself.
_ANALYSIS_SCRIPT = Path(__file__).resolve().parent.parent / "analysis.py"


def _run_analysis_cli(config_path: Path | str, tise_dir: Path, tdse_dir: Path) -> subprocess.CompletedProcess:
    """Invoke analysis.py's own CLI as a raw subprocess with all three
    required flags, returning the CompletedProcess for the caller to assert
    on. Used by TestAnalysisCliRealSubprocess below for its three scenarios
    that supply all three flags; the one scenario that omits a flag
    entirely (argparse's own required-flag enforcement) builds its argv
    inline instead, since it doesn't fit this helper's fixed three-flag
    shape.
    """
    return subprocess.run(
        [
            sys.executable,
            str(_ANALYSIS_SCRIPT),
            "--config", str(config_path),
            "--tise-dir", str(tise_dir),
            "--tdse-dir", str(tdse_dir),
        ],
        capture_output=True,
        text=True,
    )


# Real shapes, derived from config.yaml's bspline.n_nodes=51, order=12 (C++
# constants aren't importable from Python, so these are duplicated here as
# plain ints -- kept as named constants rather than inlined magic numbers so
# a future config/shape change reads as an obvious, deliberate two-sided
# edit rather than a silent drift):
#   nBSplines = n_nodes + order - 2 = 51 + 12 - 2 = 61
#   nEn       = nBSplines - 2       = 59
# Per ADR-0007 (docs/adr/0007-defer-bound-state-filtering-tise-eigenvalue-
# output.md), eigenvalues.dat/eigenvectors.dat contain ALL nEn computed
# states -- no bound-state filtering is applied by tise_solver.
_NUM_EIGENSTATES = 59  # eigenvalues.dat rows; eigenvectors.dat cols; hamiltonian.dat/overlap.dat cols
_BASIS_SIZE = 61  # eigenvectors.dat rows
_BANDWIDTH = 12  # hamiltonian.dat/overlap.dat rows == bspline.order

# tise.continuum.n_pts (real config.yaml's default) -- now genuinely
# config-driven (tise::writeContinuumInfo's npts parameter), not a fixed
# stub literal. Tests below that enable continuum reuse the real config's
# n_pts unmodified, so this constant matches it.
_CONTINUUM_STATE_POINTS_PER_FILE = 500

# tise.n_pts_eigenstate (real config.yaml's default) -- eigenstate_NNN.dat
# points per file, written unconditionally regardless of continuum settings.
_EIGENSTATE_POINTS_PER_FILE = 301


def _write_no_continuum_config(tmp_config: Path, tmp_path: Path) -> Path:
    """Build a config variant with tise.continuum.enabled: false, regardless
    of the real config.yaml's own current value -- config.yaml has had
    continuum enabled by default since commit 99b7683, so
    TestNoContinuumRealRoundTrip below can no longer rely on tmp_config's
    unmodified copy to represent the disabled-continuum case. Same
    load/mutate/dump pattern as _write_continuum_enabled_config below."""
    with open(tmp_config) as f:
        cfg = yaml.safe_load(f)
    cfg["tise"]["continuum"]["enabled"] = False
    out = tmp_path / "config_no_continuum.yaml"
    with open(out, "w") as f:
        yaml.safe_dump(cfg, f)
    return out


def _write_continuum_enabled_config(tmp_config: Path, tmp_path: Path, n_energies: int) -> Path:
    """Build a config variant with tise.continuum.enabled: true and the
    given n_energies, following the exact load-YAML/mutate-a-copy/write-to-a-
    NEW-file pattern test_controller_integration.py's own bad-config tests
    use inline -- tmp_config's own file on disk is never opened for writing
    here, only read; the mutated dict is written to a fresh path under
    tmp_path instead."""
    with open(tmp_config) as f:
        cfg = yaml.safe_load(f)
    cfg["tise"]["continuum"]["enabled"] = True
    cfg["tise"]["continuum"]["n_energies"] = n_energies
    out = tmp_path / "config_continuum_enabled.yaml"
    with open(out, "w") as f:
        yaml.safe_dump(cfg, f)
    return out


# ─── Mandatory test 1: real config.yaml default (continuum disabled) ───────


class TestNoContinuumRealRoundTrip:
    """The real config.yaml's hydrogen-like potential with continuum
    explicitly disabled (see _write_no_continuum_config -- config.yaml's own
    default has had continuum enabled since commit 99b7683): run the real
    tise_solver binary, then read its real output back with
    analysis.read_tise_output, and confirm the two genuinely agree on
    shape."""

    def test_real_binary_output_matches_real_reader_shape(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _write_no_continuum_config(tmp_config, tmp_path)
        tise_dir = tmp_path / "data" / "tise"

        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        # eigenvalues.dat: _NUM_EIGENSTATES (=nEn) rows, 2 fields each
        # (index, E_n), index ascending from 0 -- checked per-row rather than
        # just by count, so a reader that scrambled order or mis-typed a
        # field would be caught here too.
        assert len(data.eigenvalues) == _NUM_EIGENSTATES
        for expected_index, row in enumerate(data.eigenvalues):
            assert isinstance(row, EigenvalueRow)
            assert row.index == expected_index
            assert isinstance(row.energy, float)

        # Physics-sanity checks (not just shape): EigenResult's own ascending
        # contract, and at least one genuine bound state (E<0) must exist for
        # this potential. The real config.yaml's potential, "-1/x + 1/x^2",
        # is a hydrogen L=1 effective potential (see config.yaml's own
        # comment) -- its lowest few eigenvalues have an exact analytic
        # answer, E_n = -1/(2*(n+1)^2) for n=1,2,3,..., giving -1/8, -1/18,
        # -1/32, independent of this project's B-spline numerics. Pinning
        # these catches gross regressions in the solve path itself, not just
        # in file I/O shape.
        energies = [row.energy for row in data.eigenvalues]
        assert energies == sorted(energies)
        assert energies[0] < 0
        assert energies[0] == pytest.approx(-1 / 8, abs=1e-6)
        assert energies[1] == pytest.approx(-1 / 18, abs=1e-6)
        assert energies[2] == pytest.approx(-1 / 32, abs=1e-6)

        # eigenvectors.dat: _BASIS_SIZE (=nBSplines) rows x _NUM_EIGENSTATES (=nEn) cols.
        assert len(data.eigenvectors) == _BASIS_SIZE
        assert len(data.eigenvectors[0]) == _NUM_EIGENSTATES
        assert all(len(row) == _NUM_EIGENSTATES for row in data.eigenvectors)

        # hamiltonian.dat / overlap.dat: _BANDWIDTH (=order) rows x
        # _NUM_EIGENSTATES (=nEn) cols -- NOT _BASIS_SIZE; nBSplines and nEn
        # are distinct quantities under real physics (the placeholder stub
        # coincidentally reused one constant for both).
        assert len(data.hamiltonian) == _BANDWIDTH
        assert len(data.hamiltonian[0]) == _NUM_EIGENSTATES
        assert all(len(row) == _NUM_EIGENSTATES for row in data.hamiltonian)

        assert len(data.overlap) == _BANDWIDTH
        assert len(data.overlap[0]) == _NUM_EIGENSTATES
        assert all(len(row) == _NUM_EIGENSTATES for row in data.overlap)

        # tise.continuum.enabled: false (forced by _write_no_continuum_config
        # above) -- the real binary writes neither phase_shifts.dat nor any
        # continuum_state_*.dat file, and Analysis must degrade to empty
        # lists for both, not raise.
        assert data.phase_shifts == []
        assert data.continuum_states == []

        # eigenstate_NNN.dat: written unconditionally (regardless of
        # tise.continuum.enabled, unlike continuum_state_NNN.dat above) --
        # one per bound state, each with n_pts_eigenstate (real config.yaml:
        # 301) points, indexed 1.._NUM_EIGENSTATES.
        assert len(data.eigenstates) == _NUM_EIGENSTATES
        indices = [index for index, _points in data.eigenstates]
        assert indices == list(range(1, _NUM_EIGENSTATES + 1))
        for _index, points in data.eigenstates:
            assert len(points) == _EIGENSTATE_POINTS_PER_FILE
            assert all(isinstance(point, EigenstatePoint) for point in points)


# ─── Mandatory test 2: continuum-enabled config variant ────────────────────


class TestContinuumEnabledRealRoundTrip:
    """A local config variant with tise.continuum.enabled: true and a small
    n_energies (3), built via _write_continuum_enabled_config() above --
    tmp_config's own file/fixture is reused read-only as the base, never
    mutated on disk."""

    N_ENERGIES = 3

    def test_real_binary_output_matches_real_reader_shape(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        continuum_config = _write_continuum_enabled_config(tmp_config, tmp_path, self.N_ENERGIES)
        tise_dir = tmp_path / "data" / "tise"

        run_tise_solver(str(continuum_config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        # phase_shifts.dat: one row per n_energies, 3 fields per row
        # (eps_i, delta_i, ddelta_dE_i).
        assert len(data.phase_shifts) == self.N_ENERGIES
        for row in data.phase_shifts:
            assert isinstance(row, PhaseShiftRow)
            assert len(row) == 3

        # One continuum_state_NNN.dat per n_energies; tise::writeContinuumInfo
        # names them starting at i=1 (not 0), so ascending indices are
        # [1, 2, 3].
        assert len(data.continuum_states) == self.N_ENERGIES
        indices = [index for index, _points in data.continuum_states]
        assert indices == [1, 2, 3]

        for _index, points in data.continuum_states:
            assert isinstance(points, list)
            # tise.continuum.n_pts (real config.yaml's default, unmodified
            # by _write_continuum_enabled_config) -- genuinely config-driven
            # now, not a fixed stub literal.
            assert len(points) == _CONTINUUM_STATE_POINTS_PER_FILE
            assert all(isinstance(point, ContinuumPoint) for point in points)

        # eigenstate_NNN.dat is written unconditionally -- confirms this
        # holds with continuum ALSO enabled, not just in the no-continuum
        # case TestNoContinuumRealRoundTrip above already covers in detail.
        assert len(data.eigenstates) == _NUM_EIGENSTATES


# ─── Optional test 3: repeated run against the same tise_dir ───────────────


class TestRepeatedRunOverwritesConsistently:
    """Extra confidence check beyond the two mandatory scenarios above:
    running the solver a second time against the SAME tise_dir (so the
    second invocation overwrites the first run's files in place) still
    yields identical, correctly-shaped data on re-read. This guards against
    a regression class neither mandatory test above would catch on a single
    run -- e.g. output files opened in append mode instead of truncating,
    which would only surface as doubled/corrupted rows after a SECOND run
    against the same directory. Also implicitly confirms the real solve is
    deterministic (same config in, byte-identical eigenvalues/eigenvectors
    out) across repeated invocations, not just idempotent file I/O."""

    def test_second_run_overwrites_first_with_identical_shape_and_values(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        continuum_config = _write_continuum_enabled_config(tmp_config, tmp_path, n_energies=3)
        tise_dir = tmp_path / "data" / "tise"

        run_tise_solver(str(continuum_config), tise_dir, binary=tise_solver_binary)
        first = read_tise_output(tise_dir)

        run_tise_solver(str(continuum_config), tise_dir, binary=tise_solver_binary)
        second = read_tise_output(tise_dir)

        assert second == first
        assert len(second.eigenvalues) == _NUM_EIGENSTATES
        assert len(second.phase_shifts) == 3
        assert len(second.continuum_states) == 3


# ─── Physics-correctness: free particle & harmonic oscillator vs. analytic ──
#
# Every quantity checked below has a closed-form analytic answer independent
# of this project's numerics, so these are true end-to-end regression tests,
# not shape/count checks. This is the test suite for the writeContinuumInfo
# coefficient-basis bug (see project memory continuum_state_coefficient_bug):
# a free particle's V=0 makes matchAsymptotic's flat-asymptote formula exact
# everywhere (not just asymptotically), so eigenvalues, phase shifts, AND
# the continuum wavefunction shape are all independently checkable.


def _write_free_particle_config(tmp_config: Path, tmp_path: Path, n_energies: int, e_max: float) -> Path:
    """Build a config variant with a V=0 potential spanning the real
    config's own bspline.domain unchanged, continuum enabled over
    [0, e_max]. Same load/mutate/dump pattern as
    _write_continuum_enabled_config above."""
    with open(tmp_config) as f:
        cfg = yaml.safe_load(f)
    domain = cfg["bspline"]["domain"]
    cfg["potential"] = [f"{{'domain': '[{domain[0]}, {domain[1]}]', 'function': '0'}}"]
    cfg["tise"]["continuum"]["enabled"] = True
    cfg["tise"]["continuum"]["E_threshold"] = 0.0
    cfg["tise"]["continuum"]["E_max"] = e_max
    cfg["tise"]["continuum"]["n_energies"] = n_energies
    out = tmp_path / "config_free_particle.yaml"
    with open(out, "w") as f:
        yaml.safe_dump(cfg, f)
    return out


def _write_harmonic_oscillator_config(tmp_config: Path, tmp_path: Path) -> Path:
    """Build a config variant for a 1D harmonic oscillator V(x)=0.5*omega^2*x^2
    on a symmetric box -- overrides bspline.n_nodes/domain (a genuinely
    different basis geometry from the real config's radial [0,100] hydrogen
    setup; n_nodes=81 over [-20,20] was verified during the 2026-08-28
    diagnosis session to match the first 10 analytic eigenvalues to ~3e-5)
    and disables continuum: a confining potential has no true continuum
    (docs/planning/architecture-06-20.md's Case-1 classification -- all
    states are box-confined pseudostates), so only bound states are
    physically meaningful to check here."""
    with open(tmp_config) as f:
        cfg = yaml.safe_load(f)
    cfg["bspline"]["n_nodes"] = 81
    cfg["bspline"]["domain"] = [-20.0, 20.0]
    cfg["potential"] = ["{'domain': '[-20, 20]', 'function': '0.5 * 1.0 * x^2'}"]
    cfg["tise"]["continuum"]["enabled"] = False
    out = tmp_path / "config_harmonic_oscillator.yaml"
    with open(out, "w") as f:
        yaml.safe_dump(cfg, f)
    return out


class TestFreeParticleContinuumPhysics:
    """A free particle (V=0) on the real config's own [0,100] box: bound
    eigenvalues are the infinite-square-well spectrum, the continuum phase
    shift is exactly zero (mod pi) at every energy, and the continuum
    wavefunction is an exact standing wave A*sin(kx) -- no confounding
    potential-shape assumptions anywhere, unlike the real config.yaml's
    Coulomb-tailed potential."""

    N_ENERGIES = 5
    E_MAX = 0.5  # comfortably inside the basis's own E_acc ceiling (~1.23 for
    # n_nodes=51/order=12/domain width=100, per the 2026-08-28 diagnosis run)

    def test_bound_eigenvalues_match_infinite_square_well(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _write_free_particle_config(tmp_config, tmp_path, self.N_ENERGIES, self.E_MAX)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        L = 100.0  # real config.yaml's bspline.domain width
        for row in data.eigenvalues[:20]:
            n = row.index + 1  # eigenvalues.dat is 0-indexed; box quantum number starts at 1
            expected = (n**2) * (math.pi**2) / (2 * L**2)
            assert row.energy == pytest.approx(expected, abs=1e-9)

    def test_phase_shift_matches_zero_scattering(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _write_free_particle_config(tmp_config, tmp_path, self.N_ENERGIES, self.E_MAX)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert len(data.phase_shifts) == self.N_ENERGIES
        for row in data.phase_shifts:
            # delta is only physically meaningful mod pi (matchAsymptotic's
            # own atan(...)-k*R construction bakes in an arbitrary -n*pi
            # winding), so a true delta=0 free particle satisfies
            # sin(delta) == 0 exactly, independent of winding number.
            assert math.sin(row.delta) == pytest.approx(0.0, abs=1e-3)

    def test_continuum_wavefunction_matches_analytic_sine(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _write_free_particle_config(tmp_config, tmp_path, self.N_ENERGIES, self.E_MAX)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert len(data.continuum_states) == self.N_ENERGIES
        for (index, points), phase_row in zip(data.continuum_states, data.phase_shifts):
            E = phase_row.energy
            k = math.sqrt(2 * E)
            amplitude = math.sqrt(2 / (math.pi * k))

            # A continuum state's overall sign is not physically meaningful
            # (the same freedom an eigenvector has, up to which sign the
            # underlying LAPACK solve happens to return) -- +sin(kx) and
            # -sin(kx) are equally valid. Pick the sign via the aggregate
            # correlation across every point (robust to any single point
            # landing near a node), then require every point to agree under
            # that ONE consistent sign -- a genuinely wrong shape (as under
            # the pre-fix coefficient-basis bug) cannot pass this regardless
            # of which sign is picked.
            correlation = sum(p.psi * amplitude * math.sin(k * p.x) for p in points)
            sign = 1.0 if correlation >= 0 else -1.0

            for point in points:
                expected = sign * amplitude * math.sin(k * point.x)
                assert point.psi == pytest.approx(expected, abs=5e-3), (
                    f"continuum state {index} at x={point.x}: got {point.psi}, expected {expected}"
                )
            # Dirichlet wall at rMin=0: structurally forced to exactly zero
            # (physical B-spline #1 excluded from the basis), not merely a
            # coincidental zero of sin(kx) -- worth its own assertion.
            assert points[0].psi == pytest.approx(0.0, abs=1e-9)


class TestHarmonicOscillatorBoundStates:
    """A 1D harmonic oscillator has no true continuum (Case 1: confining
    potential, box-confined pseudostates only -- confirmed during the
    2026-08-28 diagnosis session), so only its closed-form bound-state
    spectrum E_n=(n+1/2)*omega is physically meaningful to validate here."""

    def test_bound_eigenvalues_match_harmonic_oscillator_spectrum(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _write_harmonic_oscillator_config(tmp_config, tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        omega = 1.0
        for row in data.eigenvalues[:10]:
            n = row.index  # eigenvalues.dat is 0-indexed; HO quantum number also starts at 0
            expected = (n + 0.5) * omega
            assert row.energy == pytest.approx(expected, abs=1e-3)


# ─── Known-solution reference configs, loaded by name from tests/*.yaml ────
#
# tests/{free_particle,finite_square_well,harmonic_oscillator,hydrogen,
# interior_singularity}.yaml previously existed only as manual-inspection
# fixtures -- never loaded by any automated test (the classes above build
# equivalent configs inline instead). This section drives the real
# tise_solver binary against the ACTUAL committed files (only output_dir
# redirected into tmp_path), closing that gap and doubling as the
# regression backstop for Part A's solveTISE/tise_solver orchestration
# unification (docs/planning/tise-release-readiness-plan.md).


def _load_known_solution_config(name: str, tmp_path: Path) -> Path:
    """Load tests/<name>.yaml (the actual committed file, not an
    inline-built equivalent), redirect run.output_dir into tmp_path, and
    write the result to a fresh path under tmp_path -- the source file
    itself is never opened for writing."""
    with open(_TESTS_DIR / f"{name}.yaml") as f:
        cfg = yaml.safe_load(f)
    cfg["run"]["output_dir"] = str(tmp_path / "data")
    out = tmp_path / f"{name}_loaded.yaml"
    with open(out, "w") as f:
        yaml.safe_dump(cfg, f)
    return out


class TestKnownSolutionConfigFilesLoadAndRun:
    """free_particle.yaml and harmonic_oscillator.yaml already have rich
    physics assertions above (via inline-built equivalent configs) -- this
    just confirms the actual committed FILES load and run correctly too,
    closing the "orphaned fixture" gap without duplicating those checks."""

    def test_free_particle_yaml_runs_and_produces_expected_output(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _load_known_solution_config("free_particle", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert len(data.eigenvalues) > 0
        assert len(data.continuum_states) == 5  # free_particle.yaml's n_energies
        assert len(data.eigenstates) == len(data.eigenvalues)

    def test_harmonic_oscillator_yaml_runs_and_produces_expected_output(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _load_known_solution_config("harmonic_oscillator", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert data.continuum_states == []  # continuum disabled -- no true continuum
        omega = 1.0
        for row in data.eigenvalues[:10]:
            expected = (row.index + 0.5) * omega
            assert row.energy == pytest.approx(expected, abs=1e-3)


class TestFiniteSquareWellRealSubprocess:
    """finite_square_well.yaml had ZERO test coverage before this -- first
    real-subprocess physics validation. V0=1.0, a=10.0 attractive well;
    closed-form s-wave phase shift (same formula TISE/tests/test_tise.cpp's
    SquareWellPhaseShiftTest validates matchAsymptotic against directly),
    exercising genuinely non-trivial scattering unlike free particle's
    flat delta=0."""

    V0 = 1.0
    A = 10.0

    @staticmethod
    def _analytic_phase_shift(E: float, V0: float, a: float) -> float:
        k = math.sqrt(2 * E)
        kappa = math.sqrt(2 * (E + V0))
        return -k * a + math.atan((k / kappa) * math.tan(kappa * a))

    @staticmethod
    def _wrap_mod_pi(delta: float) -> float:
        return delta - math.pi * round(delta / math.pi)

    def test_has_genuine_bound_states(self, tise_solver_binary: Path, tmp_path: Path):
        config = _load_known_solution_config("finite_square_well", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        bound = [row.energy for row in data.eigenvalues if row.energy < 0]
        assert len(bound) > 0

    def test_phase_shift_matches_analytic_square_well_formula(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _load_known_solution_config("finite_square_well", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert len(data.phase_shifts) > 0
        for row in data.phase_shifts:
            expected = self._analytic_phase_shift(row.energy, self.V0, self.A)
            got = self._wrap_mod_pi(row.delta)
            want = self._wrap_mod_pi(expected)
            # Tightened from the original 5e-2: after the strategic-knot fix
            # this basis actually achieves <=2.4e-4 rad
            # (docs/tests/reports/8236239/finite_square_well.md) -- the
            # looser bound didn't discriminate the pre-fix ~0.03-0.6 rad
            # error from the current accuracy.
            assert got == pytest.approx(want, abs=5e-4)

    def test_ddelta_dE_matches_analytic_derivative(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        # Closes H2 (docs/planning/tise-known-solution-followup-plan.md):
        # dDeltaDE previously finite-differenced sin(2*delta) across the
        # production grid and divided by 2*cos(2*delta), giving e.g. +2.4
        # where the true slope is -19.4. The reference here is a fine
        # central difference of the SAME closed-form phase-shift formula
        # (not reusing any solver code), independent of the solver's own
        # internal fine-step mechanism.
        config = _load_known_solution_config("finite_square_well", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert len(data.phase_shifts) > 0
        h = 1e-5
        for row in data.phase_shifts:
            expected = (
                self._analytic_phase_shift(row.energy + h, self.V0, self.A)
                - self._analytic_phase_shift(row.energy - h, self.V0, self.A)
            ) / (2 * h)
            # abs=0.05 comfortably covers the largest observed deviation
            # (~0.023, at E=0.5 -- the grid point closest to a confined
            # eigenvalue, a known ADR-0007 box-discretization effect) while
            # remaining many orders of magnitude tighter than the old
            # formula's errors (which ran to double digits).
            assert row.ddelta_dE == pytest.approx(expected, abs=0.05)


class TestHydrogenConfigFileRealSubprocess:
    """hydrogen.yaml: same potential as the real config.yaml, but with a
    conservative E_max (inside this basis's own reported E_acc reliability
    ceiling, unlike the real config.yaml's E_max=10) and
    visualization.eigenstates enabled -- a genuinely different, complementary
    scenario from TestNoContinuumRealRoundTrip above, not a duplicate."""

    def test_bound_eigenvalues_match_analytic_hydrogen(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _load_known_solution_config("hydrogen", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        energies = [row.energy for row in data.eigenvalues]
        assert energies[0] == pytest.approx(-1 / 8, abs=1e-6)
        assert energies[1] == pytest.approx(-1 / 18, abs=1e-6)
        assert energies[2] == pytest.approx(-1 / 32, abs=1e-6)

    def test_phase_shift_is_approximately_zero_for_pure_coulomb_tail(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        # hydrogen.yaml's potential (-1/x + 1/x^2, l=1) is an EXACT pure
        # radial Coulomb equation -- the 1/x^2 term is exactly the l=1
        # centrifugal term l(l+1)/2=1, not an extra short-range correction.
        # A phase shift measures the DEVIATION from pure Coulomb
        # scattering; there is none here, so delta should be ~0 (mod pi).
        # This closes the Coulomb-tail-continuum-matching gap (ADR-0009
        # supersedes ADR-0010): the potential's domain is now declared
        # unbounded ('(0, inf)', not '(0, 100]') so classifyAsymptote
        # actually recognizes the Coulomb tail, and tise.continuum.l: 1 is
        # set so the matching uses the correct Coulomb wave functions.
        # Before this fix, the flat-asymptote formula gave delta of order
        # 1 rad here (docs/tests/reports/8236239/hydrogen.md) -- wrong by
        # the entire physical phase, not a small approximation error.
        config = _load_known_solution_config("hydrogen", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert len(data.phase_shifts) == 5
        for row in data.phase_shifts:
            want = 0.0
            got = row.delta - math.pi * round((row.delta - want) / math.pi)
            # abs=3e-3: comfortably above the largest observed residual
            # (~8e-4 at E=0.1) and consistent with the Coulomb-matching
            # implementation's own validated worst-case accuracy
            # (~2.5e-3 rad, docs/planning/coulomb-tail-continuum-matching.md)
            # -- box-truncation/numerical residual, not exactly zero, but
            # many orders of magnitude tighter than the pre-fix ~1 rad error.
            assert got == pytest.approx(want, abs=3e-3)

    def test_ddelta_dE_no_longer_diverges(self, tise_solver_binary: Path, tmp_path: Path):
        # Regression pin for H2 (docs/planning/tise-known-solution-followup-
        # plan.md): the old sin(2*delta)/cos(2*delta) formula gave +211.4 at
        # E=0.1 (docs/tests/reports/8236239/hydrogen.md), from dividing by
        # cos(2*delta)~=-0.044 -- an inherent 0/0 conditioning problem at
        # delta=pi/4 (mod pi/2), not fixable by a finer step. The fixed
        # fine-local-difference formula gives ~-19.8 there.
        config = _load_known_solution_config("hydrogen", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert len(data.phase_shifts) > 0
        for row in data.phase_shifts:
            assert abs(row.ddelta_dE) < 50


class TestInteriorSingularityRealSubprocess:
    """interior_singularity.yaml: a genuine INTERIOR Singular potential join
    (not at a domain edge like hydrogen's), run through the real tise_solver
    binary. Before Part A of the original release-readiness plan,
    tise_solver_main.cpp had no A4b singular-B-spline-removal path at all
    and would have silently produced degraded eigenvalues with EXIT_SUCCESS
    and zero warning. Before the known-solution-verification follow-up
    plan's Part A, the A4b removal that WAS wired dropped every B-spline
    touching the singular point on a simple-knot grid -- forcing a
    multi-bohr dead zone and producing eigenvalues 41% too high on the
    field-free side (docs/tests/reports/8236239/interior_singularity.md).
    Fixed via knot multiplicity + single-B-spline drop, exactly mirroring
    domain-edge treatment; the exact reference spectrum below is the split-
    domain union of box states on [0,20] and repulsive-Coulomb states on
    [20,40] (report section 3)."""

    def test_all_eigenvalues_and_eigenvectors_finite(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _load_known_solution_config("interior_singularity", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert len(data.eigenvalues) > 0
        for row in data.eigenvalues:
            assert math.isfinite(row.energy)
        for evec_row in data.eigenvectors:
            for value in evec_row:
                assert math.isfinite(value)

    def test_bspline_removal_uses_knot_multiplicity_not_classic_dropset(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        # nBSplines = n_nodes+order-2 = 41+8-2 = 47 on the OLD simple-knot
        # grid; the classic dropset ({1} only) would give nEnBound = 45. The
        # knot-multiplicity fix grows the basis by 6 (one interior join,
        # extraMultiplicity=order-2=6) and drops only 2 physical indices
        # (left edge + the single spline non-zero at the join), giving
        # nEnBound = (47+6)-2-1 = 50 -- the concrete number that would catch
        # either a regression back to the pre-fix cluster-drop (36) or an
        # accidental fall-back to the classic dropset entirely (45).
        config = _load_known_solution_config("interior_singularity", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert len(data.eigenvalues) == 50

    def test_eigenvalues_match_exact_split_domain_spectrum(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _load_known_solution_config("interior_singularity", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        energies = [row.energy for row in data.eigenvalues]
        assert energies[0] == pytest.approx(0.0123370055, abs=1e-9)  # box n=1
        assert energies[1] == pytest.approx(0.0493480220, abs=1e-9)  # box n=2
        assert energies[2] == pytest.approx(0.1005405216, abs=1e-8)  # Coulomb root 1


class TestRightEdgeSingularWarningRealSubprocess:
    """right_edge_singularity.yaml: potential singular exactly AT x=rMax
    (1/(100-x)) -- distinct from hydrogen's LEFT-edge singularity. Before
    Part A of the original release-readiness plan, tise_solver_main.cpp had
    no access to this diagnostic at all (only solveTISE/H-BoundStates, via
    a bare stderr print, could detect it). Before the known-solution-
    verification follow-up plan's Part D, continuum construction completed
    anyway and wrote phase_shifts.dat/continuum_state_NNN.dat even though
    the result is not a solution of the problem (psi_E(rMax) must vanish
    physically; the construction instead produces essentially B_N alone,
    30x larger at the wall than anywhere in the interior --
    docs/tests/reports/8236239/right_edge_singularity.md section 4.3).
    Continuum construction is now refused entirely for this config."""

    def test_warnings_json_contains_right_edge_singular_warning(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _load_known_solution_config("right_edge_singularity", tmp_path)
        tise_dir = tmp_path / "data" / "tise"

        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        warnings = json.loads((tise_dir / "warnings.json").read_text())

        messages = [w["message"] for w in warnings]
        assert any("singular at the right domain edge" in m for m in messages)

    def test_bound_state_solve_is_unaffected(self, tise_solver_binary: Path, tmp_path: Path):
        # The refusal is scoped to continuum construction only -- the bound
        # states (17 below E=0.2, per the report) are still solved and
        # written normally.
        config = _load_known_solution_config("right_edge_singularity", tmp_path)
        tise_dir = tmp_path / "data" / "tise"

        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert len(data.eigenvalues) > 0
        assert data.eigenvalues[0].energy == pytest.approx(0.0149922127, abs=1e-9)

    def test_no_continuum_output_files_are_written(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _load_known_solution_config("right_edge_singularity", tmp_path)
        tise_dir = tmp_path / "data" / "tise"

        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert not (tise_dir / "phase_shifts.dat").exists()
        assert not list(tise_dir.glob("continuum_state_*.dat"))
        assert data.phase_shifts == []
        assert data.continuum_states == []


class TestBoundStateDiagnosticsRealSubprocess:
    """classifyBoundStates/checkWellContainment (Part B): previously fully
    implemented and unit-tested but never called from tise_solver_main.cpp
    -- production output never labeled the bound/continuum-artifact split
    or flagged a wall-colliding bound state. The real config.yaml's
    hydrogen potential in its R=100 box is a naturally-occurring case for
    both: 7 of 59 states are below E=0 (per ADR-0007's own count), and its
    3 highest-n bound states (larger radial extent) are only marginally
    contained by R=100."""

    def test_warnings_json_reports_bound_state_count(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(tmp_config), tise_dir, binary=tise_solver_binary)
        warnings = json.loads((tise_dir / "warnings.json").read_text())

        messages = [w["message"] for w in warnings]
        assert any("7 of 59 computed states are below E=0.0" in m for m in messages)

    def test_warnings_json_flags_wall_colliding_bound_states_only(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(tmp_config), tise_dir, binary=tise_solver_binary)
        warnings = json.loads((tise_dir / "warnings.json").read_text())

        messages = [w["message"] for w in warnings]
        colliding = [m for m in messages if "colliding with the outer wall" in m]
        # States 5, 6, 7 (E<0, larger radial extent) are only marginally
        # contained by R=100; states 1-4 are not, and no unbound/continuum-
        # like state (index 8+, E>=0) should ever be flagged -- that check
        # is deliberately restricted to E<0 states (see tise_solver_main.cpp).
        assert any("bound state 5 " in m for m in colliding)
        assert any("bound state 6 " in m for m in colliding)
        assert any("bound state 7 " in m for m in colliding)
        assert not any("bound state 1 " in m for m in colliding)
        assert not any("bound state 8 " in m for m in colliding)


def _load_case3_config_with_continuum_enabled(tmp_path: Path) -> Path:
    """case3_irregular_tail.yaml as committed, but with tise.continuum
    enabled -- exercises the "taper actually applied" branch (the
    committed file itself has continuum disabled, exercising the
    "taper gated off" branch instead -- see
    docs/tests/reports/8236239/case3_irregular_tail.md recommendation 1)."""
    with open(_TESTS_DIR / "case3_irregular_tail.yaml") as f:
        cfg = yaml.safe_load(f)
    cfg["run"]["output_dir"] = str(tmp_path / "data")
    cfg["tise"]["continuum"] = {
        "enabled": True, "E_threshold": 0.0, "E_max": 0.1, "n_energies": 2, "n_pts": 50,
    }
    out = tmp_path / "case3_irregular_tail_continuum_enabled.yaml"
    with open(out, "w") as f:
        yaml.safe_dump(cfg, f)
    return out


class TestCase3RemediationRealSubprocess:
    """case3_irregular_tail.yaml: "1/x^1.5" on (0.1,inf), same
    proven Case-3-triggering potential TISE/tests/test_tise.cpp's
    ClassifyAsymptoteTest.Case3IrregularForPowerLawOneAndHalf and
    FillBandedMatricesCase3WindowTest validate directly. Before Part B of
    the original release-readiness plan, classifyAsymptote could DETECT
    this and warn, but nothing ever tapered the potential -- detection
    without remediation. Before the known-solution-verification follow-up
    plan's Part C, the taper (once wired) was applied UNCONDITIONALLY, even
    with continuum disabled (as this committed YAML has it) -- shifting
    bound-state energies by 1.5e-5 to 1.5e-4 Ha for zero benefit, since the
    taper exists solely to satisfy matchAsymptotic's flat-asymptote
    assumption for continuum matching
    (docs/tests/reports/8236239/case3_irregular_tail.md)."""

    def test_all_eigenvalues_finite(self, tise_solver_binary: Path, tmp_path: Path):
        config = _load_known_solution_config("case3_irregular_tail", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert len(data.eigenvalues) > 0
        for row in data.eigenvalues:
            assert math.isfinite(row.energy)

    def test_ground_state_matches_untapered_reference_when_continuum_disabled(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        # Richardson-extrapolated finite-difference reference for the RAW
        # (untapered) 1/x^1.5 potential, report section 4, E^raw column,
        # state 1. With continuum disabled (as committed), no taper should
        # be applied at all -- this is the number that would have caught
        # the pre-fix unconditional-taper bug (which gave 0.0092216090127,
        # 1.5e-5 lower).
        config = _load_known_solution_config("case3_irregular_tail", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        assert data.eigenvalues[0].energy == pytest.approx(0.0092368363, abs=1e-7)

    def test_warnings_json_reports_no_taper_when_continuum_disabled(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _load_known_solution_config("case3_irregular_tail", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        warnings = json.loads((tise_dir / "warnings.json").read_text())

        messages = [w["message"] for w in warnings]
        assert any("Irregular" in m and "tapering" not in m for m in messages)

    def test_warnings_json_reports_tapering_when_continuum_enabled(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _load_case3_config_with_continuum_enabled(tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        warnings = json.loads((tise_dir / "warnings.json").read_text())

        messages = [w["message"] for w in warnings]
        assert any("Irregular" in m and "tapering" in m for m in messages)


def _write_flat_tail_config(tmp_config: Path, tmp_path: Path, tail_value: float) -> Path:
    """A well on [0,5) plus a CONSTANT (hence trivially Flat-classified --
    zero successive differences, no power-law-fit tolerance to reason
    about) tail of `tail_value` on [5, inf), box domain [0,20]. Continuum
    enabled with E_max=0.1/n_energies=2, so the smallest requested energy
    is 0.05 -- deliberately chosen so tail_value=0.05 gives ratio 1.0 (Part
    B's R-validity check should fire) and tail_value=0.0001 gives ratio
    0.002 (should not). The tail must extend to infinity, not stop exactly
    at rMax=20, or tise_solver_main.cpp's own "is there really an unbounded
    side to classify" probe skips classification entirely."""
    with open(tmp_config) as f:
        cfg = yaml.safe_load(f)
    cfg["bspline"]["n_nodes"] = 41
    cfg["bspline"]["order"] = 8
    cfg["bspline"]["domain"] = [0.0, 20.0]
    cfg["potential"] = [
        "{'domain': '[0, 5)', 'function': '-1.0'}",
        f"{{'domain': '[5, inf)', 'function': '{tail_value}'}}",
    ]
    cfg["tise"]["continuum"] = {
        "enabled": True, "E_threshold": 0.0, "E_max": 0.1, "n_energies": 2, "n_pts": 50,
    }
    out = tmp_path / f"config_flat_tail_{tail_value}.yaml"
    with open(out, "w") as f:
        yaml.safe_dump(cfg, f)
    return out


class TestRValidityWarningRealSubprocess:
    """Part B of the Coulomb-tail/boundary-validation follow-up plan:
    classifyAsymptote only checks the tail's functional SHAPE far past
    rMax (geometrically outward), never whether V has already become
    small AT rMax specifically relative to the requested continuum
    energies -- previously nothing warned a user whose box was too small
    for their energy range. Scoped to the Flat sub-case only (a Coulomb
    tail's V(rMax) is expected to still be non-negligible; that's a
    separate, not-yet-implemented concern -- see ADR-0010)."""

    def test_warnings_json_flags_non_negligible_v_at_rmax(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _write_flat_tail_config(tmp_config, tmp_path, tail_value=0.05)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        warnings = json.loads((tise_dir / "warnings.json").read_text())

        messages = [w["message"] for w in warnings]
        assert any("not negligible" in m and "V(rMax)" in m for m in messages)

    def test_warnings_json_does_not_flag_a_negligible_v_at_rmax(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        config = _write_flat_tail_config(tmp_config, tmp_path, tail_value=0.0001)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        warnings = json.loads((tise_dir / "warnings.json").read_text())

        messages = [w["message"] for w in warnings]
        assert not any("not negligible" in m and "V(rMax)" in m for m in messages)

    def test_hydrogen_yaml_coulomb_tail_is_not_flagged_by_this_check(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        # hydrogen.yaml's tail is Coulomb (V(100)=-0.01, genuinely
        # non-negligible against E=0.1 -- see docs/tests/reports/9d0f04b/
        # hydrogen.md), but this check is deliberately scoped to Flat only;
        # a Coulomb tail needs Coulomb-aware matching (ADR-0010), not a
        # bigger box, so it must not trigger this particular warning.
        config = _load_known_solution_config("hydrogen", tmp_path)
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(config), tise_dir, binary=tise_solver_binary)
        warnings = json.loads((tise_dir / "warnings.json").read_text())

        messages = [w["message"] for w in warnings]
        assert not any("not negligible" in m and "V(rMax)" in m for m in messages)


# ─── analysis.py's OWN CLI/subprocess contract (Sec 7.2.3), invoked directly ─


class TestAnalysisCliRealSubprocess:
    """Four scenarios invoking analysis.py's CLI directly as a raw
    subprocess: success with a real data/tise/ directory and a missing
    --tdse-dir (exit 0), a missing/unpopulated --tise-dir (exit 1), a
    missing --config (exit 1), and a missing required flag -- argparse's
    own enforcement (exit 2). See this file's module docstring above for
    why this class exists and how it differs from
    test_controller_integration.py's TestRunAnalysisStageRealSubprocess.
    """

    def test_success_missing_tdse_dir_exits_zero_with_info_note_on_stderr(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        # Build a real data/tise/ directory first (exactly as
        # TestNoContinuumRealRoundTrip does above) -- Analysis's required
        # input (Sec 7.2.2). run_tise_solver is Controller-side plumbing
        # used only to produce this fixture; the subprocess actually under
        # test below is analysis.py itself, invoked directly.
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(tmp_config), tise_dir, binary=tise_solver_binary)

        # tdse_dir is deliberately a path that was never created -- Sec
        # 5.4.4's tolerated case (run_tdse: false is expected until Phase
        # 8). Proves that tolerance survives the real subprocess/CLI
        # boundary, invoked with no controller.py anywhere in the loop.
        tdse_dir = tmp_path / "data" / "tdse"

        result = _run_analysis_cli(tmp_config, tise_dir, tdse_dir)

        assert result.returncode == 0
        assert result.stdout == ""
        assert "info" in result.stderr
        assert "--tdse-dir" in result.stderr
        assert "not found or not a directory" in result.stderr

    def test_missing_tise_dir_exits_one_with_eigenvalues_error_on_stderr(
        self, tmp_config: Path, tmp_path: Path
    ):
        # tise_dir deliberately never created/populated -- read_tise_output's
        # first call, read_eigenvalues(), fails the moment it can't read
        # eigenvalues.dat (Sec 7.2.2: the four core files are non-optional),
        # run() lets the resulting TiseOutputError propagate, and main()
        # catches it (AnalysisError) and exits 1 -- proven here directly
        # against the CLI, without controller.py or its
        # run_analysis_stage/SolverStageError wrapping anywhere in the loop.
        tise_dir = tmp_path / "data" / "tise"
        tdse_dir = tmp_path / "data" / "tdse"

        result = _run_analysis_cli(tmp_config, tise_dir, tdse_dir)

        assert result.returncode == 1
        assert result.stdout == ""
        assert "eigenvalues" in result.stderr

    def test_missing_config_exits_one_with_config_error_on_stderr(self, tmp_path: Path):
        # A --config path that was never created -- exercises analysis.py's
        # own load_config()/ConfigError path directly, deliberately NOT
        # using the tmp_config fixture (which always points at a real,
        # valid file). tise_dir/tdse_dir are never created either, but
        # load_config() runs first inside run(), so this never reaches
        # read_tise_output().
        nonexistent_config = tmp_path / "nonexistent_config.yaml"
        tise_dir = tmp_path / "data" / "tise"
        tdse_dir = tmp_path / "data" / "tdse"

        result = _run_analysis_cli(nonexistent_config, tise_dir, tdse_dir)

        assert result.returncode == 1
        assert result.stdout == ""
        assert "not found" in result.stderr

    def test_missing_required_flag_exits_two_via_argparse(self, tmp_config: Path, tmp_path: Path):
        # --tdse-dir omitted entirely. This is argparse's OWN
        # required-argument enforcement (all three of analysis.py main()'s
        # parser.add_argument(..., required=True) flags), firing before
        # analysis.run() -- and therefore before load_config()/
        # read_tise_output() -- is ever reached. Not something to change,
        # only to confirm: exit code 2 and argparse's own usage/error text
        # are Python's standard argparse behavior, not analysis.py's own
        # error-handling code path.
        tise_dir = tmp_path / "data" / "tise"

        result = subprocess.run(
            [
                sys.executable,
                str(_ANALYSIS_SCRIPT),
                "--config", str(tmp_config),
                "--tise-dir", str(tise_dir),
                # --tdse-dir deliberately omitted
            ],
            capture_output=True,
            text=True,
        )

        assert result.returncode == 2
        assert result.stdout == ""
        assert "usage" in result.stderr
        assert "--tdse-dir" in result.stderr
