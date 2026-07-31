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

import subprocess
import sys
from pathlib import Path

import pytest
import yaml

from analysis import ContinuumPoint, EigenvalueRow, PhaseShiftRow, read_tise_output
from controller import run_tise_solver

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


# Mirrors TISE/tise_solver_main.cpp's own placeholder-dimension constants
# (kPlaceholderNumEigenstates, kPlaceholderBasisSize, kPlaceholderBandwidth,
# as of commit be3a403, confirmed by manually running the built binary
# before writing this file). C++ constants aren't importable from Python,
# so these are duplicated here as plain ints -- kept as named constants
# rather than inlined magic numbers so a future shape change in the stub
# reads as an obvious, deliberate two-sided edit rather than a silent drift.
# If the real binary's output shape ever changes without a matching update
# here, these tests fail loudly instead of the two sides quietly
# disagreeing.
_NUM_EIGENSTATES = 5  # kPlaceholderNumEigenstates: eigenvalues.dat rows; eigenvectors.dat cols
_BASIS_SIZE = 8  # kPlaceholderBasisSize: eigenvectors.dat rows; hamiltonian.dat/overlap.dat cols
_BANDWIDTH = 3  # kPlaceholderBandwidth: hamiltonian.dat/overlap.dat rows

# writeContinuumOutputs() calls writeIndexValueFile(..., 2) for every
# continuum_state_NNN.dat file -- 2 data rows is a fixed literal in the
# stub's own source, independent of n_energies (which only controls the
# NUMBER of such files, not each one's row count).
_CONTINUUM_STATE_POINTS_PER_FILE = 2


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
    """The real config.yaml's default (tise.continuum.enabled: false),
    unmodified, via the tmp_config fixture: run the real tise_solver binary,
    then read its real output back with analysis.read_tise_output, and
    confirm the two genuinely agree on shape."""

    def test_real_binary_output_matches_real_reader_shape(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        tise_dir = tmp_path / "data" / "tise"

        run_tise_solver(str(tmp_config), tise_dir, binary=tise_solver_binary)
        data = read_tise_output(tise_dir)

        # eigenvalues.dat: kPlaceholderNumEigenstates rows, 2 fields each
        # (index, E_n), index ascending from 0 -- checked per-row rather than
        # just by count, so a reader that scrambled order or mis-typed a
        # field would be caught here too.
        assert len(data.eigenvalues) == _NUM_EIGENSTATES
        for expected_index, row in enumerate(data.eigenvalues):
            assert isinstance(row, EigenvalueRow)
            assert row.index == expected_index
            assert isinstance(row.energy, float)

        # eigenvectors.dat: kPlaceholderBasisSize rows x kPlaceholderNumEigenstates cols.
        assert len(data.eigenvectors) == _BASIS_SIZE
        assert len(data.eigenvectors[0]) == _NUM_EIGENSTATES
        assert all(len(row) == _NUM_EIGENSTATES for row in data.eigenvectors)

        # hamiltonian.dat / overlap.dat: kPlaceholderBandwidth rows x kPlaceholderBasisSize cols.
        assert len(data.hamiltonian) == _BANDWIDTH
        assert len(data.hamiltonian[0]) == _BASIS_SIZE
        assert all(len(row) == _BASIS_SIZE for row in data.hamiltonian)

        assert len(data.overlap) == _BANDWIDTH
        assert len(data.overlap[0]) == _BASIS_SIZE
        assert all(len(row) == _BASIS_SIZE for row in data.overlap)

        # tise.continuum.enabled: false in the real config.yaml -- the real
        # binary writes neither phase_shifts.dat nor any continuum_state_*.dat
        # file, and Analysis must degrade to empty lists for both, not raise.
        assert data.phase_shifts == []
        assert data.continuum_states == []


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

        # One continuum_state_NNN.dat per n_energies; writeContinuumOutputs()
        # names them starting at i=1 (not 0), so ascending indices are
        # [1, 2, 3].
        assert len(data.continuum_states) == self.N_ENERGIES
        indices = [index for index, _points in data.continuum_states]
        assert indices == [1, 2, 3]

        for _index, points in data.continuum_states:
            assert isinstance(points, list)
            # writeIndexValueFile(..., 2): a fixed literal in
            # writeContinuumOutputs(), independent of n_energies.
            assert len(points) == _CONTINUUM_STATE_POINTS_PER_FILE
            assert all(isinstance(point, ContinuumPoint) for point in points)


# ─── Optional test 3: repeated run against the same tise_dir ───────────────


class TestRepeatedRunOverwritesConsistently:
    """Extra confidence check beyond the two mandatory scenarios above:
    running the solver a second time against the SAME tise_dir (so the
    second invocation overwrites the first run's files in place) still
    yields identical, correctly-shaped data on re-read. This guards against
    a regression class neither mandatory test above would catch on a single
    run -- e.g. a stub that opened its output files in append mode instead
    of truncating, which would only surface as doubled/corrupted rows after
    a SECOND run against the same directory. Cheap to include: the stub
    writes only placeholder content, so a second invocation costs one more
    fast subprocess call, not a second real solve."""

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
