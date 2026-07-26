"""Real-subprocess integration tests for the Controller<->TISE Solver
contract (docs/SDD.md Sec 7.2.1, Sec 9.2).

Unlike test_controller_unit.py -- which mocks subprocess.run and never
spawns a real process -- every test in this module drives
controller.run_tise_solver (or controller.run) against the REAL, compiled
tise_solver binary, reading and writing REAL files on disk. This is the
"real subprocess boundary (real binaries, real files)" that SDD Sec 9.2
calls for in addition to the pure-mock unit-test suite; unit tests with
mocked subprocess.run cannot prove the real CLI flags or file contract
actually line up between the two real programs.

CRITICAL test-hygiene note: the real config.yaml has run.output_dir:
"./data". NOTHING in this module ever runs the solver against that raw
config.yaml path -- every test goes through the `tmp_config` fixture
(tests/conftest.py), which copies config.yaml and redirects run.output_dir
into pytest's tmp_path, so a real solver invocation can never write into the
actual project's data/ directory.

Marked `integration` throughout (see pytest.ini) so this slower,
real-subprocess suite can be selected/excluded independently of the fast
mocked unit suite, e.g. `pytest -m "not integration"` to skip it.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from controller import SolverStageError, run, run_tise_solver

pytestmark = pytest.mark.integration


# ─── run_tise_solver against the real binary ────────────────────────────────


class TestRunTiseSolverRealSubprocess:
    """Both tests here invoke the REAL tise_solver binary (via the
    session-scoped tise_solver_binary fixture, built once per test session)
    -- controller.subprocess.run is never patched in this module."""

    def test_success_writes_expected_files_and_no_continuum_output(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        tise_dir = tmp_path / "data" / "tise"

        run_tise_solver(str(tmp_config), tise_dir, binary=tise_solver_binary)

        for name in ("eigenvalues.dat", "eigenvectors.dat", "hamiltonian.dat", "overlap.dat", "warnings.json"):
            assert (tise_dir / name).is_file(), f"expected output file missing: {name}"

        # The real config.yaml has tise.continuum.enabled: false -- the
        # continuum-only outputs must NOT be written.
        assert not (tise_dir / "phase_shifts.dat").exists()
        assert list(tise_dir.glob("continuum_state_*.dat")) == []

        warnings = json.loads((tise_dir / "warnings.json").read_text())
        assert warnings == []

    def test_failure_raises_solver_stage_error_and_leaves_no_partial_files(
        self, tise_solver_binary: Path, tmp_path: Path
    ):
        tise_dir = tmp_path / "data" / "tise"

        with pytest.raises(SolverStageError) as excinfo:
            run_tise_solver("/nonexistent/config.yaml", tise_dir, binary=tise_solver_binary)

        assert "TISE solver" in str(excinfo.value)

        # run_tise_solver's OWN mkdir() call runs unconditionally, before the
        # subprocess is even spawned, so tise_dir may legitimately exist --
        # but the real tise_solver binary fails at YAML::LoadFile (config
        # doesn't exist) before it creates its own output dir or writes any
        # file, so tise_dir must hold nothing regardless.
        assert not tise_dir.exists() or not any(tise_dir.iterdir())


# ─── Optional: full controller.run() orchestration, end to end ─────────────


class TestRunEndToEnd:
    """One extra end-to-end check beyond the two mandatory tests above: runs
    the full config-load -> validate -> subprocess -> warnings pipeline
    through controller.run() itself (as main() would), still against the
    real binary and real tmp_path files. controller.run() doesn't take a
    `binary=` override, so it resolves controller.DEFAULT_TISE_SOLVER --
    which happens to already point at TISE/build/tise_solver, the same path
    tise_solver_binary builds/locates. Requesting tise_solver_binary here
    still matters: it guarantees the binary is built before run() looks for
    it, independent of test execution order.
    """

    def test_run_end_to_end_success(self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path, capsys):
        run(str(tmp_config))

        tise_dir = tmp_path / "data" / "tise"
        assert (tise_dir / "eigenvalues.dat").is_file()
        warnings = json.loads((tise_dir / "warnings.json").read_text())
        assert warnings == []

        # read_warnings/print_warnings ran for real too: an empty warnings
        # list means nothing gets printed to stderr.
        captured = capsys.readouterr()
        assert captured.err == ""
