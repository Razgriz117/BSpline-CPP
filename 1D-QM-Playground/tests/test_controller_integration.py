"""Real-subprocess integration tests for controller.py's own use of two
Controller-side inter-component contracts: Controller<->TISE Solver
(docs/SDD.md Sec 7.2.1) and Controller<->Analysis (Sec 7.2.3), per Sec 9.2.

Unlike test_controller_unit.py -- which mocks subprocess.run and never
spawns a real process -- every test in this module drives
controller.run_tise_solver, controller.run_analysis_stage, or controller.run
against the REAL, compiled tise_solver binary and/or the REAL analysis.py
script, reading and writing REAL files on disk. This is the "real subprocess
boundary (real binaries, real files)" that SDD Sec 9.2 calls for in addition
to the pure-mock unit-test suite; unit tests with mocked subprocess.run
cannot prove the real CLI flags or file contract actually line up between
the real programs involved. (analysis.py's OWN CLI behavior, independent of
controller.py, is proven separately by test_analysis_integration.py -- this
module only proves controller.py's correct usage of each contract.)

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
import yaml

from controller import SolverStageError, run, run_analysis_stage, run_tise_solver

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

    def test_continuum_enabled_missing_n_energies_raises_and_leaves_no_partial_files(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        """Regression guard for the fix validating tise.continuum.n_energies
        before any output is written (docs/SDD.md Sec 7.2.1).

        main() used to call writeCoreOutputs() unconditionally and only
        validate tise.continuum.n_energies afterwards, while extracting the
        continuum settings -- so continuum.enabled: true with a missing (or
        otherwise invalid) n_energies exited non-zero but still left the 5
        core files (eigenvalues.dat, etc.) sitting in the output directory.
        n_energies is now validated before writeCoreOutputs() runs, so this
        config must both raise SolverStageError AND leave tise_dir holding
        zero files -- the same "non-zero exit == no partial output"
        guarantee the config-missing test above checks for a different
        failure path.
        """
        with open(tmp_config) as f:
            cfg = yaml.safe_load(f)
        cfg["tise"]["continuum"]["enabled"] = True
        cfg["tise"]["continuum"].pop("n_energies", None)
        bad_config = tmp_path / "config_continuum_missing_n_energies.yaml"
        with open(bad_config, "w") as f:
            yaml.safe_dump(cfg, f)

        tise_dir = tmp_path / "data" / "tise"

        with pytest.raises(SolverStageError) as excinfo:
            run_tise_solver(str(bad_config), tise_dir, binary=tise_solver_binary)

        assert "TISE solver" in str(excinfo.value)

        # This is the actual regression guard: before the fix, main() wrote
        # the 5 core files unconditionally, before ever validating
        # n_energies, so this assertion would have failed (tise_dir held 5
        # files, not 0) against the pre-fix binary.
        assert not tise_dir.exists() or not any(tise_dir.iterdir())


# ─── run_analysis_stage against the real analysis.py script ────────────────


class TestRunAnalysisStageRealSubprocess:
    """Proves controller.py's real-subprocess use of the Controller<->Analysis
    contract (docs/SDD.md Sec 7.2.3), mirroring TestRunTiseSolverRealSubprocess
    above but for the Analysis stage. controller.run_analysis_stage's `script=`
    parameter is never overridden below, so it resolves
    controller.DEFAULT_ANALYSIS_SCRIPT -- the real analysis.py sitting next to
    controller.py -- meaning every subprocess.run call in this class launches
    the genuine analysis.py CLI, exactly as controller.run() would in
    production. No mocking anywhere, matching this module's own convention.
    """

    def test_success_with_real_tise_output_and_missing_tdse_dir(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        # Run the real tise_solver first (exactly as
        # TestRunTiseSolverRealSubprocess does) to produce a genuine
        # data/tise/ directory -- Analysis's required input (Sec 7.2.2).
        tise_dir = tmp_path / "data" / "tise"
        run_tise_solver(str(tmp_config), tise_dir, binary=tise_solver_binary)

        # tdse_dir is deliberately a path that was never created. Sec 5.4.4:
        # an absent --tdse-dir is analysis.py's own tolerated case (run_tdse:
        # false is expected until Phase 8) -- it must NOT raise, and this
        # proves that tolerance survives the real subprocess boundary, not
        # just a direct in-process call to analysis.run().
        tdse_dir = tmp_path / "data" / "tdse"
        run_analysis_stage(str(tmp_config), tise_dir, tdse_dir)

    def test_failure_on_unpopulated_tise_dir_raises_solver_stage_error(self, tmp_config: Path, tmp_path: Path):
        # tise_solver_binary is deliberately NOT requested here: tise_dir
        # must genuinely never have been populated by a TISE run, so this
        # exercises a real, unmocked failure path -- signal a fully-mocked
        # test suite structurally cannot provide (per prior review feedback
        # on this phase). analysis.py's own read_tise_output() raises
        # TiseOutputError the moment it can't read eigenvalues.dat (Sec
        # 7.2.2: the four core files are non-optional), main() catches that
        # AnalysisError and exits 1, and controller._run_stage re-raises the
        # resulting CalledProcessError as a clean SolverStageError.
        tise_dir = tmp_path / "data" / "tise"  # never created; never populated
        tdse_dir = tmp_path / "data" / "tdse"

        with pytest.raises(SolverStageError) as excinfo:
            run_analysis_stage(str(tmp_config), tise_dir, tdse_dir)

        # "Analysis" is _run_stage's own stage_name for this call (see
        # run_analysis_stage's body) -- confirms Controller correctly
        # surfaced the real analysis.py subprocess's exit-1 failure as a
        # SolverStageError, not a raw/uncaught exception.
        assert "Analysis" in str(excinfo.value)

    def test_end_to_end_run_with_tise_and_analysis_enabled(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        # Build a run-analysis-enabled config variant, following the exact
        # load-YAML/mutate-a-copy/write-to-a-NEW-file pattern
        # test_continuum_enabled_missing_n_energies_raises_and_leaves_no_partial_files
        # above already uses inline -- tmp_config's own file on disk is never
        # opened for writing here, only read.
        with open(tmp_config) as f:
            cfg = yaml.safe_load(f)
        # run_tise: true and run_tdse: false are already the real config.yaml's
        # own defaults (tmp_config only overrides run.output_dir) -- only
        # run_analysis needs flipping to true here.
        cfg["run"]["run_analysis"] = True
        analysis_config = tmp_path / "config_run_analysis_enabled.yaml"
        with open(analysis_config, "w") as f:
            yaml.safe_dump(cfg, f)

        # Full real pipeline through controller.run() itself -- the actual
        # CLI entry point, not directly-called Python functions: TISE runs
        # for real against the real binary and writes real data/tise/
        # output, then Analysis runs for real (via DEFAULT_ANALYSIS_SCRIPT)
        # and reads that real TISE output back. Must not raise.
        run(str(analysis_config))

        tise_dir = tmp_path / "data" / "tise"
        assert (tise_dir / "eigenvalues.dat").is_file()


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
