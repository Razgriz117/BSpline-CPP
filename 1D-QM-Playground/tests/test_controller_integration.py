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
controller.py, is proven separately by test_analysis_integration.py's
TestAnalysisCliRealSubprocess class, which invokes analysis.py directly as
a raw subprocess with no controller.py involvement at all; this module
only proves controller.py's correct *usage* of each contract, not
analysis.py's own CLI behavior under malformed/incomplete argv.)

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
        # config.yaml's own default has had tise.continuum.enabled: true
        # since commit 99b7683 -- force it off here so this test still
        # exercises the no-continuum-output contract it's named for.
        with open(tmp_config) as f:
            cfg = yaml.safe_load(f)
        cfg["tise"]["continuum"]["enabled"] = False
        no_continuum_config = tmp_path / "config_no_continuum.yaml"
        with open(no_continuum_config, "w") as f:
            yaml.safe_dump(cfg, f)

        tise_dir = tmp_path / "data" / "tise"

        run_tise_solver(str(no_continuum_config), tise_dir, binary=tise_solver_binary)

        for name in ("eigenvalues.dat", "eigenvectors.dat", "hamiltonian.dat", "overlap.dat", "warnings.json"):
            assert (tise_dir / name).is_file(), f"expected output file missing: {name}"

        # continuum disabled above -- the continuum-only outputs must NOT be written.
        assert not (tise_dir / "phase_shifts.dat").exists()
        assert list(tise_dir.glob("continuum_state_*.dat")) == []

        # warnings.json is no longer unconditionally empty here: Part B
        # (docs/planning/tise-release-readiness-plan.md) wired
        # classifyBoundStates/checkWellContainment to run regardless of
        # tise.continuum.enabled -- the real hydrogen config's own 7 bound
        # states (3 of them only marginally contained by R=100) always
        # produce these two diagnostics now. What this test still confirms:
        # no CONTINUUM-specific warning (E_acc, pole-proximity) leaks in
        # when continuum is genuinely disabled.
        warnings = json.loads((tise_dir / "warnings.json").read_text())
        messages = [w["message"] for w in warnings]
        assert any("computed states are below E=0.0" in m for m in messages)
        assert not any("E_acc" in m or "discretization artifact" in m for m in messages)

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

    def test_overlapping_potential_pieces_raises_and_leaves_no_partial_files(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        """validateNoOverlappingPotentialPieces (Part C,
        docs/planning/tise-release-readiness-plan.md): two pieces both
        covering x=50 (both inclusive at that shared boundary) previously
        resolved silently via std::map's own domain-string sort order, with
        no error at all -- now raises before any output is written,
        matching the same "non-zero exit == no partial output" contract the
        n_energies test above guards."""
        with open(tmp_config) as f:
            cfg = yaml.safe_load(f)
        cfg["potential"] = [
            "{'domain': '[0, 50]', 'function': '0'}",
            "{'domain': '[50, 100]', 'function': '1'}",
        ]
        bad_config = tmp_path / "config_overlapping_potential.yaml"
        with open(bad_config, "w") as f:
            yaml.safe_dump(cfg, f)

        tise_dir = tmp_path / "data" / "tise"

        with pytest.raises(SolverStageError) as excinfo:
            run_tise_solver(str(bad_config), tise_dir, binary=tise_solver_binary)

        assert "TISE solver" in str(excinfo.value)
        assert not tise_dir.exists() or not any(tise_dir.iterdir())

    def test_non_unity_mass_raises_and_leaves_no_partial_files(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        """physics.mass/physics.hbar guard-rail (Part C): both fields are
        documented in config.yaml but were never consumed anywhere --
        fillBandedMatrices hardcodes mass=1 internally, so a user setting
        physics.mass to anything else got silent wrong physics with no
        error. Guard-rail only (not full mass/hbar generalization, which
        would touch k=sqrt(2E)/computeEAcc/the kinetic-energy term
        throughout): an honest config error instead."""
        with open(tmp_config) as f:
            cfg = yaml.safe_load(f)
        cfg["physics"]["mass"] = 2.0
        bad_config = tmp_path / "config_non_unity_mass.yaml"
        with open(bad_config, "w") as f:
            yaml.safe_dump(cfg, f)

        tise_dir = tmp_path / "data" / "tise"

        with pytest.raises(SolverStageError) as excinfo:
            run_tise_solver(str(bad_config), tise_dir, binary=tise_solver_binary)

        assert "TISE solver" in str(excinfo.value)
        assert not tise_dir.exists() or not any(tise_dir.iterdir())


# ─── run_analysis_stage against the real analysis.py script ────────────────


class TestRunAnalysisStageRealSubprocess:
    """Proves controller.py's real-subprocess use of the Controller<->Analysis
    contract (docs/SDD.md Sec 7.2.3), mirroring TestRunTiseSolverRealSubprocess
    above but for the Analysis stage. controller.run_analysis_stage's `script=`
    parameter is never overridden below, so it resolves
    controller.DEFAULT_ANALYSIS_SCRIPT -- the real analysis.py sitting next to
    controller.py -- meaning every run_analysis_stage/run() call below
    launches the genuine analysis.py CLI, exactly as controller.run() would
    in production. (The first two tests also subprocess-invoke the real
    tise_solver binary first, to produce genuine fixture input -- see each
    test's own comments.) No mocking anywhere, matching this module's own
    convention.
    """

    def test_success_with_real_tise_output_and_missing_tdse_dir(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path, capsys
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

        # run_analysis_stage is called IN-PROCESS here (this test never
        # shells out to controller.py itself), so capsys observes THIS
        # process's stderr. The real analysis.py subprocess prints its own
        # "--tdse-dir ... not found" info note to ITS stderr -- Analysis's
        # only channel back to a user at all, since it produces no output
        # artifact of any kind (ADR-0005) -- and run_analysis_stage must
        # relay that captured stderr onward rather than silently
        # discarding it. This is the real-subprocess proof that the relay
        # fix actually works end to end, not just against a mocked
        # CompletedProcess (see test_controller_unit.py's
        # TestRunAnalysisStage.test_success_relays_stderr_to_controller_stderr
        # for the mocked equivalent).
        captured = capsys.readouterr()
        assert "--tdse-dir" in captured.err
        assert "not found or not a directory" in captured.err

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

        # The one artifact this whole pipeline exists to produce -- plots --
        # was never actually asserted on before (only eigenvalues.dat's
        # existence was checked): the real config.yaml has both
        # tise.continuum.enabled and visualization.eigenstates on by
        # default, so a real run_analysis: true pass through this exact
        # config must produce both continuum and eigenstate PNGs, not just
        # write eigenvalues.dat and silently skip plotting.
        continuum_png = tise_dir / "continuum_1.png"
        eigenstate_png = tise_dir / "eigenstate_1.png"
        assert continuum_png.is_file() and continuum_png.stat().st_size > 0
        assert eigenstate_png.is_file() and eigenstate_png.stat().st_size > 0

    def test_run_analysis_false_produces_no_analysis_output(
        self, tmp_config: Path, tise_solver_binary: Path, tmp_path: Path
    ):
        # Explicit, dedicated real-subprocess counterpart to the test above:
        # run.run_analysis: false (the real config.yaml's own default,
        # tmp_config leaves it untouched) must produce zero Analysis-stage
        # output -- no PNGs of any kind -- not just "the test happens to
        # never check for them" as the pre-existing TestRunEndToEnd.
        # test_run_end_to_end_success below does incidentally.
        run(str(tmp_config))

        tise_dir = tmp_path / "data" / "tise"
        assert (tise_dir / "eigenvalues.dat").is_file()  # TISE stage still ran
        assert list(tise_dir.glob("*.png")) == []


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
        # config.yaml's own default has had tise.continuum.enabled: true
        # since commit 99b7683, which (with its E_max=10 default) now
        # legitimately produces continuum-specific physics warnings (E_max
        # exceeds the basis's own accuracy ceiling) -- force continuum off
        # here so this test isolates the non-continuum-specific diagnostics
        # (Part B's classifyBoundStates/checkWellContainment, which now run
        # regardless of continuum -- see TestRunTiseSolverRealSubprocess.
        # test_success_writes_expected_files_and_no_continuum_output above)
        # rather than conflating them with continuum-only warnings.
        with open(tmp_config) as f:
            cfg = yaml.safe_load(f)
        cfg["tise"]["continuum"]["enabled"] = False
        no_continuum_config = tmp_path / "config_no_continuum.yaml"
        with open(no_continuum_config, "w") as f:
            yaml.safe_dump(cfg, f)

        run(str(no_continuum_config))

        tise_dir = tmp_path / "data" / "tise"
        assert (tise_dir / "eigenvalues.dat").is_file()
        warnings = json.loads((tise_dir / "warnings.json").read_text())
        messages = [w["message"] for w in warnings]
        assert any("computed states are below E=0.0" in m for m in messages)
        assert not any("E_acc" in m or "discretization artifact" in m for m in messages)

        # read_warnings/print_warnings ran for real too: a non-empty
        # warnings list means each entry gets printed to stderr, one line
        # per entry, "[category] message".
        captured = capsys.readouterr()
        assert "[physics]" in captured.err
        assert "computed states are below E=0.0" in captured.err
