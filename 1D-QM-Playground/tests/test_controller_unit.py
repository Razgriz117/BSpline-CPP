"""Unit tests for controller.py (docs/SDD.md Sec 7.2.1, Sec 10.2 Phase 1).

Pure unit tests only: no real subprocess invocation of tise_solver, and no
real filesystem writes that stand in for the solver's own output. Wherever
controller.py talks to a subprocess (`_run_stage`, `run_tise_solver`),
`subprocess.run` is mocked via `unittest.mock.patch`; where run_tise_solver's
own `Path.mkdir` call needs to fail, `Path.mkdir` itself is patched instead
of touching a real binary. Real-subprocess integration tests against the
actual tise_solver binary are Task 5's job, in a separate test module in
this same tests/ directory.

Organized into one test class per controller.py function, in roughly the
same order those functions appear in controller.py.
"""

from __future__ import annotations

import io
import json
import subprocess
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from controller import (
    ConfigValidationError,
    ControllerError,
    Interval,
    SolverStageError,
    _run_stage,
    load_config,
    main,
    parse_interval,
    parse_potential_piece,
    print_warnings,
    read_warnings,
    run,
    run_tise_solver,
    validate_config,
    validate_potential_tiling,
)


# ─── Shared helpers ─────────────────────────────────────────────────────────


def _piece(domain: str, function: str = "0") -> dict:
    """A piece dict shaped like parse_potential_piece's return value (i.e.
    what validate_potential_tiling actually consumes), for testing tiling
    logic directly without going through the string-parsing layer."""
    return {"domain": domain, "function": function}


def _minimal_cfg(run_tise=False, run_tdse=False, run_analysis=False, output_dir="out") -> dict:
    """A config dict that passes validate_config: one potential piece that
    exactly tiles the bspline domain [0, 10]."""
    return {
        "run": {
            "run_tise": run_tise,
            "run_tdse": run_tdse,
            "run_analysis": run_analysis,
            "output_dir": output_dir,
        },
        "bspline": {"domain": [0.0, 10.0]},
        "potential": ["{'domain': '[0, 10]', 'function': '0'}"],
    }


def _write_yaml(path: Path, data: dict) -> Path:
    path.write_text(yaml.safe_dump(data))
    return path


# ─── parse_interval ─────────────────────────────────────────────────────────


class TestParseInterval:
    def test_closed_interval(self):
        assert parse_interval("[0,100]") == Interval(lo=0.0, lo_inclusive=True, hi=100.0, hi_inclusive=True)

    def test_open_interval(self):
        assert parse_interval("(0,100)") == Interval(lo=0.0, lo_inclusive=False, hi=100.0, hi_inclusive=False)

    def test_half_open_lo_exclusive_hi_inclusive(self):
        assert parse_interval("(0, 100]") == Interval(lo=0.0, lo_inclusive=False, hi=100.0, hi_inclusive=True)

    def test_half_open_lo_inclusive_hi_exclusive(self):
        assert parse_interval("[0, 100)") == Interval(lo=0.0, lo_inclusive=True, hi=100.0, hi_inclusive=False)

    def test_inf_upper_bound(self):
        result = parse_interval("[0, inf)")
        assert result == Interval(lo=0.0, lo_inclusive=True, hi=float("inf"), hi_inclusive=False)

    def test_inf_lower_bound(self):
        result = parse_interval("(-inf, 0]")
        assert result == Interval(lo=float("-inf"), lo_inclusive=False, hi=0.0, hi_inclusive=True)

    def test_both_bounds_infinite(self):
        result = parse_interval("(-inf, inf)")
        assert result == Interval(lo=float("-inf"), lo_inclusive=False, hi=float("inf"), hi_inclusive=False)

    def test_whitespace_around_bounds_tolerated(self):
        assert parse_interval(" [ 0 , 100 ] ") == Interval(lo=0.0, lo_inclusive=True, hi=100.0, hi_inclusive=True)

    @pytest.mark.parametrize(
        "bad",
        [
            "",              # too short, no brackets at all
            "0,100",         # missing both brackets
            "[0 100]",       # missing the comma
            "[0,50,100]",    # too many commas
            "[a,100]",       # non-numeric bound
            "[100,0]",       # lo > hi
            "[0,100",        # missing closing bracket
            "0,100]",        # missing opening bracket
        ],
    )
    def test_malformed_strings_raise(self, bad):
        with pytest.raises(ConfigValidationError):
            parse_interval(bad)

    def test_non_string_input_raises(self):
        with pytest.raises(ConfigValidationError):
            parse_interval(123)


# ─── parse_potential_piece ──────────────────────────────────────────────────


class TestParsePotentialPiece:
    def test_valid_single_quoted_dict_parses(self):
        piece = "{'domain': '(0, 100]', 'function': '-1/x + 1/x^2'}"
        assert parse_potential_piece(piece) == {"domain": "(0, 100]", "function": "-1/x + 1/x^2"}

    def test_trailing_comma_tolerated(self):
        # A trailing comma is valid Python dict-literal syntax.
        piece = "{'domain': '[0,5)', 'function': '0',}"
        assert parse_potential_piece(piece) == {"domain": "[0,5)", "function": "0"}

    def test_json_bareword_true_is_rejected(self):
        # Valid JSON, but NOT valid Python literal syntax: `true` parses as
        # a bare name reference (not the Python literal `True`), which
        # ast.literal_eval refuses to evaluate. Confirms this function uses
        # ast.literal_eval, not json.loads (which would accept this fine).
        piece = '{"domain": "[0,5)", "function": "0", "enabled": true}'
        with pytest.raises(ConfigValidationError):
            parse_potential_piece(piece)

    def test_malformed_syntax_raises(self):
        with pytest.raises(ConfigValidationError):
            parse_potential_piece("{'domain': '[0,5)', 'function': }")

    def test_parses_but_not_a_dict_raises(self):
        with pytest.raises(ConfigValidationError):
            parse_potential_piece("[1, 2, 3]")

    def test_missing_domain_key_raises(self):
        with pytest.raises(ConfigValidationError):
            parse_potential_piece("{'function': '0'}")

    def test_missing_function_key_raises(self):
        with pytest.raises(ConfigValidationError):
            parse_potential_piece("{'domain': '[0,5)'}")

    def test_non_string_input_raises(self):
        with pytest.raises(ConfigValidationError):
            parse_potential_piece(123)


# ─── validate_potential_tiling ──────────────────────────────────────────────


class TestValidatePotentialTiling:
    """The trickiest logic in controller.py -- see validate_potential_tiling's
    own docstring for the exact tiling rules. Each test is named for the
    scenario it exercises."""

    def test_single_piece_whole_domain_real_config_case(self):
        # The REAL config.yaml's exact case: potential is a single
        # "(0, 100]" piece (open at x=0 because -1/x is singular there)
        # against bspline.domain = [0.0, 100.0] -- open at the domain's OWN
        # left edge. The first piece is deliberately not required to be
        # inclusive at domain_lo. Load-bearing regression test: must PASS.
        validate_potential_tiling([_piece("(0, 100]")], (0.0, 100.0))

    def test_multi_piece_correct_tiling(self):
        # Every internal boundary has exactly one inclusive side.
        pieces = [_piece("[0,5)"), _piece("[5,6]"), _piece("(6,10]")]
        validate_potential_tiling(pieces, (0.0, 10.0))

    def test_genuine_gap_raises(self):
        # [0,5) then [6,10]: nothing covers (5, 6).
        pieces = [_piece("[0,5)"), _piece("[6,10]")]
        with pytest.raises(ConfigValidationError):
            validate_potential_tiling(pieces, (0.0, 10.0))

    def test_genuine_overlap_both_inclusive_raises(self):
        # [0,5] then [5,10]: both inclusive at the shared point x=5 --
        # ambiguous which piece's function value applies there.
        pieces = [_piece("[0,5]"), _piece("[5,10]")]
        with pytest.raises(ConfigValidationError):
            validate_potential_tiling(pieces, (0.0, 10.0))

    def test_touching_boundary_both_exclusive_tolerated(self):
        # [0,5) then (5,10]: both EXCLUSIVE at the shared point x=5. A
        # measure-zero gap, not numerically meaningful for a continuous
        # potential -- must be TOLERATED (pass), not rejected.
        pieces = [_piece("[0,5)"), _piece("(5,10]")]
        validate_potential_tiling(pieces, (0.0, 10.0))

    def test_pieces_out_of_order_still_validate(self):
        # Deliberately not given in sorted order; validate_potential_tiling
        # sorts internally by lower bound before checking adjacency.
        pieces = [_piece("[5,6]"), _piece("[0,5)"), _piece("(6,10]")]
        validate_potential_tiling(pieces, (0.0, 10.0))

    def test_domain_not_covered_at_start_raises(self):
        pieces = [_piece("[1,10]")]
        with pytest.raises(ConfigValidationError):
            validate_potential_tiling(pieces, (0.0, 10.0))

    def test_domain_not_covered_at_end_raises(self):
        pieces = [_piece("[0,9]")]
        with pytest.raises(ConfigValidationError):
            validate_potential_tiling(pieces, (0.0, 10.0))

    def test_infinite_domain_single_piece(self):
        pieces = [_piece("(-inf, inf)")]
        validate_potential_tiling(pieces, (float("-inf"), float("inf")))

    def test_empty_pieces_list_raises(self):
        with pytest.raises(ConfigValidationError):
            validate_potential_tiling([], (0.0, 10.0))


# ─── load_config ────────────────────────────────────────────────────────────


class TestLoadConfig:
    def test_valid_config_loads(self, tmp_path):
        cfg = _minimal_cfg()
        path = _write_yaml(tmp_path / "config.yaml", cfg)
        assert load_config(str(path)) == cfg

    def test_missing_file_raises_with_clear_message(self, tmp_path):
        missing = tmp_path / "does_not_exist.yaml"
        with pytest.raises(ConfigValidationError, match="not found"):
            load_config(str(missing))

    def test_unparseable_yaml_raises_with_clear_message(self, tmp_path):
        path = tmp_path / "bad.yaml"
        path.write_text("key: [1, 2\n")  # unclosed flow sequence -> YAMLError
        with pytest.raises(ConfigValidationError, match="failed to parse"):
            load_config(str(path))

    def test_non_mapping_yaml_raises(self, tmp_path):
        path = tmp_path / "list.yaml"
        path.write_text("- 1\n- 2\n")  # parses fine, but to a list, not a dict
        with pytest.raises(ConfigValidationError, match="did not parse to a mapping"):
            load_config(str(path))

    def test_unreadable_path_raises_via_generic_oserror_branch(self, tmp_path):
        # A directory is not FileNotFoundError -- open() raises
        # IsADirectoryError, a *different* OSError subclass, exercising the
        # generic `except OSError` branch (distinct from the
        # FileNotFoundError-specific one above).
        with pytest.raises(ConfigValidationError, match="could not read config file"):
            load_config(str(tmp_path))


# ─── validate_config ────────────────────────────────────────────────────────


class TestValidateConfig:
    def test_valid_config_passes(self):
        validate_config(_minimal_cfg())  # must not raise

    def test_missing_potential_key_raises(self):
        cfg = _minimal_cfg()
        del cfg["potential"]
        with pytest.raises(ConfigValidationError, match="potential"):
            validate_config(cfg)

    def test_potential_not_a_list_raises(self):
        cfg = _minimal_cfg()
        cfg["potential"] = "not a list"
        with pytest.raises(ConfigValidationError):
            validate_config(cfg)

    def test_potential_empty_list_raises(self):
        cfg = _minimal_cfg()
        cfg["potential"] = []
        with pytest.raises(ConfigValidationError):
            validate_config(cfg)

    def test_missing_bspline_domain_raises(self):
        cfg = _minimal_cfg()
        del cfg["bspline"]["domain"]
        with pytest.raises(ConfigValidationError, match="bspline.domain"):
            validate_config(cfg)

    def test_domain_not_a_lo_hi_pair_raises(self):
        cfg = _minimal_cfg()
        cfg["bspline"]["domain"] = "nope"
        with pytest.raises(ConfigValidationError):
            validate_config(cfg)

    def test_bad_tiling_propagates_as_config_validation_error(self):
        cfg = _minimal_cfg()
        # Only covers [0, 5] of the [0.0, 10.0] bspline domain.
        cfg["potential"] = ["{'domain': '[0, 5]', 'function': '0'}"]
        with pytest.raises(ConfigValidationError):
            validate_config(cfg)


# ─── _run_stage (subprocess.run mocked) ─────────────────────────────────────


class TestRunStage:
    """subprocess.run is ALWAYS mocked in this class -- no real process is
    ever spawned by these tests."""

    def test_success_returns_completed_process(self):
        with patch("controller.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(
                args=["tise_solver"], returncode=0, stdout="ok", stderr=""
            )
            result = _run_stage("TISE solver", ["tise_solver", "--config", "c.yaml"])

        assert result.returncode == 0
        mock_run.assert_called_once_with(
            ["tise_solver", "--config", "c.yaml"], check=True, capture_output=True, text=True
        )

    def test_called_process_error_reraised_as_solver_stage_error(self):
        with patch("controller.subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(
                returncode=42, cmd=["tise_solver"], output="", stderr="physics blew up: NaN encountered"
            )
            with pytest.raises(SolverStageError) as excinfo:
                _run_stage("TISE solver", ["tise_solver"])

        message = str(excinfo.value)
        assert "TISE solver" in message
        assert "physics blew up: NaN encountered" in message

    def test_file_not_found_reraised_as_solver_stage_error(self):
        with patch("controller.subprocess.run") as mock_run:
            mock_run.side_effect = FileNotFoundError("No such file or directory: 'tise_solver'")
            with pytest.raises(SolverStageError) as excinfo:
                _run_stage("TISE solver", ["tise_solver"])

        message = str(excinfo.value)
        assert "TISE solver" in message
        assert "No such file or directory" in message

    def test_generic_os_error_reraised_as_solver_stage_error(self):
        # Hardening added after an earlier code review: any bare OSError
        # (e.g. a permission error hitting the binary) must ALSO become a
        # clean SolverStageError, not propagate as a raw OSError traceback.
        with patch("controller.subprocess.run") as mock_run:
            mock_run.side_effect = OSError("Permission denied")
            with pytest.raises(SolverStageError) as excinfo:
                _run_stage("TISE solver", ["tise_solver"])

        message = str(excinfo.value)
        assert "TISE solver" in message
        assert "Permission denied" in message


# ─── run_tise_solver (subprocess.run mocked) ────────────────────────────────


class TestRunTiseSolver:
    def test_success_creates_output_dir_and_invokes_stage(self, tmp_path):
        out_dir = tmp_path / "tise"
        with patch("controller.subprocess.run") as mock_run:
            mock_run.return_value = subprocess.CompletedProcess(args=["x"], returncode=0, stdout="", stderr="")
            run_tise_solver("c.yaml", out_dir, binary=Path("/fake/tise_solver"))

        assert out_dir.is_dir()
        mock_run.assert_called_once_with(
            ["/fake/tise_solver", "--config", "c.yaml", "--output-dir", str(out_dir)],
            check=True,
            capture_output=True,
            text=True,
        )

    def test_mkdir_oserror_becomes_clean_solver_stage_error(self):
        # run_tise_solver's OWN mkdir call, not _run_stage's subprocess
        # handling. Must become a clean SolverStageError, not a propagated
        # OSError traceback. subprocess.run is also mocked here so that,
        # even if the mkdir patch somehow failed to raise, this test still
        # could not fall through to a real subprocess invocation.
        with patch.object(Path, "mkdir", side_effect=OSError("Permission denied")), \
                patch("controller.subprocess.run") as mock_run:
            with pytest.raises(SolverStageError) as excinfo:
                run_tise_solver("c.yaml", Path("/fake/out"))

        assert "Permission denied" in str(excinfo.value)
        mock_run.assert_not_called()


# ─── read_warnings ───────────────────────────────────────────────────────────


class TestReadWarnings:
    """Real file I/O against tmp_path is fine here: warnings.json is plain
    JSON, not the tise_solver binary itself."""

    def test_valid_warnings_json_parses(self, tmp_path):
        warnings = [{"category": "info", "message": "hello"}]
        (tmp_path / "warnings.json").write_text(json.dumps(warnings))
        assert read_warnings(tmp_path) == warnings

    def test_missing_file_returns_empty_list_and_notes_stderr(self, tmp_path, capsys):
        result = read_warnings(tmp_path)  # no warnings.json written at all
        assert result == []
        captured = capsys.readouterr()
        assert "warnings.json" in captured.err
        assert captured.out == ""

    def test_malformed_json_returns_empty_list_and_notes_stderr(self, tmp_path, capsys):
        (tmp_path / "warnings.json").write_text("{not valid json")
        result = read_warnings(tmp_path)
        assert result == []
        captured = capsys.readouterr()
        assert "warnings.json" in captured.err

    def test_non_list_json_returns_empty_list_and_notes_stderr(self, tmp_path, capsys):
        (tmp_path / "warnings.json").write_text(json.dumps({"not": "a list"}))
        result = read_warnings(tmp_path)
        assert result == []
        captured = capsys.readouterr()
        assert "warnings.json" in captured.err


# ─── print_warnings ─────────────────────────────────────────────────────────


class TestPrintWarnings:
    def test_omitted_stream_is_captured_by_capsys(self, capsys):
        # Regression test for a real bug: the original signature defaulted
        # to `stream=sys.stderr`, a value bound ONCE at function-definition
        # (module-import) time. That stale reference keeps getting written
        # to even after capsys reassigns sys.stderr for a given test, so
        # capsys.readouterr().err would come back empty -- the test would
        # have FAILED against that old default. It was fixed to
        # `stream=None`, resolved inside the function body on every call,
        # which capsys does observe. Calling print_warnings with NO
        # explicit stream= and asserting on capsys.readouterr().err is
        # therefore the only way this test actually exercises the fix.
        warnings = [
            {"category": "info", "message": "hi"},
            {"category": "warn", "message": "careful"},
        ]
        print_warnings(warnings)
        captured = capsys.readouterr()
        assert "[info] hi" in captured.err
        assert "[warn] careful" in captured.err
        assert captured.out == ""

    def test_explicit_stream_is_honored(self):
        buf = io.StringIO()
        print_warnings([{"category": "info", "message": "hi"}], stream=buf)
        assert "[info] hi" in buf.getvalue()

    def test_missing_keys_fall_back_to_defaults(self, capsys):
        print_warnings([{}])
        captured = capsys.readouterr()
        assert "[unknown]" in captured.err

    def test_empty_list_prints_nothing(self, capsys):
        print_warnings([])
        captured = capsys.readouterr()
        assert captured.err == ""


# ─── run() orchestration ────────────────────────────────────────────────────


class TestRunOrchestrationGuards:
    def test_run_tdse_true_raises_controller_error(self, tmp_path):
        cfg = _minimal_cfg(run_tdse=True)
        path = _write_yaml(tmp_path / "config.yaml", cfg)
        with pytest.raises(ControllerError, match="run_tdse"):
            run(str(path))

    def test_run_analysis_true_raises_controller_error(self, tmp_path):
        cfg = _minimal_cfg(run_analysis=True)
        path = _write_yaml(tmp_path / "config.yaml", cfg)
        with pytest.raises(ControllerError, match="run_analysis"):
            run(str(path))

    def test_all_stages_disabled_is_a_silent_no_op(self, tmp_path):
        cfg = _minimal_cfg()  # run_tise/run_tdse/run_analysis all False
        path = _write_yaml(tmp_path / "config.yaml", cfg)
        assert run(str(path)) is None

    def test_run_tise_dispatches_to_solver_and_warnings(self, tmp_path):
        # run_tise_solver itself is mocked here -- no real subprocess call --
        # so this exercises run()'s wiring (which functions get called, with
        # what arguments) without touching the real tise_solver binary.
        output_dir = tmp_path / "outdir"
        cfg = _minimal_cfg(run_tise=True, output_dir=str(output_dir))
        path = _write_yaml(tmp_path / "config.yaml", cfg)

        with patch("controller.run_tise_solver") as mock_solver, \
                patch("controller.read_warnings") as mock_read, \
                patch("controller.print_warnings") as mock_print:
            mock_read.return_value = [{"category": "info", "message": "hi"}]
            run(str(path))

        mock_solver.assert_called_once_with(str(path), output_dir / "tise")
        mock_read.assert_called_once_with(output_dir / "tise")
        mock_print.assert_called_once_with([{"category": "info", "message": "hi"}])


# ─── main() CLI wrapper ──────────────────────────────────────────────────────


class TestMain:
    """run() itself is mocked here -- these tests only exercise argparse
    wiring and the ControllerError -> exit-code-1 translation, not any real
    orchestration or subprocess behavior."""

    def test_success_returns_zero(self):
        with patch("controller.run") as mock_run:
            rc = main(["--config", "my_config.yaml"])
        assert rc == 0
        mock_run.assert_called_once_with("my_config.yaml")

    def test_controller_error_returns_one_and_prints_to_stderr(self, capsys):
        with patch("controller.run") as mock_run:
            mock_run.side_effect = ControllerError("boom")
            rc = main(["--config", "my_config.yaml"])
        assert rc == 1
        captured = capsys.readouterr()
        assert "boom" in captured.err

    def test_default_config_argument_is_config_yaml(self):
        with patch("controller.run") as mock_run:
            main([])
        mock_run.assert_called_once_with("config.yaml")
