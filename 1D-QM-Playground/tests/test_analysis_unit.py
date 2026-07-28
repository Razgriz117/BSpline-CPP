"""Unit tests for analysis.py (docs/SDD.md Sec 7.2.2, Sec 10.2 Phase 2).

Pure unit tests only: every fixture here is a hand-written .dat file written
directly into pytest's tmp_path, shaped to match the file formats documented
in docs/SDD.md Sec 6.3 (Persistent Storage Format) and analysis.py's own
reader docstrings. No real subprocess invocation, and no real tise_solver
binary is built or run anywhere in this file -- real-subprocess integration
tests against actual TISE solver output are a separate, later task's job, in
tests/test_analysis_integration.py.

Organized into one test class per analysis.py function, in the same order
those functions appear in analysis.py. read_eigenvectors/read_hamiltonian/
read_overlap all share identical underlying parsing logic (_read_data_rows,
_check_row_width) -- see TestReadEigenvectors's docstring for how the
error-path testing economy is handled across those three classes, mirroring
test_controller_unit.py's own handling of shared-implementation testing (its
TestRunStage tests _run_stage's error paths exhaustively once, rather than
re-testing them through every caller).
"""

from __future__ import annotations

from pathlib import Path

import pytest

from analysis import (
    AnalysisError,
    ContinuumPoint,
    EigenvalueRow,
    PhaseShiftRow,
    TiseData,
    TiseOutputError,
    read_continuum_states,
    read_eigenvalues,
    read_eigenvectors,
    read_hamiltonian,
    read_overlap,
    read_phase_shifts,
    read_tise_output,
)


# ─── Shared helpers ─────────────────────────────────────────────────────────


def _write_dat(path: Path, *lines: str) -> Path:
    """Write `path` as a .dat file: each positional arg is one raw line,
    already formatted exactly as the caller wants -- the '#' comment marker,
    field spacing, and any blank lines are all under the caller's control.
    Lines are joined with a newline plus a trailing newline at the end."""
    path.write_text("\n".join(lines) + "\n")
    return path


def _write_required_tise_files(tise_dir: Path) -> None:
    """Write minimal-but-valid eigenvalues/eigenvectors/hamiltonian/overlap
    fixtures directly into `tise_dir` (which must already exist) -- the four
    files read_tise_output requires regardless of tise.continuum.enabled.
    See TestReadTiseOutput."""
    _write_dat(tise_dir / "eigenvalues.dat", "# eigenvalues.dat: index, E_n", "0  -1.0", "1  -0.25")
    _write_dat(
        tise_dir / "eigenvectors.dat",
        "# eigenvectors.dat: columns are c_n coefficient vectors",
        "0.1  0.2",
        "0.3  0.4",
    )
    _write_dat(tise_dir / "hamiltonian.dat", "# hamiltonian.dat: H matrix (banded)", "1.0  0.0", "0.0  1.0")
    _write_dat(tise_dir / "overlap.dat", "# overlap.dat: S matrix (banded)", "1.0  0.0", "0.0  1.0")


# ─── read_eigenvalues ───────────────────────────────────────────────────────


class TestReadEigenvalues:
    def test_happy_path_parses_index_and_energy_with_field_name_access(self, tmp_path):
        path = _write_dat(
            tmp_path / "eigenvalues.dat",
            "# eigenvalues.dat: index, E_n",
            "0  -13.6",
            "1  -3.4",
            "2  -1.51",
        )

        result = read_eigenvalues(path)

        assert result == [
            EigenvalueRow(0, -13.6),
            EigenvalueRow(1, -3.4),
            EigenvalueRow(2, -1.51),
        ]
        assert result[0].index == 0
        assert result[0].energy == -13.6
        assert result[2].index == 2
        assert result[2].energy == -1.51

    def test_missing_file_raises_tise_output_error(self, tmp_path):
        path = tmp_path / "eigenvalues.dat"  # never written
        with pytest.raises(TiseOutputError) as excinfo:
            read_eigenvalues(path)
        # TiseOutputError is a subclass of AnalysisError -- confirmed here
        # once rather than in every other test, since it's a static fact
        # about the exception hierarchy, not something that varies by call.
        assert isinstance(excinfo.value, AnalysisError)

    def test_ragged_row_width_raises(self, tmp_path):
        # Row 1 has 2 fields (establishing the file's width), row 2 has 3 --
        # an internal inconsistency caught by _read_data_rows itself, distinct
        # from _check_row_width's fixed-expected-width check exercised below
        # by TestReadPhaseShifts's wrong-column-count test.
        path = _write_dat(
            tmp_path / "eigenvalues.dat",
            "# eigenvalues.dat: index, E_n",
            "0  -13.6",
            "1  -3.4  99.0",
        )
        with pytest.raises(TiseOutputError):
            read_eigenvalues(path)

    def test_non_numeric_field_raises(self, tmp_path):
        path = _write_dat(tmp_path / "eigenvalues.dat", "# eigenvalues.dat: index, E_n", "0  abc")
        with pytest.raises(TiseOutputError):
            read_eigenvalues(path)

    def test_nan_field_raises(self, tmp_path):
        # float("nan") parses without raising in plain Python -- this is the
        # explicit math.isfinite() hardening added in the fix cycle.
        path = _write_dat(tmp_path / "eigenvalues.dat", "# eigenvalues.dat: index, E_n", "1  nan")
        with pytest.raises(TiseOutputError):
            read_eigenvalues(path)

    def test_infinity_field_raises(self, tmp_path):
        # Likewise float("inf") parses fine in plain Python; must be rejected.
        path = _write_dat(tmp_path / "eigenvalues.dat", "# eigenvalues.dat: index, E_n", "1  inf")
        with pytest.raises(TiseOutputError):
            read_eigenvalues(path)

    def test_fractional_index_raises(self, tmp_path):
        # A non-whole-number index (e.g. a corrupted file) must raise rather
        # than being silently truncated by int().
        path = _write_dat(tmp_path / "eigenvalues.dat", "# eigenvalues.dat: index, E_n", "2.9  -1.0")
        with pytest.raises(TiseOutputError):
            read_eigenvalues(path)

    def test_whole_number_float_index_parses_to_int(self, tmp_path):
        # "2.0" is a whole number written with a decimal point -- a legal
        # on-disk representation that must still parse to the Python int 2,
        # not be left as (or rejected for being) a float.
        path = _write_dat(tmp_path / "eigenvalues.dat", "# eigenvalues.dat: index, E_n", "2.0  -1.0")

        result = read_eigenvalues(path)

        assert result == [EigenvalueRow(2, -1.0)]
        assert isinstance(result[0].index, int)


# ─── read_eigenvectors ──────────────────────────────────────────────────────


class TestReadEigenvectors:
    """read_eigenvectors, read_hamiltonian, and read_overlap all delegate
    straight to the same _read_data_rows() with no extra per-function logic
    on top (unlike read_eigenvalues, which layers its own index-integer
    check on the result) -- so the missing-file/invalid-encoding/
    ragged-row/non-numeric/NaN/Infinity error paths are byte-for-byte
    identical across all three, and are tested exhaustively HERE, once,
    rather than being re-tested per file. This mirrors
    test_controller_unit.py's TestRunStage, which tests _run_stage's error
    handling exhaustively in one place rather than re-testing it through
    every caller that goes through it.
    TestReadHamiltonian/TestReadOverlap below still each get their own
    happy-path test, since eigenvectors.dat/hamiltonian.dat/overlap.dat are
    semantically distinct files.
    """

    def test_happy_path_parses_rectangular_grid(self, tmp_path):
        path = _write_dat(
            tmp_path / "eigenvectors.dat",
            "# eigenvectors.dat: columns are c_n coefficient vectors",
            "0.1  0.2",
            "0.3  0.4",
            "0.5  0.6",
        )

        assert read_eigenvectors(path) == [
            [0.1, 0.2],
            [0.3, 0.4],
            [0.5, 0.6],
        ]

    def test_missing_file_raises(self, tmp_path):
        with pytest.raises(TiseOutputError):
            read_eigenvectors(tmp_path / "eigenvectors.dat")

    def test_ragged_row_raises(self, tmp_path):
        path = _write_dat(
            tmp_path / "eigenvectors.dat",
            "# eigenvectors.dat",
            "0.1  0.2",
            "0.3  0.4  0.5",
        )
        with pytest.raises(TiseOutputError):
            read_eigenvectors(path)

    def test_non_numeric_field_raises(self, tmp_path):
        path = _write_dat(tmp_path / "eigenvectors.dat", "# eigenvectors.dat", "0.1  abc")
        with pytest.raises(TiseOutputError):
            read_eigenvectors(path)

    def test_nan_field_raises(self, tmp_path):
        path = _write_dat(tmp_path / "eigenvectors.dat", "# eigenvectors.dat", "0.1  nan")
        with pytest.raises(TiseOutputError):
            read_eigenvectors(path)

    def test_infinity_field_raises(self, tmp_path):
        path = _write_dat(tmp_path / "eigenvectors.dat", "# eigenvectors.dat", "0.1  inf")
        with pytest.raises(TiseOutputError):
            read_eigenvectors(path)

    def test_invalid_utf8_bytes_raises_tise_output_error(self, tmp_path):
        # 0xFF is never a valid UTF-8 leading byte, so path.read_text(encoding=
        # "utf-8") is guaranteed to raise UnicodeDecodeError on this fixture --
        # exercising _read_data_rows's `except UnicodeDecodeError` branch
        # specifically (a binary/corrupted file that can't be read as text),
        # distinct from the `except OSError` branch already covered by
        # test_missing_file_raises above. Must surface as TiseOutputError, not
        # let the raw UnicodeDecodeError propagate.
        path = tmp_path / "eigenvectors.dat"
        path.write_bytes(b"\xff\xfe\x00\x01\x02\x03")
        with pytest.raises(TiseOutputError):
            read_eigenvectors(path)

    def test_blank_lines_and_comments_skipped_anywhere_not_just_at_top(self, tmp_path):
        # _read_data_rows's docstring is explicit that blank/comment lines
        # are skipped WHEREVER they occur, not just as a leading header --
        # a regression to "only skip line 1" would silently misparse this.
        path = _write_dat(
            tmp_path / "eigenvectors.dat",
            "# eigenvectors.dat: columns are c_n coefficient vectors",
            "0.1  0.2",
            "",
            "# a stray mid-file comment",
            "0.3  0.4",
        )

        assert read_eigenvectors(path) == [[0.1, 0.2], [0.3, 0.4]]


# ─── read_hamiltonian ───────────────────────────────────────────────────────


class TestReadHamiltonian:
    """See TestReadEigenvectors's docstring: read_hamiltonian shares that
    same _read_data_rows parsing/error-handling verbatim, so only a
    happy-path test is needed here to confirm hamiltonian.dat's own file is
    wired up correctly."""

    def test_happy_path_parses_rectangular_grid(self, tmp_path):
        path = _write_dat(
            tmp_path / "hamiltonian.dat",
            "# hamiltonian.dat: H matrix (banded)",
            "1.0  -0.5  0.0",
            "-0.5  1.0  -0.5",
        )

        assert read_hamiltonian(path) == [
            [1.0, -0.5, 0.0],
            [-0.5, 1.0, -0.5],
        ]


# ─── read_overlap ───────────────────────────────────────────────────────────


class TestReadOverlap:
    """See TestReadEigenvectors's docstring: read_overlap shares that same
    _read_data_rows parsing/error-handling verbatim, so only a happy-path
    test is needed here to confirm overlap.dat's own file is wired up
    correctly."""

    def test_happy_path_parses_rectangular_grid(self, tmp_path):
        path = _write_dat(
            tmp_path / "overlap.dat",
            "# overlap.dat: S matrix (banded)",
            "1.0  0.2",
            "0.2  1.0",
        )

        assert read_overlap(path) == [
            [1.0, 0.2],
            [0.2, 1.0],
        ]


# ─── read_phase_shifts ──────────────────────────────────────────────────────


class TestReadPhaseShifts:
    def test_happy_path_parses_energy_delta_ddelta_with_field_name_access(self, tmp_path):
        path = _write_dat(
            tmp_path / "phase_shifts.dat",
            "# phase_shifts.dat: eps_i, delta(eps_i), d(delta)/dE",
            "0.1  0.05  1.2",
            "0.2  0.09  1.1",
        )

        result = read_phase_shifts(path)

        assert result == [
            PhaseShiftRow(0.1, 0.05, 1.2),
            PhaseShiftRow(0.2, 0.09, 1.1),
        ]
        assert result[0].energy == 0.1
        assert result[0].delta == 0.05
        assert result[0].ddelta_dE == 1.2

    def test_missing_file_returns_empty_list_not_raise(self, tmp_path):
        # tise.continuum.enabled: false -- phase_shifts.dat simply isn't
        # written, and that is a normal, successful outcome, not an error.
        path = tmp_path / "phase_shifts.dat"  # never written
        assert read_phase_shifts(path) == []

    def test_wrong_column_count_raises(self, tmp_path):
        # All rows agree with each other (so _read_data_rows sees no
        # raggedness) but the shared width is 2, not the 3 phase_shifts.dat
        # requires -- this is what exercises _check_row_width's own raise,
        # distinct from _read_data_rows's internal-consistency raise.
        path = _write_dat(
            tmp_path / "phase_shifts.dat",
            "# phase_shifts.dat: eps_i, delta(eps_i), d(delta)/dE",
            "0.1  0.05",
        )
        with pytest.raises(TiseOutputError):
            read_phase_shifts(path)

    def test_non_numeric_field_raises(self, tmp_path):
        path = _write_dat(
            tmp_path / "phase_shifts.dat",
            "# phase_shifts.dat: eps_i, delta(eps_i), d(delta)/dE",
            "0.1  abc  1.2",
        )
        with pytest.raises(TiseOutputError):
            read_phase_shifts(path)

    def test_nan_field_raises(self, tmp_path):
        path = _write_dat(
            tmp_path / "phase_shifts.dat",
            "# phase_shifts.dat: eps_i, delta(eps_i), d(delta)/dE",
            "0.1  nan  1.2",
        )
        with pytest.raises(TiseOutputError):
            read_phase_shifts(path)

    def test_infinity_field_raises(self, tmp_path):
        path = _write_dat(
            tmp_path / "phase_shifts.dat",
            "# phase_shifts.dat: eps_i, delta(eps_i), d(delta)/dE",
            "0.1  inf  1.2",
        )
        with pytest.raises(TiseOutputError):
            read_phase_shifts(path)


# ─── read_continuum_states ──────────────────────────────────────────────────


class TestReadContinuumStates:
    def test_happy_path_multiple_files_ascending_order_and_ignores_unrelated_files(self, tmp_path):
        # An unrelated required-file also present, exactly like a real
        # tise_dir would have -- proves the continuum_state_NNN.dat regex
        # match is selective, not "every file in the directory."
        _write_dat(tmp_path / "eigenvalues.dat", "# eigenvalues.dat: index, E_n", "0  -1.0")
        _write_dat(
            tmp_path / "continuum_state_001.dat",
            "# continuum_state_001.dat: x, psi",
            "0.0  0.0",
            "1.0  0.5",
        )
        _write_dat(
            tmp_path / "continuum_state_002.dat",
            "# continuum_state_002.dat: x, psi",
            "0.0  0.0",
            "1.0  0.7",
        )

        result = read_continuum_states(tmp_path)

        assert result == [
            (1, [ContinuumPoint(0.0, 0.0), ContinuumPoint(1.0, 0.5)]),
            (2, [ContinuumPoint(0.0, 0.0), ContinuumPoint(1.0, 0.7)]),
        ]
        # Field-name access, not just positional/tuple equality.
        assert result[0][1][1].x == 1.0
        assert result[0][1][1].psi == 0.5

    def test_no_matching_files_in_existing_dir_returns_empty_list(self, tmp_path):
        # tise.continuum.enabled: false -- the directory exists (as every
        # tmp_path does) and holds other real files, but none matching
        # continuum_state_NNN.dat. Must return [], not raise.
        (tmp_path / "eigenvalues.dat").write_text("# eigenvalues.dat: index, E_n\n0  -1.0\n")
        assert read_continuum_states(tmp_path) == []

    def test_nonexistent_dir_returns_empty_list(self, tmp_path):
        missing_dir = tmp_path / "does_not_exist"
        assert read_continuum_states(missing_dir) == []

    def test_malformed_continuum_file_raises_and_names_the_file(self, tmp_path):
        # Unlike the other five readers (which each pass a fixed literal
        # description string, e.g. "eigenvalues.dat"), read_continuum_states
        # passes each file's OWN path.name as the description -- so the
        # raised error should name the specific offending file.
        _write_dat(
            tmp_path / "continuum_state_001.dat",
            "# continuum_state_001.dat: x, psi",
            "0.0  0.0  0.0",  # 3 fields; continuum files require 2
        )
        with pytest.raises(TiseOutputError, match="continuum_state_001.dat"):
            read_continuum_states(tmp_path)

    def test_sorts_by_parsed_integer_not_filename_string(self, tmp_path):
        # Deliberately unpadded/mixed-width filenames written in a shuffled
        # order. Lexicographic ("string") order would compare
        # "continuum_state_1.dat" < "continuum_state_10.dat" <
        # "continuum_state_2.dat" (since '.' sorts before '0', and '1'
        # sorts before '2'), giving [1, 10, 2]. The correct numeric-integer
        # sort is [1, 2, 10]. A regression to a naive string sort key would
        # make this assertion fail.
        _write_dat(tmp_path / "continuum_state_2.dat", "# continuum_state_2.dat: x, psi", "0.0  0.2")
        _write_dat(tmp_path / "continuum_state_10.dat", "# continuum_state_10.dat: x, psi", "0.0  0.10")
        _write_dat(tmp_path / "continuum_state_1.dat", "# continuum_state_1.dat: x, psi", "0.0  0.1")

        result = read_continuum_states(tmp_path)

        indices = [index for index, _points in result]
        assert indices == [1, 2, 10]

    def test_directory_matching_the_filename_pattern_is_ignored(self, tmp_path):
        # entry.is_file() guards against a same-named subdirectory being
        # mistaken for a continuum-state file (Path.iterdir() yields
        # directories too). Without that guard this would instead raise
        # (read_text() on a directory is an OSError), so this genuinely
        # exercises the guard rather than trivially passing either way.
        (tmp_path / "continuum_state_099.dat").mkdir()
        assert read_continuum_states(tmp_path) == []


# ─── read_tise_output (aggregator) ──────────────────────────────────────────


class TestReadTiseOutput:
    def test_full_happy_path_all_six_files_present_continuum_enabled(self, tmp_path):
        _write_required_tise_files(tmp_path)
        _write_dat(
            tmp_path / "phase_shifts.dat",
            "# phase_shifts.dat: eps_i, delta(eps_i), d(delta)/dE",
            "0.1  0.05  1.2",
            "0.2  0.09  1.1",
        )
        _write_dat(
            tmp_path / "continuum_state_001.dat",
            "# continuum_state_001.dat: x, psi",
            "0.0  0.0",
            "1.0  0.5",
        )
        _write_dat(
            tmp_path / "continuum_state_002.dat",
            "# continuum_state_002.dat: x, psi",
            "0.0  0.0",
            "1.0  0.7",
        )

        result = read_tise_output(tmp_path)

        assert result == TiseData(
            eigenvalues=[EigenvalueRow(0, -1.0), EigenvalueRow(1, -0.25)],
            eigenvectors=[[0.1, 0.2], [0.3, 0.4]],
            hamiltonian=[[1.0, 0.0], [0.0, 1.0]],
            overlap=[[1.0, 0.0], [0.0, 1.0]],
            phase_shifts=[PhaseShiftRow(0.1, 0.05, 1.2), PhaseShiftRow(0.2, 0.09, 1.1)],
            continuum_states=[
                (1, [ContinuumPoint(0.0, 0.0), ContinuumPoint(1.0, 0.5)]),
                (2, [ContinuumPoint(0.0, 0.0), ContinuumPoint(1.0, 0.7)]),
            ],
        )

    def test_continuum_absent_case_returns_empty_phase_shifts_and_continuum_states(self, tmp_path):
        # Only the four required files present -- tise.continuum.enabled:
        # false. Must NOT raise; phase_shifts/continuum_states degrade to
        # [] while the four required fields are still fully populated.
        _write_required_tise_files(tmp_path)

        result = read_tise_output(tmp_path)

        assert result.phase_shifts == []
        assert result.continuum_states == []
        assert result.eigenvalues == [EigenvalueRow(0, -1.0), EigenvalueRow(1, -0.25)]
        assert result.eigenvectors == [[0.1, 0.2], [0.3, 0.4]]
        assert result.hamiltonian == [[1.0, 0.0], [0.0, 1.0]]
        assert result.overlap == [[1.0, 0.0], [0.0, 1.0]]

    def test_missing_required_file_raises(self, tmp_path):
        # hamiltonian.dat is deliberately never written; eigenvalues.dat and
        # eigenvectors.dat ARE written first so the failure is unambiguously
        # at the hamiltonian step (source order: eigenvalues, eigenvectors,
        # hamiltonian, overlap, phase_shifts, continuum), not coincidentally
        # from some earlier missing file too.
        _write_dat(tmp_path / "eigenvalues.dat", "# eigenvalues.dat: index, E_n", "0  -1.0")
        _write_dat(
            tmp_path / "eigenvectors.dat",
            "# eigenvectors.dat: columns are c_n coefficient vectors",
            "0.1  0.2",
        )

        with pytest.raises(TiseOutputError, match="hamiltonian"):
            read_tise_output(tmp_path)

    def test_nonexistent_tise_dir_raises(self, tmp_path):
        # No special-casing exists for a nonexistent tise_dir: read_eigenvalues
        # runs first and naturally raises because its own file can't be found.
        # match="eigenvalues" confirms it's specifically that first reader
        # failing, not some other unrelated error.
        missing_dir = tmp_path / "does_not_exist"
        with pytest.raises(TiseOutputError, match="eigenvalues"):
            read_tise_output(missing_dir)
