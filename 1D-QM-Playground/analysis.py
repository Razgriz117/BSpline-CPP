"""analysis.py -- TISE-output reading layer for the 1D QM Playground.

Reads the files the TISE solver writes under data/tise/ and parses them
into plain Python containers, per docs/SDD.md Sec 7.2.2 (TISE Solver to
Analysis) and Sec 6.3 (Persistent Storage Format).

Phase 2 scope only (docs/SDD.md Sec 10.2): this module implements ONLY the
data-reading side of the TISE-to-Analysis contract -- six reader functions
plus read_tise_output(), the aggregator that bundles their results into one
TiseData. Deliberately NOT implemented here (all later phases):
  - No CLI, no main(), no argparse, no `if __name__ == "__main__":` block
    -- the Controller-to-Analysis subprocess invocation (Sec 7.2.3,
    `analysis.py --config <path> --tise-dir <path> --tdse-dir <path>`) is
    a future phase.
  - No reading of data/tise/warnings.json -- that sidecar is a
    Controller-facing artifact (Sec 8), not one of Sec 7.2.2's Analysis
    inputs.
  - No binary-format support for hamiltonian.dat/overlap.dat -- Sec 6.3
    leaves "plain text or binary" undecided; only plain text exists today.
  - No plotting.
  - No computation of any REQ-F-060 quantity (bound-state populations,
    expectation values, spectral distributions, etc.) -- those require
    TDSE output (Sec 7.2.5) too and are a much later phase.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

# Matches continuum_state_NNN.dat filenames, capturing the numeric index.
# Deliberately \d+ rather than \d{3}: Sec 6.3 documents zero-padded 3-digit
# indices as the normal case, but read_continuum_states() must sort by the
# PARSED INTEGER, not the filename string (see its docstring) -- so parsing
# is kept lenient about width rather than assuming exactly 3 digits.
_CONTINUUM_STATE_RE = re.compile(r"^continuum_state_(\d+)\.dat$")


class AnalysisError(Exception):
    """Base class for all errors raised by the Analysis module."""


class TiseOutputError(AnalysisError):
    """A TISE output file is missing, unreadable, or malformed.

    Raised unconditionally by the four always-present-file readers
    (read_eigenvalues, read_eigenvectors, read_hamiltonian, read_overlap).
    Raised by read_phase_shifts/read_continuum_states only when their
    file(s) are PRESENT but malformed -- both tolerate their file(s) being
    simply absent (see those functions' docstrings) by returning [] instead.
    """


@dataclass(frozen=True)
class TiseData:
    """Bundled result of reading every file under one data/tise/ directory.

    phase_shifts and continuum_states are [] when tise.continuum.enabled
    was false for the run that produced tise_dir (see
    read_phase_shifts/read_continuum_states) -- an empty list here is a
    normal, successful result, not a sentinel for "not read yet".
    """

    eigenvalues: list[tuple[int, float]]
    eigenvectors: list[list[float]]
    hamiltonian: list[list[float]]
    overlap: list[list[float]]
    phase_shifts: list[tuple[float, float, float]]
    continuum_states: list[tuple[int, list[tuple[float, float]]]]


# ─── Generic .dat row parsing ───────────────────────────────────────────────


def _read_data_rows(path: Path, description: str) -> list[list[float]]:
    """Read `path` and parse its data rows into a list of float lists.

    Every `.dat` file (docs/SDD.md Sec 6.3) has one '#'-prefixed comment
    header line, then data rows of whitespace-separated numeric fields --
    fields may be separated by any run of whitespace, not necessarily
    exactly one space. Blank lines and '#'-prefixed lines are skipped
    wherever they occur (not just as the first line); every other line is
    treated as one data row.

    Raises TiseOutputError if `path` doesn't exist, can't be read as text,
    has a data row with a different field count than the first data row,
    or has a data row containing a non-numeric field. A single except
    clause covers "missing" and "otherwise unreadable" together (both are
    OSError subclasses, e.g. FileNotFoundError/PermissionError/
    IsADirectoryError) since str(e) already carries the OS's own reason.
    """
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as e:
        raise TiseOutputError(f"could not read {description} at {path}: {e}") from e
    except UnicodeDecodeError as e:
        raise TiseOutputError(f"{description} at {path} is not valid text: {e}") from e

    rows: list[list[float]] = []
    ncols: int | None = None
    for lineno, line in enumerate(text.splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        fields = stripped.split()
        if ncols is None:
            ncols = len(fields)
        elif len(fields) != ncols:
            raise TiseOutputError(
                f"malformed {description} at {path}: line {lineno} has {len(fields)} "
                f"field(s), expected {ncols} (from the first data row)"
            )

        try:
            rows.append([float(field) for field in fields])
        except ValueError as e:
            raise TiseOutputError(
                f"malformed {description} at {path}: line {lineno} has a non-numeric field: {e}"
            ) from e

    return rows


def _check_row_width(rows: list[list[float]], width: int, path: Path, description: str) -> None:
    """Raise TiseOutputError if non-empty `rows` doesn't have `width` fields
    per row.

    `rows` from _read_data_rows() is already self-consistent (every row has
    the same field count as the first), so checking rows[0] suffices.
    Vacuously fine on an empty file (zero data rows) -- there is no row to
    be the wrong width.
    """
    if rows and len(rows[0]) != width:
        raise TiseOutputError(
            f"malformed {description} at {path}: rows have {len(rows[0])} field(s), expected {width}"
        )


# ─── Required-file readers (always present after a successful TISE run) ────


def read_eigenvalues(path: Path) -> list[tuple[int, float]]:
    """Read eigenvalues.dat: 2 columns per row, (index, E_n).

    Raises TiseOutputError if `path` is missing, unreadable, or malformed.
    """
    rows = _read_data_rows(path, "eigenvalues.dat")
    _check_row_width(rows, 2, path, "eigenvalues.dat")
    return [(int(row[0]), row[1]) for row in rows]


def read_eigenvectors(path: Path) -> list[list[float]]:
    """Read eigenvectors.dat: a rectangular grid of floats, no index column.

    Each row is one basis-coefficient row; each column is one eigenstate's
    coefficient. Raises TiseOutputError if `path` is missing, unreadable,
    or malformed (rows of inconsistent width, or a non-numeric field).
    """
    return _read_data_rows(path, "eigenvectors.dat")


def read_hamiltonian(path: Path) -> list[list[float]]:
    """Read hamiltonian.dat: a rectangular grid of floats, no index column.

    Banded-storage layout is not yet defined anywhere in docs/SDD.md (Sec
    6.3) -- this parses whatever rectangular grid is in the plain-text
    file, with no attempt to interpret bandwidth/layout semantics. Raises
    TiseOutputError if `path` is missing, unreadable, or malformed.
    """
    return _read_data_rows(path, "hamiltonian.dat")


def read_overlap(path: Path) -> list[list[float]]:
    """Read overlap.dat: a rectangular grid of floats, no index column.

    See read_hamiltonian()'s docstring for the banded-storage caveat,
    which applies identically here. Raises TiseOutputError if `path` is
    missing, unreadable, or malformed.
    """
    return _read_data_rows(path, "overlap.dat")


# ─── Continuum-tolerant readers (absent when tise.continuum.enabled: false) ─


def read_phase_shifts(path: Path) -> list[tuple[float, float, float]]:
    """Read phase_shifts.dat: 3 columns per row, (eps_i, delta_i, ddelta_dE_i).

    Only written when tise.continuum.enabled: true. Returns [] (does NOT
    raise) if `path` doesn't exist -- this is what satisfies docs/SDD.md
    Sec 7.2.2's Contract ("Analysis must tolerate tise.continuum.enabled:
    false"). If `path` DOES exist but is malformed, this raises
    TiseOutputError exactly like the required-file readers above --
    absence is tolerated, corruption is not.
    """
    if not path.exists():
        return []

    rows = _read_data_rows(path, "phase_shifts.dat")
    _check_row_width(rows, 3, path, "phase_shifts.dat")
    return [(row[0], row[1], row[2]) for row in rows]


def read_continuum_states(tise_dir: Path) -> list[tuple[int, list[tuple[float, float]]]]:
    """Read every continuum_state_NNN.dat file under `tise_dir`.

    Each matching file holds 2 columns per row, (x, psi). Returns a list of
    (index, [(x, psi), ...]) pairs, one per file, sorted by ASCENDING
    NUMERIC index parsed from the filename's NNN -- an explicit int() sort,
    not string/lexicographic order (zero-padding makes the two agree for
    up to 999 states, but the sort key here is the parsed integer either
    way) -- and never from any config value; this module has no dependency
    on config.yaml/tise.continuum.n_energies.

    Returns [] (does NOT raise) if `tise_dir` doesn't exist, or exists but
    contains no matching files -- this is what satisfies docs/SDD.md Sec
    7.2.2's Contract ("Analysis must tolerate tise.continuum.enabled:
    false"), mirroring read_phase_shifts() above. If a matching file IS
    present but malformed, this raises TiseOutputError -- absence is
    tolerated, corruption is not.
    """
    if not tise_dir.is_dir():
        return []

    matches: list[tuple[int, Path]] = []
    for entry in tise_dir.iterdir():
        if not entry.is_file():
            continue
        m = _CONTINUUM_STATE_RE.match(entry.name)
        if m:
            matches.append((int(m.group(1)), entry))

    matches.sort(key=lambda pair: pair[0])

    result: list[tuple[int, list[tuple[float, float]]]] = []
    for index, path in matches:
        rows = _read_data_rows(path, path.name)
        _check_row_width(rows, 2, path, path.name)
        result.append((index, [(row[0], row[1]) for row in rows]))

    return result


# ─── Aggregator ──────────────────────────────────────────────────────────────


def read_tise_output(tise_dir: Path) -> TiseData:
    """Read every file under one data/tise/ directory into one TiseData.

    Calls all six reader functions above. No special-casing is needed for
    a nonexistent `tise_dir`: read_eigenvalues() runs first and naturally
    raises TiseOutputError (its file can't be found under a nonexistent
    directory) before any of the other five readers run.
    """
    eigenvalues = read_eigenvalues(tise_dir / "eigenvalues.dat")
    eigenvectors = read_eigenvectors(tise_dir / "eigenvectors.dat")
    hamiltonian = read_hamiltonian(tise_dir / "hamiltonian.dat")
    overlap = read_overlap(tise_dir / "overlap.dat")
    phase_shifts = read_phase_shifts(tise_dir / "phase_shifts.dat")
    continuum_states = read_continuum_states(tise_dir)

    return TiseData(
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        hamiltonian=hamiltonian,
        overlap=overlap,
        phase_shifts=phase_shifts,
        continuum_states=continuum_states,
    )
