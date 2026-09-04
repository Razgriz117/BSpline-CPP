"""analysis.py -- TISE-output reading layer and CLI for the 1D QM Playground.

Reads the files the TISE solver writes under data/tise/ and parses them
into plain Python containers, per docs/SDD.md Sec 7.2.2 (TISE Solver to
Analysis) and Sec 6.3 (Persistent Storage Format). Also implements the
Controller-to-Analysis subprocess CLI (Sec 7.2.3):

    analysis.py --config <config.yaml> --tise-dir <data/tise/> --tdse-dir <data/tdse/>

Phase 3 scope (docs/SDD.md Sec 10.2, Sec 7.2.3): main()/run()/load_config()
implement that CLI's subprocess/exit-code contract -- load config.yaml
independently of controller.py, read the required data/tise/ output (the
six reader functions plus read_tise_output(), from Phase 2, unchanged
here), and tolerate --tdse-dir not existing (Sec 5.4.4: run_tdse: false is
expected until Phase 8). Deliberately still NOT implemented here (later
phases):
  - No reading of data/tdse/ *content* -- only --tdse-dir's existence is
    checked; actually reading snapshot_NNNNN.dat/observables.dat (Sec
    7.2.5) is Phase 6's job.
  - No reading of data/tise/warnings.json -- that sidecar is a
    Controller-facing artifact (Sec 8), not one of Sec 7.2.2's Analysis
    inputs.
  - No binary-format support for hamiltonian.dat/overlap.dat -- Sec 6.3
    leaves "plain text or binary" undecided; only plain text exists today.
  - No plotting.
  - No computation of any REQ-F-060 quantity (bound-state populations,
    expectation values, spectral distributions, etc.) -- those require
    TDSE output (Sec 7.2.5) too and are Phase 8's job.
  - No output artifact of any kind: no plot files, no derived-data files,
    no data/analysis/ directory, no placeholder stdout payload. run()
    computes and writes nothing -- deliberately deferred per ADR-0005
    (docs/adr/0005-defer-analysis-output-artifact-format.md), not an
    oversight.
"""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import NamedTuple
import matplotlib.pyplot as plt

import yaml

# Matches continuum_state_NNN.dat filenames, capturing the numeric index.
# Deliberately \d+ rather than \d{3}: Sec 6.3 documents zero-padded 3-digit
# indices as the normal case, but read_continuum_states() must sort by the
# PARSED INTEGER, not the filename string (see its docstring) -- so parsing
# is kept lenient about width rather than assuming exactly 3 digits.
_CONTINUUM_STATE_RE = re.compile(r"^continuum_state_(\d+)\.dat$")

# Same leniency-about-width reasoning as _CONTINUUM_STATE_RE above.
_EIGENSTATE_RE = re.compile(r"^eigenstate_(\d+)\.dat$")


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


class ConfigError(AnalysisError):
    """config.yaml could not be loaded for Analysis (missing, unreadable,
    unparseable, or not a mapping).

    Deliberately narrower than controller.ConfigValidationError: this
    module's load_config() performs none of controller.py's solver-specific
    schema checks (potential-piece tiling, bspline.domain, run.run_tise/
    output_dir) -- Phase 3 doesn't read or act on the `analysis`/
    `visualization` blocks' contents yet (their REQ-F-060 quantities all
    require TDSE output, unavailable until Phase 8 -- see this module's
    top docstring). That's a separate point from ADR-0005
    (docs/adr/0005-defer-analysis-output-artifact-format.md), which defers
    only Phase 3's own output-artifact format/location, not input
    validation.
    """


class EigenvalueRow(NamedTuple):
    """One eigenvalues.dat row: a basis-state index and its energy E_n."""

    index: int
    energy: float


class PhaseShiftRow(NamedTuple):
    """One phase_shifts.dat row: eps_i, delta(eps_i), and d(delta)/dE."""

    energy: float
    delta: float
    ddelta_dE: float


class ContinuumPoint(NamedTuple):
    """One row of a continuum_state_NNN.dat file: a grid point x, psi_eps(x)."""

    x: float
    psi: float


class EigenstatePoint(NamedTuple):
    """One row of an eigenstate_NNN.dat file: a grid point x, phi_n(x)."""

    x: float
    phi: float


@dataclass(frozen=True)
class TiseData:
    """Bundled result of reading every file under one data/tise/ directory.

    phase_shifts and continuum_states are [] when tise.continuum.enabled
    was false for the run that produced tise_dir (see
    read_phase_shifts/read_continuum_states) -- an empty list here is a
    normal, successful result, not a sentinel for "not read yet".
    eigenstates, by contrast, is unconditional -- tise_solver writes an
    eigenstate_NNN.dat for every computed bound state regardless of
    tise.continuum.enabled (see read_eigenstates); [] there means the run
    predates that output existing, or the directory is otherwise empty.

    Note: frozen=True stops a field from being REASSIGNED (e.g.
    `some_tise_data.eigenvalues = [...]` raises FrozenInstanceError), but
    does not deep-freeze its contents -- every field here is a plain,
    still-mutable Python list (only the per-row tuples are the immutable,
    self-documenting NamedTuple types below). Nothing in this module
    mutates them after construction, but a caller technically could.
    """

    eigenvalues: list[EigenvalueRow]
    eigenvectors: list[list[float]]
    hamiltonian: list[list[float]]
    overlap: list[list[float]]
    phase_shifts: list[PhaseShiftRow]
    continuum_states: list[tuple[int, list[ContinuumPoint]]]
    eigenstates: list[tuple[int, list[EigenstatePoint]]]


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
    has a data row containing a non-numeric field, or has a field that
    parses but isn't finite (NaN/Infinity, e.g. Python's own float("nan")/
    float("inf") both parse without raising, so this needs an explicit
    math.isfinite() check of its own -- a solver that failed to converge
    could plausibly write one of these once Phase 4 lands real numerics,
    and it should be caught here rather than flow silently downstream). A
    single except clause covers "missing" and "otherwise unreadable"
    together (both are OSError subclasses, e.g. FileNotFoundError/
    PermissionError/IsADirectoryError) since str(e) already carries the
    OS's own reason.
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

        row: list[float] = []
        for field in fields:
            try:
                value = float(field)
            except ValueError as e:
                raise TiseOutputError(
                    f"malformed {description} at {path}: line {lineno} has a non-numeric field: {e}"
                ) from e
            if not math.isfinite(value):
                raise TiseOutputError(
                    f"malformed {description} at {path}: line {lineno} has a non-finite field "
                    f"(NaN/Infinity not permitted): {field!r}"
                )
            row.append(value)
        rows.append(row)

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


def read_eigenvalues(path: Path) -> list[EigenvalueRow]:
    """Read eigenvalues.dat: 2 columns per row, (index, E_n).

    The index field must be a whole number (its file representation may
    still be e.g. "2" or "2.0" -- either parses fine). A fractional index
    (e.g. "2.9") is treated as malformed and raises, rather than being
    silently truncated by int(), consistent with this reader's "corruption
    is not tolerated" stance on its other fields.

    Raises TiseOutputError if `path` is missing, unreadable, or malformed.
    """
    rows = _read_data_rows(path, "eigenvalues.dat")
    _check_row_width(rows, 2, path, "eigenvalues.dat")

    result: list[EigenvalueRow] = []
    for row in rows:
        index_value = row[0]
        if not index_value.is_integer():
            raise TiseOutputError(
                f"malformed eigenvalues.dat at {path}: index field {index_value!r} is not a whole number"
            )
        result.append(EigenvalueRow(int(index_value), row[1]))
    return result


def read_eigenvectors(path: Path) -> list[list[float]]:
    """Read eigenvectors.dat: a rectangular grid of floats, no index column.

    Each row is one basis-coefficient row; each column is one eigenstate's
    coefficient. Raises TiseOutputError if `path` is missing, unreadable,
    or malformed (rows of inconsistent width, a non-numeric field, or a
    non-finite field).
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


def read_phase_shifts(path: Path) -> list[PhaseShiftRow]:
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
    return [PhaseShiftRow(row[0], row[1], row[2]) for row in rows]


def read_continuum_states(tise_dir: Path) -> list[tuple[int, list[ContinuumPoint]]]:
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

    result: list[tuple[int, list[ContinuumPoint]]] = []
    for index, path in matches:
        rows = _read_data_rows(path, path.name)
        _check_row_width(rows, 2, path, path.name)
        result.append((index, [ContinuumPoint(row[0], row[1]) for row in rows]))

    return result


def read_eigenstates(tise_dir: Path) -> list[tuple[int, list[EigenstatePoint]]]:
    """Read every eigenstate_NNN.dat file under `tise_dir`.

    Each matching file holds 2 columns per row, (x, phi). Returns a list of
    (index, [(x, phi), ...]) pairs, one per file, sorted by ASCENDING
    NUMERIC index parsed from the filename's NNN -- same int()-not-string
    sort as read_continuum_states() above.

    Returns [] (does NOT raise) if `tise_dir` doesn't exist, or exists but
    contains no matching files. Unlike continuum_state_NNN.dat, tise_solver
    writes these unconditionally (one per bound state, regardless of
    tise.continuum.enabled), so [] here means either an older tise_dir that
    predates this output, or an otherwise-empty directory -- not a
    config-driven "disabled" case. If a matching file IS present but
    malformed, this raises TiseOutputError -- absence is tolerated,
    corruption is not.
    """
    if not tise_dir.is_dir():
        return []

    matches: list[tuple[int, Path]] = []
    for entry in tise_dir.iterdir():
        if not entry.is_file():
            continue
        m = _EIGENSTATE_RE.match(entry.name)
        if m:
            matches.append((int(m.group(1)), entry))

    matches.sort(key=lambda pair: pair[0])

    result: list[tuple[int, list[EigenstatePoint]]] = []
    for index, path in matches:
        rows = _read_data_rows(path, path.name)
        _check_row_width(rows, 2, path, path.name)
        result.append((index, [EigenstatePoint(row[0], row[1]) for row in rows]))

    return result


# ─── Aggregator ──────────────────────────────────────────────────────────────


def read_tise_output(tise_dir: Path) -> TiseData:
    """Read every file under one data/tise/ directory into one TiseData.

    Calls all seven reader functions above. No special-casing is needed for
    a nonexistent `tise_dir`: read_eigenvalues() runs first and naturally
    raises TiseOutputError (its file can't be found under a nonexistent
    directory) before any of the other six readers run.
    """
    eigenvalues = read_eigenvalues(tise_dir / "eigenvalues.dat")
    eigenvectors = read_eigenvectors(tise_dir / "eigenvectors.dat")
    hamiltonian = read_hamiltonian(tise_dir / "hamiltonian.dat")
    overlap = read_overlap(tise_dir / "overlap.dat")
    phase_shifts = read_phase_shifts(tise_dir / "phase_shifts.dat")
    continuum_states = read_continuum_states(tise_dir)
    eigenstates = read_eigenstates(tise_dir)

    return TiseData(
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        hamiltonian=hamiltonian,
        overlap=overlap,
        phase_shifts=phase_shifts,
        continuum_states=continuum_states,
        eigenstates=eigenstates,
    )


# ─── Controller-to-Analysis CLI (Sec 7.2.3) ─────────────────────────────────


def load_config(config_path: str) -> dict:
    """Load and parse config.yaml, independently of controller.py.

    Raises ConfigError (never a raw exception) if the file doesn't exist,
    can't be read, fails to parse as YAML, or doesn't parse to a mapping.
    """
    try:
        with open(config_path, "r") as f:
            cfg = yaml.safe_load(f)
    except FileNotFoundError as e:
        raise ConfigError(f"config file not found: {config_path}") from e
    except OSError as e:
        raise ConfigError(f"could not read config file {config_path}: {e}") from e
    except UnicodeDecodeError as e:
        raise ConfigError(f"config file {config_path} is not valid text: {e}") from e
    except yaml.YAMLError as e:
        raise ConfigError(f"failed to parse config file {config_path}: {e}") from e

    if not isinstance(cfg, dict):
        raise ConfigError(
            f"config file {config_path} did not parse to a mapping (got {type(cfg).__name__})"
        )
    return cfg


def plot_tise(tise_data: TiseData, tise_dir: str) -> None:
    # plot continuum states
    continuum = tise_data.continuum_states
    for idx, continuum_state in continuum:
        plt.xlim(0, continuum_state[-1].x)
        plt.xlabel("r")
        plt.ylabel("psi")
        # reformat continuum state data for matplotlib
        xs = [p.x for p in continuum_state]
        psis = [p.psi for p in continuum_state]
        plt.plot(xs, psis)
        plt.title(fr"Continuum state {idx}", pad=20)
        plt.savefig(f"{tise_dir}/continuum_{idx}.png")
        plt.close()


def plot_eigenstates(tise_data: TiseData, tise_dir: str) -> None:
    """Plot |phi_n(x)|^2 for each bound eigenstate (docs/SDD.md Sec 6.1's
    visualization.eigenstates toggle) -- squared, unlike plot_tise's raw
    psi, per that same spec entry."""
    for idx, eigenstate in tise_data.eigenstates:
        xs = [p.x for p in eigenstate]
        densities = [p.phi**2 for p in eigenstate]
        plt.xlim(eigenstate[0].x, eigenstate[-1].x)
        plt.xlabel("x")
        plt.ylabel(r"$|\phi_n(x)|^2$")
        plt.plot(xs, densities)
        plt.title(fr"Eigenstate {idx}", pad=20)
        plt.savefig(f"{tise_dir}/eigenstate_{idx}.png")
        plt.close()


def plot_phase_shifts(tise_data: TiseData, tise_dir: str) -> None:
    """Plot delta(E) and d(delta)/dE vs. E (docs/SDD.md Sec 6.1's
    visualization.phase_shifts toggle -- TISE-only, not TDSE-gated, unlike
    most of its visualization.* neighbors). One PNG, two stacked subplots;
    writes nothing when tise.continuum.enabled was false for this run
    (phase_shifts == [], same tolerate-absence convention plot_tise's own
    continuum-state loop already follows)."""
    if not tise_data.phase_shifts:
        return

    energies = [row.energy for row in tise_data.phase_shifts]
    deltas = [row.delta for row in tise_data.phase_shifts]
    ddelta_dEs = [row.ddelta_dE for row in tise_data.phase_shifts]

    fig, (ax_delta, ax_ddelta) = plt.subplots(2, 1, sharex=True)
    ax_delta.plot(energies, deltas)
    ax_delta.set_ylabel(r"$\delta(E)$")
    ax_ddelta.plot(energies, ddelta_dEs)
    ax_ddelta.set_xlabel("E")
    ax_ddelta.set_ylabel(r"$d\delta/dE$")
    fig.suptitle("Phase shift")
    fig.savefig(f"{tise_dir}/phase_shifts.png")
    plt.close(fig)


def run(config_path: str, tise_dir: str, tdse_dir: str) -> None:
    """Run Analysis's Phase 3 scope (Sec 7.2.3, Sec 10.2 Phase 3).

    Loads `config_path`, reads the required data/tise/ output at
    `tise_dir` (raises TiseOutputError if missing/malformed -- Sec 7.2.2's
    required files are non-optional), and tolerates `tdse_dir` not
    existing (Sec 5.4.4: run_tdse: false is expected until Phase 8) with a
    stderr note rather than raising. Computes and writes nothing (ADR-0005).

    Raises AnalysisError (ConfigError or TiseOutputError) on failure;
    never a raw exception.
    """
    cfg = load_config(config_path)
    tise_output = read_tise_output(Path(tise_dir))

    plot_tise(tise_output, tise_dir)
    if cfg.get("visualization", {}).get("eigenstates", False):
        plot_eigenstates(tise_output, tise_dir)
    if cfg.get("visualization", {}).get("phase_shifts", False):
        plot_phase_shifts(tise_output, tise_dir)

    if not Path(tdse_dir).is_dir():
        print(
            f"analysis.py: info: --tdse-dir {tdse_dir} not found or not a directory; "
            f"proceeding without TDSE-derived analysis (run_tdse: false is expected until Phase 8)",
            file=sys.stderr,
        )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="analysis.py",
        description="Analysis for the 1D QM Playground: reads TISE/TDSE output per config.yaml (Sec 7.2.3).",
    )
    parser.add_argument("--config", required=True, help="path to config.yaml")
    parser.add_argument("--tise-dir", required=True, help="path to data/tise/ (required)")
    parser.add_argument("--tdse-dir", required=True, help="path to data/tdse/ (tolerated absent)")
    args = parser.parse_args(argv)

    try:
        run(args.config, args.tise_dir, args.tdse_dir)
    except AnalysisError as e:
        print(f"analysis.py: {e}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
