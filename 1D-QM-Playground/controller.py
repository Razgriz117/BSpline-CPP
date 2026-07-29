"""controller.py -- orchestration entry point for the 1D QM Playground.

Loads config.yaml, validates it, and invokes solver subprocesses per
docs/SDD.md Sec 5.1 (Controller) and Sec 7.2 (Inter-Component Interfaces).

This module implements the Controller <-> TISE Solver contract (Sec 7.2.1,
Phase 1) and the Controller <-> Analysis contract (Sec 7.2.3, Phase 3). TDSE
orchestration (Sec 7.2.4, Phase 5) is a future phase and is deliberately NOT
implemented here -- run() raises a clear error if a config asks for it
rather than silently no-op-ing.
"""

from __future__ import annotations

import argparse
import ast
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

import yaml

# Tolerance used when comparing potential-piece boundaries to each other and
# to the bspline domain edges (see validate_potential_tiling).
EPS = 1e-9

# TISE/build/tise_solver resolved relative to this file's OWN location (not
# the current working directory), so controller.py is runnable from any
# directory. This deliberately points at the real, current build location
# (TISE/build/tise_solver) rather than the SDD Sec 11.2 "target" unified
# ./build/tise_solver layout -- that directory-layout migration is a
# separate, out-of-scope concern.
DEFAULT_TISE_SOLVER = Path(__file__).resolve().parent / "TISE" / "build" / "tise_solver"

# analysis.py resolved the same way, for the same reason (see above).
DEFAULT_ANALYSIS_SCRIPT = Path(__file__).resolve().parent / "analysis.py"


class ControllerError(Exception):
    """Base class for all errors raised by the Controller."""


class ConfigValidationError(ControllerError):
    """config.yaml could not be loaded, or failed validation."""


class SolverStageError(ControllerError):
    """A solver subprocess stage failed: non-zero exit or missing binary."""


@dataclass(frozen=True)
class Interval:
    """A mathematical interval with independently open/closed endpoints.

    lo/hi may be +-inf. `lo_inclusive`/`hi_inclusive` record whether that
    endpoint is closed ("[" / "]", True) or open ("(" / ")", False).
    """

    lo: float
    lo_inclusive: bool
    hi: float
    hi_inclusive: bool


# ─── Config loading & validation ───────────────────────────────────────────


def load_config(config_path: str) -> dict:
    """Load and parse config.yaml.

    Raises ConfigValidationError (never a raw exception) if the file
    doesn't exist, can't be read, fails to parse as YAML, or doesn't parse
    to a mapping.
    """
    try:
        with open(config_path, "r") as f:
            cfg = yaml.safe_load(f)
    except FileNotFoundError as e:
        raise ConfigValidationError(f"config file not found: {config_path}") from e
    except OSError as e:
        raise ConfigValidationError(f"could not read config file {config_path}: {e}") from e
    except yaml.YAMLError as e:
        raise ConfigValidationError(f"failed to parse config file {config_path}: {e}") from e

    if not isinstance(cfg, dict):
        raise ConfigValidationError(
            f"config file {config_path} did not parse to a mapping (got {type(cfg).__name__})"
        )

    return cfg


def validate_config(cfg: dict) -> None:
    """Validate `cfg`.

    Checks that the 'run' section exists and is a mapping containing the
    'run_tise' and 'output_dir' keys -- the two run() indexes directly with
    raw `[...]` (as opposed to 'run_tdse'/'run_analysis', which run() reads
    via `.get()` and are genuinely optional) -- so a config missing either
    fails here with a clear ConfigValidationError instead of run() raising
    an uncaught KeyError later.

    Phase 1 scope: parses every `potential` piece and checks that,
    together, they tile `bspline.domain` with no gaps and no ambiguous
    overlaps (validate_potential_tiling) -- the one piece of real
    validation logic this module implements.
    """
    try:
        run_cfg = cfg["run"]
    except KeyError as e:
        raise ConfigValidationError("config missing required 'run' section") from e

    if not isinstance(run_cfg, dict):
        raise ConfigValidationError(f"'run' must be a mapping, got {type(run_cfg).__name__}")

    for key in ("run_tise", "output_dir"):
        if key not in run_cfg:
            raise ConfigValidationError(f"config missing required 'run.{key}'")

    try:
        potential_cfg = cfg["potential"]
    except KeyError as e:
        raise ConfigValidationError("config missing required 'potential' section") from e

    if not isinstance(potential_cfg, list) or not potential_cfg:
        raise ConfigValidationError("'potential' must be a non-empty list of piece strings")

    pieces = [parse_potential_piece(p) for p in potential_cfg]

    try:
        domain = cfg["bspline"]["domain"]
    except (KeyError, TypeError) as e:
        raise ConfigValidationError("config missing required 'bspline.domain'") from e

    try:
        domain_tuple = (float(domain[0]), float(domain[1]))
    except (TypeError, ValueError, IndexError) as e:
        raise ConfigValidationError(f"'bspline.domain' must be a [lo, hi] pair: {domain!r}") from e

    validate_potential_tiling(pieces, domain_tuple)


# ─── Interval / potential-piece parsing ────────────────────────────────────


def parse_interval(domain: str) -> Interval:
    """Parse interval notation, e.g. "[0,100]", "(0, 100]", "[0, inf)".

    '[' / ']' denote an inclusive (closed) bound; '(' / ')' denote an
    exclusive (open) bound. `inf` / `-inf` are valid bounds and parse via
    Python's native float() -- no special-casing needed. Raises
    ConfigValidationError on any malformed input.
    """
    if not isinstance(domain, str):
        raise ConfigValidationError(f"interval must be a string, got {type(domain).__name__}: {domain!r}")

    s = domain.strip()
    if len(s) < 3 or s[0] not in "[(" or s[-1] not in "])":
        raise ConfigValidationError(f"malformed interval (expected e.g. '[a, b)'): {domain!r}")

    lo_inclusive = s[0] == "["
    hi_inclusive = s[-1] == "]"

    parts = s[1:-1].split(",")
    if len(parts) != 2:
        raise ConfigValidationError(f"malformed interval (expected exactly one comma): {domain!r}")

    lo_str, hi_str = parts[0].strip(), parts[1].strip()
    try:
        lo = float(lo_str)
        hi = float(hi_str)
    except ValueError as e:
        raise ConfigValidationError(f"malformed interval bound in {domain!r}: {e}") from e

    if lo > hi:
        raise ConfigValidationError(f"interval lower bound exceeds upper bound: {domain!r}")

    return Interval(lo=lo, lo_inclusive=lo_inclusive, hi=hi, hi_inclusive=hi_inclusive)


def parse_potential_piece(piece: str) -> dict:
    """Parse one `potential` list entry: a string encoding a Python dict
    literal with 'domain' and 'function' keys, e.g.

        "{'domain': '(0, 100]', 'function': '-1/x + 1/x^2'}"

    Parsed with ast.literal_eval -- NOT json.loads. This single-quoted
    Python-dict-literal format is the current config.yaml schema (see
    docs/superpowers/specs/2026-06-28-config-yaml-schema-design.md and
    the real config.yaml), distinct from an older, now-retired prototype
    (formerly input.py) that used JSON-double-quoted strings.
    """
    if not isinstance(piece, str):
        raise ConfigValidationError(f"potential piece must be a string, got {type(piece).__name__}: {piece!r}")

    try:
        parsed = ast.literal_eval(piece)
    except (ValueError, SyntaxError) as e:
        raise ConfigValidationError(f"malformed potential piece (expected a dict literal): {piece!r}: {e}") from e

    if not isinstance(parsed, dict):
        raise ConfigValidationError(f"potential piece did not parse to a dict: {piece!r}")

    if "domain" not in parsed or "function" not in parsed:
        raise ConfigValidationError(f"potential piece missing 'domain' or 'function' key: {piece!r}")

    return parsed


def validate_potential_tiling(pieces: list[dict], domain: tuple[float, float]) -> None:
    """Check that `pieces` collectively tile `domain` = (lo, hi) with no
    gaps and no ambiguous overlaps.

    Pieces are sorted by lower bound. Only *coverage* of the domain edges
    is required: the lowest piece's lo <= domain[0] + EPS, the highest
    piece's hi >= domain[1] - EPS. The first piece is deliberately NOT
    required to be inclusive at domain[0] -- e.g. the real config.yaml's
    single Coulomb piece is "(0, 100]" against a [0.0, 100.0] bspline
    domain, open at x=0 because -1/x is singular there.

    Adjacent pieces must meet with no numeric gap (abs(cur.hi - nxt.lo) <=
    EPS). At a shared boundary point, it is only a genuine (ambiguous)
    overlap if BOTH sides are inclusive there -- both sides exclusive at a
    shared point is tolerated as a measure-zero gap, not numerically
    meaningful for a continuous potential.
    """
    if not pieces:
        raise ConfigValidationError("'potential' must contain at least one piece")

    domain_lo, domain_hi = domain

    entries = [(parse_interval(piece["domain"]), piece) for piece in pieces]
    entries.sort(key=lambda entry: entry[0].lo)

    first_interval, first_piece = entries[0]
    if first_interval.lo > domain_lo + EPS:
        raise ConfigValidationError(
            f"potential does not cover domain lower bound {domain_lo}: "
            f"lowest piece {first_piece['domain']!r} starts at {first_interval.lo}"
        )

    last_interval, last_piece = entries[-1]
    if last_interval.hi < domain_hi - EPS:
        raise ConfigValidationError(
            f"potential does not cover domain upper bound {domain_hi}: "
            f"highest piece {last_piece['domain']!r} ends at {last_interval.hi}"
        )

    for (cur_interval, cur_piece), (nxt_interval, nxt_piece) in zip(entries, entries[1:]):
        if abs(cur_interval.hi - nxt_interval.lo) > EPS:
            raise ConfigValidationError(
                f"gap or overlap between potential pieces {cur_piece['domain']!r} and "
                f"{nxt_piece['domain']!r}: {cur_interval.hi} does not meet {nxt_interval.lo}"
            )
        if cur_interval.hi_inclusive and nxt_interval.lo_inclusive:
            raise ConfigValidationError(
                f"overlapping potential pieces {cur_piece['domain']!r} and "
                f"{nxt_piece['domain']!r}: both are inclusive at shared boundary {cur_interval.hi}"
            )


# ─── TISE solver subprocess ────────────────────────────────────────────────


def _run_stage(stage_name: str, cmd: list[str]) -> subprocess.CompletedProcess:
    try:
        return subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        raise SolverStageError(f"{stage_name} failed (exit {e.returncode}).\n"
                                f"command: {' '.join(cmd)}\n--- stderr ---\n{e.stderr}") from e
    except OSError as e:
        # Covers a missing binary (FileNotFoundError), a binary that exists
        # but isn't executable (PermissionError), and other exec-time
        # failures (e.g. NotADirectoryError, "Exec format error") -- all
        # OSError subclasses raised by subprocess.run() itself before the
        # child process even starts, as opposed to CalledProcessError,
        # which is raised after the child ran and exited non-zero. str(e)
        # already carries the OS's own reason (e.g. "No such file or
        # directory" vs "Permission denied"), so it reads sensibly for
        # either case.
        raise SolverStageError(f"{stage_name}: could not execute {cmd[0]}: {e}") from e


def run_tise_solver(config_path: str, tise_output_dir: Path, binary: Path = DEFAULT_TISE_SOLVER) -> None:
    """Invoke the TISE solver subprocess per docs/SDD.md Sec 7.2.1."""
    try:
        tise_output_dir.mkdir(parents=True, exist_ok=True)
    except OSError as e:
        raise SolverStageError(f"TISE solver: could not create output directory {tise_output_dir}: {e}") from e

    _run_stage(
        "TISE solver",
        [str(binary), "--config", config_path, "--output-dir", str(tise_output_dir)],
    )


# ─── Analysis subprocess ───────────────────────────────────────────────────


def run_analysis_stage(
    config_path: str,
    tise_dir: Path,
    tdse_dir: Path,
    script: Path = DEFAULT_ANALYSIS_SCRIPT,
) -> None:
    """Invoke the Analysis subprocess per docs/SDD.md Sec 7.2.3.

    Unlike run_tise_solver, this does NOT create tise_dir/tdse_dir -- both
    are Analysis's INPUT paths (read, not written; see ADR-0005,
    docs/adr/0005-defer-analysis-output-artifact-format.md), and a missing
    one is Analysis's own contract to handle (Sec 5.4.4: absent tdse_dir
    tolerated; absent/malformed tise_dir a hard TiseOutputError), not
    something the Controller should paper over.
    """
    _run_stage(
        "Analysis",
        [sys.executable, str(script), "--config", config_path,
         "--tise-dir", str(tise_dir), "--tdse-dir", str(tdse_dir)],
    )


def read_warnings(tise_output_dir: Path) -> list[dict]:
    """Best-effort read of tise_output_dir/warnings.json.

    Degrades gracefully on any problem -- missing file, unreadable file,
    malformed JSON, or unexpected JSON shape -- printing a note to stderr
    and returning [] rather than raising. This sidecar is best-effort:
    nothing downstream depends on it yet, so a missing/malformed file
    should not turn an otherwise-successful physics run into a Controller
    failure.
    """
    path = tise_output_dir / "warnings.json"
    try:
        with path.open("r") as f:
            data = json.load(f)
    except OSError:
        print(f"controller.py: warning: {path} not found; skipping", file=sys.stderr)
        return []
    except ValueError as e:
        # json.JSONDecodeError is a ValueError subclass -- the file exists
        # and was read, but its contents aren't valid JSON.
        print(f"controller.py: warning: {path} has invalid contents ({e}); skipping", file=sys.stderr)
        return []

    if not isinstance(data, list):
        print(
            f"controller.py: warning: {path} has invalid contents "
            f"(expected a JSON array, got {type(data).__name__}); skipping",
            file=sys.stderr,
        )
        return []

    if not all(isinstance(w, dict) for w in data):
        print(
            f"controller.py: warning: {path} has invalid contents "
            f"(expected a JSON array of objects); skipping",
            file=sys.stderr,
        )
        return []

    return data


def print_warnings(warnings: list[dict], stream=None) -> None:
    """Print each {"category": ..., "message": ...} warning dict to `stream`.

    `stream` defaults to None (resolved to sys.stderr inside the function
    body) rather than defaulting to `sys.stderr` directly in the signature:
    a default argument value is bound once, at function-definition time, so
    `stream=sys.stderr` in the signature would capture whatever object
    `sys.stderr` pointed to at import time and keep writing to it even
    after something (e.g. pytest's capsys, contextlib.redirect_stderr)
    reassigns `sys.stderr` later. Resolving inside the body does a fresh
    lookup on every call, so redirection works as expected.
    """
    if stream is None:
        stream = sys.stderr
    for w in warnings:
        print(f"[{w.get('category', 'unknown')}] {w.get('message', '')}", file=stream)


# ─── Orchestration & CLI ───────────────────────────────────────────────────


def run(config_path: str) -> None:
    cfg = load_config(config_path)
    validate_config(cfg)

    if cfg["run"].get("run_tdse"):
        raise ControllerError("run.run_tdse is not yet supported (Phase 5); set run.run_tdse: false")

    output_dir = Path(cfg["run"]["output_dir"])
    tise_dir = output_dir / "tise"
    tdse_dir = output_dir / "tdse"

    if cfg["run"]["run_tise"]:
        run_tise_solver(config_path, tise_dir)
        print_warnings(read_warnings(tise_dir))

    # Analysis and TISE are independently-gated per SDD Sec 7.2.2/Figure 2; a
    # TISE failure above already aborts run() via the exception, so no
    # explicit "and tise succeeded" check is needed here.
    if cfg["run"].get("run_analysis"):
        run_analysis_stage(config_path, tise_dir, tdse_dir)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="controller.py",
        description="Controller for the 1D QM Playground: orchestrates solver subprocesses per config.yaml.",
    )
    parser.add_argument("--config", default="config.yaml", help="path to config.yaml (default: %(default)s)")
    args = parser.parse_args(argv)

    try:
        run(args.config)
    except ControllerError as e:
        print(f"controller.py: {e}", file=sys.stderr)
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
