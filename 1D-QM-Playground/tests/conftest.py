"""Shared pytest configuration for the 1D QM Playground test suite.

`controller.py` is a flat script at the project root, not an installed
package, so it isn't importable as `import controller` unless the project
root is on sys.path. This shim makes that import work regardless of the
directory pytest is invoked from.

This file also hosts the fixtures shared by the real-subprocess integration
suite (tests/test_controller_integration.py, Task 5): `tise_solver_binary`
(builds/locates the real tise_solver binary) and `tmp_config` (a tmp_path-
scoped copy of the real config.yaml with run.output_dir redirected away from
the real project's ./data). test_controller_unit.py does not use either
fixture -- it mocks subprocess.run and never touches a real binary or the
real config.yaml.
"""

import subprocess
import sys
from pathlib import Path

import pytest
import yaml

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))


@pytest.fixture(scope="session")
def tise_solver_binary() -> Path:
    """Build the real tise_solver binary via CMake if not already present, return its path."""
    tise_dir = PROJECT_ROOT / "TISE"
    build_dir = tise_dir / "build"
    binary = build_dir / "tise_solver"
    if not binary.exists():
        subprocess.run(["cmake", "-S", str(tise_dir), "-B", str(build_dir)], check=True)
        subprocess.run(["cmake", "--build", str(build_dir), "--target", "tise_solver"], check=True)
    return binary


@pytest.fixture
def tmp_config(tmp_path: Path) -> Path:
    """Copy the real config.yaml but redirect run.output_dir into pytest's tmp_path.

    CRITICAL test-hygiene note: the real config.yaml has run.output_dir:
    "./data". Every integration test that actually runs the solver MUST go
    through this fixture (or an equivalent tmp_path-scoped copy) rather than
    the raw config.yaml path -- otherwise a solver run would write into the
    real project's data/ directory on every test run.
    """
    with open(PROJECT_ROOT / "config.yaml") as f:
        cfg = yaml.safe_load(f)
    cfg["run"]["output_dir"] = str(tmp_path / "data")
    out = tmp_path / "config.yaml"
    with open(out, "w") as f:
        yaml.safe_dump(cfg, f)
    return out
