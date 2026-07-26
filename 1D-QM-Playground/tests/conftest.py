"""Shared pytest configuration for the 1D QM Playground test suite.

`controller.py` is a flat script at the project root, not an installed
package, so it isn't importable as `import controller` unless the project
root is on sys.path. This shim makes that import work regardless of the
directory pytest is invoked from.

NOTE: a later task (integration tests, Task 5) extends this file with
additional fixtures (e.g. building the real tise_solver binary, rewriting
config.yaml's output_dir into a tmp path). Keep unit-test-only concerns out
of this file -- it is shared by both test_controller_unit.py and the
upcoming integration-test module.
"""

import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PROJECT_ROOT))
