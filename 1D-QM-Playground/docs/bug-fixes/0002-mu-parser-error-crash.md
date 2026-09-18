# BUG-0002: `mu::ParserError` Crashed the Whole Process Instead of Failing Cleanly

- **Status:** Fixed
- **Date introduced:** 2026-07-26 (`d01efad`, "Gracedev (#7)")
- **Date fixed:** 2026-09-18 (`289a7da`, "fix(TISE): catch mu::ParserError from malformed potential expressions")
- **Report date:** 2026-09-18

## Symptom

A malformed `function` expression in a config's `potential` block (unbalanced
parens, a stray operator, an undefined identifier, etc.) crashed the whole
`tise_solver` (and separately, `H-BoundStates`) process:

```
terminate called after throwing an instance of 'mu::ParserError'
Aborted
```

instead of the intended `tise_solver: <message>` stderr line plus a clean
non-zero exit every other config error gets.

## Root Cause

`mu::Parser` is constructed in exactly one place in the entire repo —
`evaluateFunction` (`TISE/tise.cpp`). Its `SetExpr()`/`Eval()` calls throw
`mu::ParserError` (`/usr/include/muParserError.h`), a standalone class that
does **not** derive from `std::exception`. It therefore escaped every
`catch (const std::exception &)` in the codebase — the top-level one in
`tise_solver_main.cpp`, a local one wrapping an explicit `evaluateFunction`
call in the same file, and two in `main.cpp` (`H-BoundStates`) — reaching
`std::terminate()`.

The earliest point a bad expression could crash a normal `tise_solver` run
was `classifyAsymptote` (called from `buildStrategicGridAndDropSet`),
invoked right after the potential map was parsed and overlap-validated, and
before any output file was opened. A piece whose domain a given run's own
grid/asymptote probing never happened to reach could go unvalidated
entirely.

## Fix

Two parts, both in `TISE/tise.cpp`:

- `evaluateFunction`'s `SetExpr`/`Eval` calls are now wrapped in
  `try { ... } catch (const mu::ParserError &e)`, re-thrown as
  `std::runtime_error` naming the offending domain, expression, and
  muParser's own message (e.g. `malformed potential expression in domain
  '(0, inf)': 'x +* 2': Unexpected operator "*" found at position 3`). Since
  `evaluateFunction` is the sole chokepoint, this one change fixes the crash
  for both `tise_solver` and `H-BoundStates` with no `main.cpp` change
  needed.
- New `validatePotentialExpressionsParse`, wired into
  `tise_solver_main.cpp` immediately after the existing
  `validateNoOverlappingPotentialPieces` call: it evaluates one
  representative finite point per potential piece, once, right after config
  load — so a malformed expression anywhere is caught before
  `buildStrategicGridAndDropSet` (or any other solve work) begins, even in
  a piece a given run's own probing would never reach.

## Verification

- C++ (`TISE/tests/test_tise.cpp`): `EvaluateFunctionTest.MalformedExpression*`
  (proves the throw-type contract and that the error names the offending
  domain/expression) and a new `ValidatePotentialExpressionsParseTest` suite
  (covers an interior piece, an infinite-bound piece, and a narrow piece a
  solve-driven probe might miss).
- Python (`tests/test_controller_integration.py`):
  `test_malformed_potential_expression_fails_cleanly_not_sigabrt` drives the
  real compiled `tise_solver` binary and asserts the failure mode itself —
  `SolverStageError` reports `exit 1` (not a negative/signal return code)
  and stderr contains the clean `tise_solver: ` prefix but *not* the abort
  markers (`terminate called after throwing`, `Aborted`). Verified by
  temporarily reverting the C++ fix and confirming this test fails for the
  expected reason (`exit -6`, raw abort banner in stderr) before restoring
  it — a bare `pytest.raises(SolverStageError)` alone would have passed
  against either the crashing or the fixed binary, since SIGABRT's
  `returncode == -6` is nonzero too.
- Full regression: `ctest --test-dir TISE/build` (4/4 suites, 211 tests) and
  `pytest tests/` (219 tests) both green, confirming the new eager
  validation rejects nothing that was previously valid.
- Manual: both `tise_solver` and `H-BoundStates`, run directly against a
  malformed expression, now exit 1 with a clean message and no partial
  output, instead of aborting.

## Source

`TISE/tise.cpp` (`evaluateFunction`, `validatePotentialExpressionsParse`),
`TISE/tise.hpp`, `TISE/tise_solver_main.cpp`; tests in
`TISE/tests/test_tise.cpp` and `tests/test_controller_integration.py`; docs
in `docs/guides/input-file-authoring-guide.md` and `docs/SDD.md`. Introduced
by `d01efad`; fixed by `289a7da`.
