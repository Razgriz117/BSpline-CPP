# Bug-Fix Reports

This directory records notable bugs after the fact: what broke, why, how it
was fixed, and how the fix was verified. Unlike the [ADRs](../adr/), these
are not design decisions — they document root-cause diagnoses worth keeping
so a future regression to the same failure mode is recognizable.

| ID | Title | Status | Summary |
|---|---|---|---|
| [0001](0001-continuum-state-coefficient-basis.md) | `writeContinuumInfo` Evaluated Continuum States in the Wrong Coefficient Basis | Fixed | Continuum-state plots were physically wrong at every energy because `writeContinuumInfo` called `bs.eval` on confined-eigenstate-basis coefficients instead of transforming them to B-spline coefficients first, as `matchAsymptotic` already did. |
| [0002](0002-mu-parser-error-crash.md) | `mu::ParserError` Crashed the Whole Process Instead of Failing Cleanly | Fixed | A malformed `function` expression escaped every `catch (const std::exception&)` uncaught (mu::ParserError doesn't derive from std::exception), crashing the process via `std::terminate()` instead of a clean error exit. |

## Format

Each report uses the same fields: **Status**, **Date introduced**, **Date
fixed**, **Report date**, **Symptom**, **Root Cause**, **Fix**,
**Verification**, and **Source**.
