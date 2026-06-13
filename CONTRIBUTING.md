# Contributing to fortnum

fortnum is a clean-room Fortran numerical library. Every public procedure is
primal-first and derivative-ready. This document describes how to contribute,
what the milestone structure is, and what contract a new or ported module must
satisfy.

## Milestones

| Milestone | Scope |
|-----------|-------|
| M0 | Infrastructure: CMake, fpm, CI skeleton, code style, AD policy framework |
| M1 | Special-function kernels (Bessel, gamma, error function, …) |
| M2 | ODE solvers |
| M3 | Roots, interpolation, random number generation |
| M4 | Adaptive integration and quadrature |
| M5 | Hardening: edge cases, precision audit, portability |
| M6 | Differentiability: Enzyme wiring, oracle suite, AD policy review |

Work in a feature branch. Name it `m<N>/<short-description>`, e.g.
`m1/bessel-j0`. Open a PR against `main` when tests pass.

## Build

CMake is the primary build system:

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure
```

`fpm` also works via `fo`. Run `fo` (no arguments) for the full pipeline:
static analysis, build, test, lint, format check.

## Code style

- Free source form. `implicit none` in every scoping unit.
- `use <mod>, only: ...` before `implicit none`.
- `real(dp)` via `use, intrinsic :: iso_fortran_env, only: dp => real64`.
- All dummy arguments have explicit `intent`. Declarations at start of scope.
- Derived-type names end in `_t`. `allocatable` over pointers.
- No module-level `save` or global mutable state.
- `fprettify`: 88-column width, 4-space indent.
- Comments say why a choice was made, not what the next line does.

## Per-module contract

Every public procedure must satisfy the following contract before a PR merges.

### 1. Derivative policy

Each public procedure declares a derivative policy in its module-level doc
comment. The canonical policy names and their semantics are defined in
`docs/design/ad.md` (issue #37). The four policies are:

- `transparent`: no branches on active variables; Enzyme differentiates
  through the routine without annotation.
- `analytic_rule`: a hand-coded derivative is registered alongside the
  primal (e.g. a known recurrence for a special-function derivative).
- `implicit_rule`: the derivative follows from an implicit equation or
  identity (e.g. Newton step for a root solver).
- `trace_rule`: requires source transformation; the routine is tagged for
  Enzyme but needs a driver in `src/ad/`.
- `primal_only`: differentiation is not supported; must be documented with
  a justification.

Do not invent policy names. See `docs/design/ad.md` for the authoritative
definitions and the annotation syntax.

### 2. Active vs. inactive arguments

The module doc comment or the procedure doc comment must list which dummy
arguments are *active* (carry derivative information) and which are *inactive*
(parameters, indices, tolerances). Enzyme infers activity by data flow, but
the explicit declaration is the contract for reviewers and for M6 wiring.

### 3. Primal purity

Write primal routines `pure` wherever the algorithm allows it. Avoid I/O,
internal state, and pointer aliasing in the primal path. A `pure` primal is a
prerequisite for `transparent` policy; it also keeps oracle tests deterministic.

### 4. Tests

Every numerical module requires:

- **Unit tests** under `test/<domain>/` covering normal inputs and documented
  edge cases.
- **Oracle tests** under `test/oracle/` that compare the routine against an
  independent reference (scipy, DLMF tables, or an analytic result). The
  oracle test is the ground truth for both value and derivative correctness.
  Use `test-drive` assertions with explicit tolerances.

Tests run via CTest. Do not merge if any test fails.

## Filing issues

Use the issue templates in `.github/ISSUE_TEMPLATE/`. Two templates exist:

- **Module port**: for porting a known algorithm from GSL or another
  reference library.
- **Clean-room implementation**: for writing a new routine from primary
  sources (DLMF, textbooks, papers) without reference to existing code.

Both templates ask for the derivative policy up front. Fill it in, or mark it
`TBD` and open a follow-up referencing `docs/design/ad.md`.

## Pull requests

Use the pull request template. A PR must include:

- a reference to the issue it closes,
- the derivative policy for each new public procedure,
- real `ctest` output showing all tests pass.

Do not open a PR that skips, weakens, or disables existing tests.
