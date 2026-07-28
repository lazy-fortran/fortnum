# Contributing

`fortnum` accepts small, independently validated changes. Read
[AGENTS.md](AGENTS.md) and the relevant design contract before editing a
numerical module.

## Workflow

1. Select one unchecked item in [ROADMAP.md](ROADMAP.md), or one narrowly
   scoped issue.
2. Make the smallest complete change.
3. Add an independent behavioral oracle.
4. Update the affected user, design, and benchmark documentation.
5. Run focused checks while developing.
6. Run the full gate before committing.
7. Check off the ROADMAP item, commit, push `main`, and install.

The full local gate is:

```bash
fo
cmake -S . -B build
cmake --build build
ctest --test-dir build --output-on-failure
```

Code-generation changes also require:

```bash
cd tools/codegen
./check_generated.sh
```

Use `fo test <name>` and `fo exec <target>` for focused Fortran work. Do not
invoke build-tree test binaries directly.

## Numerical changes

Tests must compare behavior with an independent source:

- an exact formula or manufactured solution
- an adjoint dot-product identity
- a residual linearization
- high-precision or published reference data
- central finite differences or complex step when their assumptions hold

A check that restates the implementation is not an oracle. Cover normal
inputs, boundaries, regime changes, singular cases, and documented failure
status where applicable.

Preserve stable floating-point semantics. Fast-math, reassociation, approximate
functions, and relaxed contraction need a stated accuracy contract plus a
measured complete-workload benefit.

## Derivative changes

Use `autodiff`, `analytical`, and `hybrid` as defined in
[docs/design/ad.md](docs/design/ad.md). Register competing candidates per
product. A mechanism name never selects a winner.

Document:

- the mathematical operator
- active and inactive arguments
- requested product (`JVP`, `VJP`, gradient, or `HVP`)
- all admissible candidates
- the independent validation oracle
- workload dimensions and reusable primal state
- median wall clock, dispersion, and candidate-specific peak memory
- input, output, and direction scaling when it can change forward/reverse
  selection
- cache or device counters when they explain complete-workload wall clock

Use an analytical implicit candidate for residual-defined outputs. Keep solver
iterations, adaptive decisions, and index searches outside generic autodiff
unless a measured candidate deliberately differentiates a fixed trace.

Use `../lazy-fortran/fortsym` for symbolic algebra and generated Fortran.
Generated code handles local explicit algebra and mechanical Enzyme
boundaries. Stable recurrences, solver orchestration, traces, and independent
oracles remain hand-written. Extend and test `fortsym` in its own repository
when its supported interface is insufficient.

## Fortran style

- free-form Fortran 2018
- `implicit none` in every scope
- explicit `intent` on every dummy argument
- declarations at the start of a scope
- `use ..., only:` before `implicit none`
- `real(dp)` with `dp => real64`
- `_t` suffix for derived types
- `allocatable` in preference to pointers
- no mutable global numerical state
- no allocation inside hot loops
- 4-space indentation and an 88-column formatting target

Fortran `.and.` and `.or.` do not short-circuit. Split guards that protect
indexing, allocation, optional arguments, pointer association, or descriptors.

## Pull requests and issues

Use the repository templates. A pull request identifies its oracle, commands
run, benchmark evidence, documentation changes, and generated-source status.
Do not weaken or skip a failing test.
