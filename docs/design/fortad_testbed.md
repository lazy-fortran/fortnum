# fortad as a second `autodiff` backend

This branch adds [fortad](https://github.com/lazy-fortran/fortad) as a source
generator for derivative kernels, alongside the existing Enzyme path.

## What changes for fortnum

Nothing at build time or run time. fortad is a **build-time source generator**:
it reads plain Fortran and writes plain Fortran, the generated files are
committed under `src/generated/`, and they compile with any conforming
compiler. fortnum gains no link-time dependency, no compiler plugin, and no
LLVM version constraint.

That is the whole difference from the Enzyme path, which needs flang-new or
LFortran plus a `ClangEnzyme-NN.so` matched to an exact LLVM.

## Layout

| Path | What it is |
|---|---|
| `tools/fortad/dot_sin_kernel.f90` | the primal, unannotated |
| `tools/fortad/generate.sh` | regenerates the derivative kernels |
| `src/generated/fortnum_dot_sin_jvp.f90` | scalar JVP, generated |
| `src/generated/fortnum_dot_sin_jvp_v.f90` | vector JVP, generated |
| `src/generated/fortnum_dot_sin_vjp.f90` | VJP (gradient), generated |
| `test/fortad/test_fortad_dot_sin.f90` | acceptance test |

Regenerate with:

```bash
FORTAD_REPO=../fortad tools/fortad/generate.sh
```

## How the candidate is judged

By fortnum's existing rules, unmodified. A generated `autodiff` candidate
clears the same bar as any other:

1. **Independent oracle first.** `test/fortad/test_fortad_dot_sin.f90` checks
   the tangent against central finite differences with a step-size convergence
   requirement, and checks that vector mode reproduces the scalar tangent
   direction by direction. The gradient is checked against the analytical
   gradient of this kernel and against the adjoint identity with the already
   verified tangent. Nothing here compares fortad against another AD tool;
   agreement with Enzyme would be corroboration, never the oracle.
2. **Then a measurement.** Complete-workload wall clock decides. The
   cross-engine timings live in
   [fortad-bench](https://github.com/lazy-fortran/fortad-bench), not here,
   because they need multiple toolchains and must not gate a fortnum commit.

The mechanism's name selects nothing. If Enzyme or a hand-written analytical
kernel wins a workload, that one is selected.

## Why the vector form is generated too

fortnum's derivative products - linear UQ, sensitivity analysis, Gauss-Newton -
almost never want exactly one direction. Vector mode carries `n_dir` tangents
through one primal sweep with the direction axis leading, so it is the
contiguous one:

```fortran
real(8), intent(in), dimension(n_dir, n) :: a_d
```

Measured on the equivalent benchmark kernel, per-direction cost falls from
about 9.6 ns per element at one direction to about 0.9 at sixteen, while any
engine called once per direction stays flat. Those are numbers from one
machine, recorded in fortad-bench, not a promise about another.

## Current limits

fortad refuses what it cannot do correctly, by name. For reverse mode today
that means nonlinear loop-carried recurrences, nested loops, and branches
inside loops. Forward mode handles all of those.

Measured against Enzyme on the equivalent benchmark kernel, one machine:
forward mode about 8% faster at one direction and roughly 10x per direction at
sixteen; reverse mode 1.7-1.8x faster; build time 2.7x faster. Raw data in
fortad-bench. Selection still follows fortnum's rules - measured
complete-workload wall clock on the real workload, not these microbenchmarks.
