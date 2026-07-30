# fortnum

`fortnum` is a Fortran 2018 numerical library for special functions,
quadrature, adaptive integration, FFTs, ODEs, roots, interpolation, splines,
linear algebra, and random numbers. The library has no mutable global numerical
state. LAPACK and BLAS are its only runtime numerical dependencies.

The derivative API is product-oriented. Kernels may provide JVPs, VJPs,
gradients, or HVPs through competing `autodiff`, `analytical`, and `hybrid`
candidates. Correctness comes from independent oracles. Complete-workload wall
clock and peak memory select production candidates.

## Build

Requirements:

- a Fortran 2018 compiler
- CMake 3.22 or newer, or `fpm`
- BLAS and LAPACK

With CMake:

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build
```

With the local `fo` workflow:

```bash
fo
fo install
```

`fo` runs static checks, builds, tests, and lints the fpm configuration. See
[CONTRIBUTING.md](CONTRIBUTING.md) for the complete contributor gate.

Useful CMake options:

| Option | Default | Effect |
| --- | --- | --- |
| `FORTNUM_BUILD_TESTING` | on for a top-level build | Build the test suite |
| `FORTNUM_BUILD_EXAMPLES` | on for a top-level build | Build API examples |
| `FORTNUM_CHECK_GENERATED` | off | Regenerate and compare committed `fortsym` kernels |
| `FORTNUM_USE_BENCHMARK_SELECTION` | on | Read committed derivative selections |
| `FORTNUM_ENABLE_ENZYME` | off | Build supported CPU Enzyme tests |
| `FORTNUM_ENZYME_REQUIRED` | off | Fail if the requested Enzyme toolchain is unavailable |
| `FORTNUM_GPU_BACKEND` | `NONE` | Select `NONE`, `OPENACC`, or `OPENMP` |

Selecting an unavailable GPU backend fails configuration. GPU execution never
falls back silently to the host. OpenMP GPU builds currently require the
validated NVIDIA HPC SDK `nvfortran` compiler.

## Use

Each numerical area has its own module:

```fortran
program example
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special, only: dawson, dawson_jvp
    implicit none

    real(dp) :: value, tangent

    value = dawson(0.5_dp)
    call dawson_jvp(0.5_dp, 1.0_dp, tangent)
    print *, value, tangent
end program example
```

The main public modules are:

| Area | Modules |
| --- | --- |
| kinds and status | `fortnum_kinds`, `fortnum_status` |
| special functions | `fortnum_special` and its domain modules |
| quadrature and integration | `fortnum_quadrature`, `fortnum_integrate`, `fortnum_integrate_gk`, `fortnum_cquad`, `fortnum_levin` |
| transforms | `fortnum_fft` |
| ODEs | `fortnum_ode` and method modules |
| roots | `fortnum_roots`, `fortnum_multiroot`, `fortnum_roots_complex` |
| interpolation | `fortnum_interp`, `fortnum_polynomial`, `fortnum_bspline` |
| linear algebra | `fortnum_linalg`, `fortnum_krylov` |
| random numbers | `fortnum_rng` |

The [API guide](docs/api.md) shows the public families and common call
patterns. Production module sources remain authoritative for exact signatures.
The installed C interface is declared in `include/fortnum.h`.

## Differentiation

Public terminology is fixed:

- `autodiff` means compiler or source-transformation differentiation
- `analytical` includes explicit formulas, stable recurrences, linear
  operators, frozen traces, and implicit tangent or adjoint solves
- `hybrid` means autodiff composition with at least one analytical operator
  rule

A procedure can have several candidates for the same product. Finite
differences are validation diagnostics unless measurements justify production
use. Residual-defined outputs always admit an analytical implicit candidate.

`fortsym` owns build-time symbolic algebra and code generation. The lock file
at `tools/codegen/fortsym.lock` identifies the tested revision. Generated
production kernels are committed under `src/generated`; temporary algebraic
variants and Enzyme wrappers stay in the build tree. Generators are quiet by
default. Set `FORTNUM_CODEGEN_VERBOSE=1` for diagnostic output.

CPU support includes selected `autodiff`, `analytical`, and `hybrid`
candidates. GPU support currently uses validated `analytical` leaves:

| Target | `autodiff` | `analytical` | `hybrid` |
| --- | --- | --- | --- |
| CPU | supported selectively | supported | supported selectively |
| GPU | unsupported | supported for listed kernels | unsupported |

The exact contracts and evidence live in:

- [derivative contract](docs/design/ad.md)
- [implementation plan](docs/design/differentiation_plan.md)
- [kernel ownership inventory](docs/design/derivative_kernel_inventory.csv)
- [CPU/GPU contract](docs/design/gpu.md)
- [current evidence report](docs/design/differentiation_report.md)
- [ROADMAP.md](ROADMAP.md)

## Performance evidence

Machine-readable measurements are committed under `benchmark/reference`.
Normalized cross-kernel statistics live in
`benchmark/report/data/mechanism_tournaments.csv`. Reproducible `fortplot`
programs generate the report figures; generated PNG files are not committed.

Every selected derivative record includes validation, workload shape, compiler
and hardware identity, wall-clock time, and peak memory. Records add direction,
input/output scaling, cache counters, code size, or GPU profiling data when
those measurements affect the decision.

See [benchmark/README.md](benchmark/README.md) to reproduce the reports and add
new evidence.

## Documentation

[docs/README.md](docs/README.md) is the documentation map. It distinguishes
normative contracts, user guides, migration notes, plans, and measured
evidence. This keeps implementation history out of the README.

## License

MIT. See [LICENSE](LICENSE).
