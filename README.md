# fortnum

A Fortran 2018 numerical library: special functions, quadrature and adaptive
integration, FFTs, ODEs, root finding, interpolation and splines, linear algebra,
random numbers and sampling. It is written to replace GSL-, QUADPACK- and
FFTW-style dependencies in Fortran research codes, and provides derivative
products (JVPs, VJPs, gradients) for many of its kernels.

Status: research code under active development; interfaces may change.

## Features

- Special functions: Bessel (real and complex), gamma, error and Faddeeva
  functions, Dawson, elliptic, Legendre and toroidal, Coulomb, hypergeometric
  1F1 and others.
- Quadrature, adaptive (Gauss-Kronrod, CQUAD, Levin) integration, FFTs.
- ODE integrators (DOP853, RK8PD, VODE, DDEABM, Gauss-Radau, extrapolation,
  geometric integrators) with events and sensitivities.
- Roots and multidimensional roots, B-splines, Krylov solvers and matrix-free
  eigenproblems, Toeplitz and tensor-product operators.
- Random numbers with explicit state, NUTS sampling, Gaussian-process
  multi-fidelity prediction.
- Verified computing: interval and ball arithmetic and rigorous enclosures.
- No mutable global numerical state; BLAS and LAPACK are the only runtime
  numerical dependencies. A C interface is declared in `include/fortnum.h`.

## Building

Requirements: a Fortran 2018 compiler, CMake 3.22 or newer (or fpm), and
BLAS/LAPACK.

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build
```

With fpm (Python 3.11 or newer supplies the per-test data arguments):

```bash
fpm test --profile release --runner "python3 scripts/fpm_test_runner.py"
```

CMake options, GPU backends, quad-precision support and macOS notes are in
[docs/building.md](docs/building.md).

## Usage

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

| Area | Modules |
| --- | --- |
| kinds and status | `fortnum_kinds`, `fortnum_status` |
| special functions | `fortnum_special` and its domain modules |
| quadrature and integration | `fortnum_quadrature`, `fortnum_integrate`, `fortnum_integrate_gk`, `fortnum_cquad`, `fortnum_levin` |
| transforms | `fortnum_fft` |
| ODEs | `fortnum_ode` and method modules |
| roots | `fortnum_roots`, `fortnum_multiroot`, `fortnum_roots_complex` |
| interpolation | `fortnum_interp`, `fortnum_polynomial`, `fortnum_bspline`, `fortnum_bspline_lsq` |
| linear algebra | `fortnum_linalg`, `fortnum_cholesky`, `fortnum_krylov`, `fortnum_symmetric_eigen`, `fortnum_tensor_product`, `fortnum_toeplitz` |
| random numbers | `fortnum_rng` |
| verified computing | `fortnum_rounding`, `fortnum_interval`, `fortnum_interval_qp`, `fortnum_ball`, `fortnum_idual`, `fortnum_cfft_rigorous`, `fortnum_fseries`, `fortnum_verified_linalg` |

The [API guide](docs/api.md) shows the module families and common call patterns.

## Documentation

- [Documentation map](docs/README.md)
- [API guide](docs/api.md)
- [Building](docs/building.md)
- [Derivatives](docs/differentiation.md)
- [Migrating from libneo](docs/migration_libneo.md)
- [Benchmarks](benchmark/README.md)

## Related projects

- [fortsym](https://github.com/lazy-fortran/fortsym): symbolic algebra used at
  build time to generate kernels under `src/generated`.
- [fortgen](https://github.com/lazy-fortran/fortgen): code-generation runtime
  (build dependency).
- [libneo](https://github.com/itpplasma/libneo): plasma physics library whose
  numerical routines are migrating to fortnum.

## Contributing

See [CONTRIBUTING.md](CONTRIBUTING.md).

## License

MIT. See [LICENSE](LICENSE).
