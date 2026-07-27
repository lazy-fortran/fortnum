# fortnum

A primal-first, derivative-plural clean-room Fortran numerical library. fortnum
is the numerical backend for the itpplasma codes, providing special
functions, integration, quadrature, FFT, ODE solvers, root finding,
interpolation, and random number generation under a permissive MIT license.

"Primal-first" means every routine computes its value correctly and efficiently
before derivative concerns. "Derivative-plural" means a product can have
autodiff, analytical, and hybrid implementations. Independent validation plus
measured application runtime and peak memory select the production candidate.

See `docs/performance_optimal_differentiation.md` for the rationale and
`docs/design/differentiation_plan.md` for the implementation plan.
`fortnum` will use the unfinished `../fortsym` library for future symbolic
algebra and code generation. No integration API is fixed yet.

## Build

CMake is the primary build system with CTest integration:

```
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Debug
cmake --build build -j
ctest --test-dir build --output-on-failure
```

An `fpm.toml` is provided so the `fo` tool and `fpm` work as well:

```
fpm build
fpm test
```

## Layout

- `src/`: library modules, grouped by domain (`special/`, `integrate/`,
  `quadrature/`, `fft/`, `ode/`, `roots/`, `interp/`, `rng/`, `ad/`).
- `test/`: CTest suite, with `special/`, `oracle/`, and `ad/` subgroups.
- `benchmark/`: performance benchmarks.
- `docs/`: documentation, with `design/` for design notes.
- `cmake/`: CMake helper modules.

## License

MIT, Copyright (c) 2025 lazy-fortran.
