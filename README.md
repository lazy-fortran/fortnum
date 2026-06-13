# fortnum

A primal-first, derivative-ready clean-room Fortran numerical library. fortnum
replaces GSL as the numerical backend for the itpplasma codes, providing special
functions, integration, quadrature, FFT, ODE solvers, root finding,
interpolation, and random number generation under a permissive MIT license.

"Primal-first" means every routine is written to compute its value correctly and
efficiently first; "derivative-ready" means the implementations are structured so
that automatic differentiation (planned Enzyme integration) can produce
derivatives without rewriting the numerics.

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

- `src/` — library modules, grouped by domain (`special/`, `integrate/`,
  `quadrature/`, `fft/`, `ode/`, `roots/`, `interp/`, `rng/`, `ad/`).
- `test/` — CTest suite, with `special/`, `oracle/`, and `ad/` subgroups.
- `benchmark/` — performance benchmarks.
- `docs/` — documentation, with `design/` for design notes.
- `cmake/` — CMake helper modules.

## License

MIT, Copyright (c) 2025 lazy-fortran.
