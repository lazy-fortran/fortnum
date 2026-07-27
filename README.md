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
`fortnum` uses `fortsym` at build time for symbolic algebra and code
generation. The first pinned generator is under `tools/codegen/`; generated
production sources are committed under `src/generated/`.

Current derivative infrastructure includes generic analytical implicit JVP and
VJP boundaries for scalar and vector roots, analytical fixed-point and
linear-solve products, real autodiff/analytical hybrid Dawson and scalar-root
JVP/VJP boundaries, hybrid vector-root JVP/VJP candidates using Enzyme residual
products, reusable converged vector-root Jacobians and JVP/VJP factorizations,
opt-in vector-root JVP/VJP reliability reporting, and static benchmark-selected
candidates. The completed scalar-root tournament compares forward and reverse
autodiff through fixed Newton iterations against analytical implicit, hybrid,
and finite-difference diagnostic candidates. Measured runtime, dispersion,
peak memory, validation, hardware, and compiler evidence is committed in
`docs/design/differentiation_benchmarks.md` and `benchmark/reference/`.
`ROADMAP.md` is the authoritative implementation checklist.

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
- `tools/codegen/`: isolated `fortsym` build-time generators.

## License

MIT, Copyright (c) 2025 lazy-fortran.
