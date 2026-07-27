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
candidates. The completed scalar- and vector-root tournaments compare forward
and reverse autodiff through fixed Newton iterations against analytical
implicit, hybrid, and finite-difference diagnostic candidates. Measured
runtime, dispersion, peak memory, validation, hardware, and compiler evidence
is committed in
`docs/design/differentiation_benchmarks.md` and `benchmark/reference/`.
`ROADMAP.md` is the authoritative implementation checklist.

Fixed-bound parameterized integrals also expose an analytical JVP based on
differentiation under the integral sign, including the analytical Leibniz term
for an active lower or upper bound. Measured `hybrid` candidates use Enzyme
forward or reverse mode for integrand products and the corresponding analytical
fixed-quadrature contraction. With four active parameters and one scalar
output, one reverse VJP is about four times faster than four forward JVPs for
`analytical`, `autodiff`, and `hybrid` mechanisms; the completed fixed-rule
tournament selects `analytical` for both products. Adaptive integration has an
`analytical` frozen-subdivision JVP for QAG, QAGS, QAGP, and one-sided QAGIU
traces, a forward-`autodiff` candidate through a fixed accepted trace, and a
`hybrid` candidate using Enzyme integrand JVPs inside the analytical trace walk.
The smooth adaptive tournament selects a compact `analytical` frozen-trace
replay. For an integrable endpoint singularity, compact whole-trace `autodiff`
has the lowest median wall clock, but is within measurement noise of compact
`analytical`; both clearly beat the generic callback paths. Candidate selection
prioritizes complete-workload wall clock and records peak memory, input/output
scaling, and cache behavior where supported.

For batches of independent fixed integrals with four active inputs per scalar
output, the full-Jacobian tournament selects reverse `analytical`. At batch 16
it is 4.16 times faster than four forward analytical sweeps and scales nearly
linearly without material batch-sized memory growth.

ODE forward sensitivities can compose an Enzyme-generated RHS JVP with the
analytical frozen Cash–Karp trace. On the first scalar trajectory this `hybrid`
is 4.0% slower than an explicit RHS derivative but 1.54 times faster than
complete-solve finite differences.

Cash–Karp also has an analytical discrete adjoint over its frozen accepted-step
trace. For a two-state system and one scalar terminal objective, one reverse
adjoint is 1.74 times faster than reconstructing the VJP from two forward
tangent sweeps.

The discrete adjoint can accumulate RHS-parameter VJPs at every Runge–Kutta
stage. One forward sweep wins for one active parameter, while reverse wins by
2.54× at four parameters and 8.51× at sixteen; production selection therefore
depends on parameter count.

A checkpointed Cash–Karp adjoint can retain every fourth or sixteenth state and
recompute intervening segments. On the measured 41-step trajectory it reduces
retained trace storage by 46–55%, but is 23–26% slower and does not reduce
end-to-end peak RSS because the current primal is compressed after integration.

A recomputation-only candidate retains just the accepted schedule and initial
state. It cuts retained trace storage by 59%, but repeated forward prefixes make
the complete VJP 5.04 times slower; full trace remains selected by wall clock.

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
