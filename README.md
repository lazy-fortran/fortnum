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
depends on parameter count. This is the many-parameter trajectory tournament:
complete primal-plus-gradient wall clock selects forward at one parameter and
reverse at four and sixteen.

A checkpointed Cash–Karp adjoint can retain every fourth or sixteenth state and
recompute intervening segments. On the measured 41-step trajectory it reduces
retained trace storage by 46–55%, but is 23–26% slower and does not reduce
end-to-end peak RSS because the current primal is compressed after integration.

A recomputation-only candidate retains just the accepted schedule and initial
state. It cuts retained trace storage by 59%, but repeated forward prefixes make
the complete VJP 5.04 times slower; full trace remains selected by wall clock.

Implicit ODE stages expose a batched `analytical` tangent solve for
`stage - base - alpha*rhs = 0`. Reusing one stage-Jacobian factorization makes
the measured four-state, sixteen-direction workload 4.88 times faster than
complete-solve finite differences, with equal measured peak RSS.

The matching contracted adjoint reuses one transposed factorization and sends
RHS cotangents through a parameter-VJP callback. With four states, sixteen
parameters, and one objective it is 20.21 times faster than complete-solve
finite differences; complete wall clock selects it in every measured regime.

Transversal event times use an `analytical` implicit JVP of the event residual,
not the root-location iterations. For four parameters and sixteen directions,
one primal integration plus all event-time products takes 2.48 µs versus
77.79 µs for complete-solve finite differences, a 31.37× wall-clock win.
The event-state interface composes the fixed-time tangent with
`f_event * dt_event`. For the complete `(event time, event state)` output, the
same four-parameter, sixteen-direction workload takes 2.71 µs analytically
versus 84.22 µs diagnostically, a 31.06× wall-clock win.

Small determinant products now expose `fortsym`-generated analytical JVPs.
For a complete 3x3 determinant value plus 64 directional products, analytical
takes 1.198 µs versus 2.208 µs for central finite differences, a 1.843×
wall-clock win with no derivative-workspace allocation.
The matching scalar-output VJP returns all nine 3x3 input sensitivities in
29.01 ns versus 159.54 ns for componentwise finite differences, a 5.50×
complete-workload wall-clock win.
For callers that require an explicit inverse, the fused 3x3 inverse-plus-JVP
uses the analytical identity `dAinv = -Ainv*dA*Ainv`: 125.08 ns versus
160.00 ns for a nominal inverse plus central differences, a 1.279× wall-clock
win.
The matching inverse-plus-VJP returns all nine input sensitivities in 113.43 ns
versus 929.10 ns for componentwise central differences, an 8.19× complete
wall-clock win.
A fixed-capacity `lu_factorization_t` now owns reusable LU factors and pivots.
For the measured 16×16 solve workload its type-bound solve takes 1.118 µs
versus 1.104 µs for the raw factored API, so raw reuse remains the hot-loop
winner while the object provides safer ownership.
Analytical linear-solve JVPs also accept multiple directions over one reused
factorization. On the 16×16 workload, repeated scalar products remain faster:
12.835 µs versus 22.584 µs for 16 directions, so batching is a convenience
interface rather than the selected hot-loop implementation.
The matching multiple-cotangent adjoint interface reuses one transposed
factorization. At 16 cotangents, repeated scalar VJPs take 18.925 µs versus
20.401 µs batched, so repeated scalar remains selected by wall clock.
A real forward-Enzyme comparator differentiates a 4×4 direct elimination
kernel. For 16 value-plus-JVP directions, analytical implicit differentiation
takes 0.774 µs, finite differences 0.982 µs, and `autodiff` 2.719 µs.
For the scalar-objective VJP with 20 active solve inputs, reverse Enzyme is much
more competitive: one cotangent takes 73.50 ns versus 59.40 ns analytically and
1,025.39 ns by componentwise finite differences.
A fixed-trace iterative-solver comparator distinguishes the executed iteration
map from the converged solve. At 32 iterations and 16 directions, the analytical
tangent recurrence takes 5.226 µs, finite differences 8.143 µs, and forward
`autodiff` 19.764 µs.
No hybrid BLAS/LAPACK custom rule is selected yet: the measured analytical
direct-solver JVP and VJP already beat Enzyme by 3.51× and 1.28× respectively,
so an extra external-rule boundary is not currently justified.
The completed linear-algebra tournament covers 12 product or interface
workloads. Analytical implementations win all 12 by complete-workload wall
clock; the closest mechanism contest is the direct-solver VJP, where analytical
is 1.28× faster than reverse Enzyme.

The ODE forward sensitivity now has an explicit continuous contract: it
approximates the variational IVP at fixed terminal time on the primal's frozen
mesh. Halving the maximum step reduced the closed-form-oracle error from
`1.41e-6` to `6.40e-11` over four refinements. For sixteen sequential
directions it is 1.10× faster than complete-solve finite differences.

The same recurrence has an explicit discrete contract: it is the exact tangent
of the accepted Cash-Karp steps with their schedule fixed. Independent
fixed-step replay gives maximum JVP and VJP errors of `7.61e-12` and
`1.25e-11`. For the full 16 by 16 initial-state Jacobian, sixteen analytical
forward sweeps take 187.99 µs and sixteen analytical reverse sweeps take
230.68 µs, including the primal trace. Forward is 1.23 times faster on this
square workload. Peak RSS and cache counters are supporting evidence;
complete-workload wall clock selects the implementation.

For the short, 41-step nonstiff trajectory, the requested product determines
the winner. One parameter JVP selects the analytical forward recurrence at
9.32 µs, while one scalar-objective VJP with two active initial-state inputs
selects the full-trace analytical reverse recurrence at 20.67 µs. The reverse
VJP is 1.74 times faster than reconstruction from two forward sweeps.
Checkpointing and recomputation reduce retained trace bytes but lose on
complete wall clock and do not reduce measured application peak RSS.

On the 400-step long nonstiff trajectory, one full-trace analytical reverse VJP
takes 207.65 µs including the primal. It is 1.71 times faster than two forward
sweeps, 1.24 times faster than checkpointed reverse, and 42.23 times faster
than recomputed reverse. Checkpointing and recomputation reduce retained trace
storage by 56.76% and 59.91%, but maximum observed peak RSS remains within
12 KB across candidates. Full-trace reverse therefore remains selected by
complete wall clock.

For the current explicit fallback on a stiff system with eigenvalues `-1` and
`-1000`, Cash-Karp requires 338 accepted steps to reach `t=1`. A full-trace
analytical reverse gradient takes 177.73 µs including the primal. It is 1.72
times faster than two analytical forward sweeps and 1.45 times faster than
finite differences. This measures the available differentiable path; fortnum
still needs a differentiable stiff primal before it can claim an
application-optimal stiff solution.

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
