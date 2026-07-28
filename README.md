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
The current CPU/GPU support matrix and the strict device-leaf contract are in
`docs/design/gpu.md`. GPU differentiation is not yet a supported production
path: initial offload work is limited to validated generated `analytical`
leaves, and silent host fallback is a test failure.
GPU compilation is selected explicitly with
`-DFORTNUM_GPU_BACKEND=NONE|OPENACC|OPENMP`; `NONE` is the default, and an
unavailable selection fails configuration. OpenMP offload also requires
compiler-specific `FORTNUM_OPENMP_TARGET_FLAGS`, so host-only OpenMP cannot be
mistaken for GPU support.
The benchmark pilot uses one shared Dawson batch loop for OpenACC and OpenMP
target; both schedules call the same generated `analytical` leaf. The OpenACC
path now has proven NVIDIA execution and independent CPU validation. On the
reference RTX 5060 Ti its transfer-inclusive 1,048,576-element workload is
2.006 times faster than the CPU generated-leaf loop; this pilot is not a
general GPU-support claim.
The identical OpenMP-target pilot also proves non-host execution and is 1.6585
times faster than its CPU generated-leaf loop in the recorded run.
Keeping the Dawson batch resident reduces synchronous call time from 3.9251 ms
to 0.2098 ms with OpenACC and from 3.9800 ms to 0.2101 ms with OpenMP target.
The generated contracted analytical VJP also passes a batched adjoint identity
on both devices and takes 0.1690 ms per resident 1,048,576-element call.
Generated value/product fusion wins all 24 Dawson GPU comparisons across
launch-, throughput-, and transfer/memory-dominated batches. Resident fusion
is 1.45 to 1.96 times faster than two separate generated-kernel launches.
For the Dawson JVP, `fortsym` also proves raw, simplified, and factored
algebraic candidates equivalent before native benchmarking. The simplified
candidate reduces the post-CSE count from 30 to 10 operations and is 1.27
times faster than raw on the largest CPU workload and about 1.09 times faster
resident on either GPU. The simplified and factored GPU results are practical
ties, so CPU wall clock and the smaller 96-byte native symbol select the
simplified form. Only that generated production leaf and its symbolic
generator are committed; raw and factored losers are reproduced only in a
requested temporary tournament directory.
The generated multi-input scalar-output pilot covers 2, 4, 8, and 16 active
inputs and 1, 4, and 16 JVP directions or VJP cotangents. At 65,536 points,
resident GPU execution is 20.6 to 38.7 times faster than the pinned CPU,
depending on product and active-input count; at 256 points the CPU remains
faster. Both product costs scale linearly, and OpenACC and OpenMP target are
effectively tied while data is resident. Full wall-clock and memory evidence
is in `docs/design/gpu.md`.
The selected multi-input layout is `x(batch,active)`: it preserves coalesced
access across GPU threads, is never materially slower in resident workloads,
and beats `x(active,batch)` by 4.9 to 5.8 percent for the mid-size resident
VJP. Transfer-only differences were noisy and inconsistent, so they did not
override the resident and CPU evidence.
Released Nsight profiling reports 99.4 to 102.5 GB/s device-memory throughput,
24.3 to 24.6 percent achieved occupancy, 147 to 150 registers per thread, zero
local-memory spills, and about 153.1 MB peak device allocation for the
1,048,576-point resident workload. Compute throughput is about 84.7 percent
while DRAM throughput is about 23 percent, so this generated kernel is
compute-throughput limited rather than bandwidth limited.
The GPU evidence also has a pinned `fortplot` generator for batch and
active-input scaling, CPU/GPU crossover, achieved bandwidth, and peak device
allocation. Its normalized CSV inputs and Fortran source are committed under
`benchmark/report/`; generated PNGs remain ignored and are never committed.
For the measured eight-input workloads, OpenACC and OpenMP target are tied
within the three-percent timing threshold. A benchmark-level pre-launch lookup
therefore selects OpenACC using its 73,658-byte lower peak device allocation;
unmeasured workload keys return unavailable instead of guessing. Selection
occurs before choosing the separately compiled backend executable, never
inside the launch loop.
The first implicit GPU composition solves \(x^2-p=0\) with primal Newton
iterations, then combines `fortsym`-generated \(R\), \(R_x\), and \(R_p\,dp\)
with the shared analytical scalar-root boundary. It is `analytical`, not
`hybrid`: no autodiff participates, and solver iterations remain inactive.
Resident GPU execution is 17.0 to 21.8 times faster than CPU at 65,536 and
1,048,576 roots; CPU wins the 256-root launch-bound case.
GPU family coverage is tracked separately rather than as one oversized
milestone. Special functions are covered by the generated Dawson
value/JVP/VJP pilots, and residual kernels are covered by the generated
scalar-root residual plus analytical implicit JVP. Fixed-cell interpolation is
covered by generated Lagrange value/JVP/VJP leaves. At 1,048,576 points the
resident GPU is 15.0 times faster than CPU for JVP and 11.9 times faster for
VJP. With transfers included, GPU still wins JVP by 1.31 times, while CPU wins
VJP by 1.24 times because five adjoint arrays must return to the host. Fixed
quadrature is also covered: at 1,048,576 order-16 rules the best resident GPU
is 28.5 times faster for JVP and 26.7 times faster for VJP, but CPU wins when
transfers are included. Generated 3x3 determinant and inverse products cover
fixed-size linear algebra: resident GPU is 21.5 to 33.4 times faster at
1,048,576 matrices. Transfers keep determinant products on CPU, while inverse
products remain 1.32 to 1.36 times faster on GPU. A shared analytical
length-8 radix-2 leaf covers batched FFT JVP/VJP: resident GPU is 26.8 to 29.6
times faster at 1,048,576 transforms, while transfer-inclusive execution stays
on CPU. Fixed-trace two-state ODE tangent/adjoint recurrences are also covered:
at 1,048,576 trajectories and 64 recorded maps, resident GPU is 50.7 to 52.1
times faster and transfer-inclusive GPU is 17.7 to 17.9 times faster.
The corresponding terminal least-squares application composes forward state,
loss/cotangent, and reverse gradient inside one persistent data region. Its
complete value-plus-gradient workload is 50.6 times faster resident and 16.8
times faster with transfers at 1,048,576 trajectories; CPU remains selected
for the 256-trajectory transfer-inclusive workload.
`fortnum` uses `fortsym` at build time for symbolic algebra and code
generation. The first pinned generator is under `tools/codegen/`; generated
production sources are committed under `src/generated/`. The
`tools/codegen/fortsym.lock` file pins the generator dependency, and every
generated numerical kernel records that full revision in its banner. Generated
modules use pure procedures where valid, current `fortsym` line wrapping, and
`fo exec` regeneration commands. Run `tools/codegen/check_generated.sh` from
its directory to regenerate into a temporary directory and byte-compare all
committed numerical kernels. CTest and CI run the same check when
`FORTNUM_CHECK_GENERATED=ON`. Benchmark records for generated kernels also
capture compiler flags, native symbol sizes, structural operation counts, and
the exact `fortsym` revision. The Dawson generator now derives value-only,
contracted JVP/VJP, and fused value/product leaves from one symbolic DAG. On
the reference CPU the fused complete workload is 1.5699 times faster for JVP
and 1.6779 times faster for VJP than calling the generated value and product
leaves separately.
Generators are quiet by default to keep build logs concise; set
`FORTNUM_CODEGEN_VERBOSE=1` to print generated paths, equivalence proofs, and
operation counts.

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
The pre-refactor CPU Enzyme baseline covers all 17 fixtures: an Enzyme-only
clean build takes 3.11 s with 139.7 MB peak RSS, while sequential complete
fixture execution has a 27.5859 ms median and 3.6 MB peak RSS. The inventory
also records 13 duplicated peak-memory interfaces, nine median/MAD
implementations, and 22 raw Enzyme interfaces; these are the measured starting
point for the shared-scaffolding migration.
The migrated scalar fixtures now generate their mechanical Enzyme wrappers
with `fortsym`. Dawson uses the shared timing and custom-rule counter; its
normal-build analytical JVP changes by 1.0%, inside the 3% gate. The generated
square VJP wrapper is 6.5% faster than the equivalent raw wrapper and adds only
344 bytes to the smoke executable.
Scalar-root residual and fixed-Newton wrappers are now generated as well.
For the two-input scalar root, hybrid implicit JVP is selected at 1.891 ns
versus 2.297 ns analytically, while analytical implicit VJP is selected at
14.117 ns versus 19.490 ns for the reverse hybrid.
One internal support module now centralizes environment parsing, warmup/sample
collection, timing, median/MAD, standard output, and peak-RSS access. Its
independent plain-compiler and Flang/Enzyme tests pass; existing numerical
fixtures are intentionally unchanged until their dedicated migration items.
`fortsym` now also generates temporary Enzyme JVP/VJP wrappers for one to four
active scalar inputs and a scalar result, with an optional analytical forward
rule. All eight generated products pass real-Enzyme formula and adjoint tests;
the wrapper-plus-polynomial-kernel medians range from 7.21 to 8.30 ns.
A shared 52-byte rule counter now proves that Enzyme selected a generated
analytical forward rule without per-kernel counter state. Enabled and disabled
timings are indistinguishable within combined dispersion, and production
benchmarks disable the counter.
The cumulative report in `docs/design/differentiation_report.md` summarizes 32
mechanism tournaments and provides reproducible `fortplot` figure generation
without committing generated PNGs.
`ROADMAP.md` is the authoritative implementation checklist.

Generated-source ownership is explicit. Symbolic specifications, generators,
selected numerical kernels, provenance, and benchmark evidence are committed.
Temporary algebraic candidates and routine Enzyme wrappers are generated in
the build tree. Stable recurrences, solver orchestration, traces, and
independent validation oracles remain hand-written.

Lagrange interpolation exposes analytical JVP and VJP products for active
support-node locations as well as active evaluation points and sampled values.
At 16 nodes the active-node analytical JVP is 2.31 times faster and the VJP is
16.61 times faster than complete finite-difference diagnostics.

Within a fixed knot span, B-spline values also expose analytical products for
active breakpoint locations. At 18 breakpoints analytical is 2.04 times faster
for JVP and 2.09 times faster for VJP than complete finite differences.

Fitted B-spline coefficients defined by a collocation solve expose
factorization-reusing analytical implicit JVP and VJP products. At 16
coefficients they are respectively 4.90 and 1,224 times faster than complete
finite-difference diagnostics.

Interpolation callers can explicitly guard a directional probe with
`grid_search_derivative_status`; a changed cell reports non-smoothness through
`FORTNUM_DOMAIN_ERROR` without adding overhead to value-only searches.
B-spline callers have the matching `bspline_span_derivative_status` guard for
knot-span crossings.

Separate and fused combined-active Lagrange products are performance ties over
4--16 nodes; both APIs remain available, and callers use the fused form when
simultaneous evaluation-point and value activity is natural.
For B-splines, the fused combined-active implementation shares one basis pass
and is selected: at 16 coefficients it is 1.54 times faster for JVP and 1.45
times faster for VJP.
Inside a fixed cubic B-spline span, analytical and Enzyme products are both
validated. After generated-wrapper and shared-fixture migration, analytical is
20.1% faster for 16 JVPs and 10.2% faster for 16 VJPs. Removing string dispatch
from the measured inner loop reduces selected complete wall clock by 66.4% and
60.4%, respectively.

The modified-Bessel tournament likewise selects by region and product. For 16
products, analytical wins series JVP/VJP and asymptotic VJP; raw Enzyme wins
recurrence JVP and asymptotic JVP. After migration to shared generated
scaffolding, complete wall clock selects analytical for recurrence VJP by
8.0% over reverse autodiff. The measured hybrid JVP proves the
analytical custom-rule boundary but does not win because it evaluates both
`I0` and `I1` separately.

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
After generated-wrapper migration, the smooth Enzyme-enabled run is a close
three-way race: compact hybrid takes 1.260 µs, whole-trace autodiff 1.264 µs,
and compact analytical 1.269 µs. Analytical remains the normal-build
selection; the singular Enzyme-enabled workload selects autodiff at 2.807 µs
versus 2.836 µs analytical. The adaptive trace and its guards are not generated.
ODE forward sensitivity now generates only the local three-input Enzyme RHS
wrapper. The complete trajectory-plus-JVP workload takes 9.071 µs hybrid,
9.169 µs analytical, and 14.846 µs by complete-solve finite differences.
Analytical remains the normal-build selection; the solver and variational
equation stay explicit.
The fixed-quadrature Enzyme mechanics are generated from one five-input
integrand profile and one four-input whole-operator profile; the old hand-written
forward and reverse wrapper modules are gone. For a full four-input gradient,
four analytical JVPs take 1.223 µs while one analytical VJP takes 0.301 µs.
Whole-operator autodiff takes 1.702 µs and 0.440 µs, respectively.

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
A generated-wrapper forward-Enzyme comparator differentiates a 4×4 direct
elimination kernel. For 16 value-plus-JVP directions, analytical implicit
differentiation takes 0.643 µs, finite differences 0.816 µs, and `autodiff`
2.449 µs. For the scalar-objective VJP with 20 active solve inputs, the fused
value-plus-reverse wrapper takes 0.974 µs for 16 cotangents versus 0.637 µs
analytically and 14.037 µs by componentwise finite differences.
A fixed-trace iterative-solver comparator distinguishes the executed iteration
map from the converged solve. At 32 iterations and 16 directions, the analytical
tangent recurrence takes 4.846 µs, finite differences 7.628 µs, and forward
`autodiff` 18.092 µs after generated-wrapper migration.
The coupled two-state vector root uses seven generated Enzyme wrappers and the
shared fixture support. Complete wall clock selects the hybrid implicit JVP at
78.27 ns, effectively tied with analytical at 78.31 ns, and the analytical
implicit VJP at 78.80 ns over hybrid at 79.82 ns. Autodiff through 12 Newton
iterations takes 218.70 ns for JVP and 206.30 ns for VJP; complete-solve
finite-difference diagnostics take 341.65 ns and 343.62 ns.
No hybrid BLAS/LAPACK custom rule is selected yet: the measured analytical
direct-solver JVP and VJP already beat Enzyme by 3.81× and 1.53× respectively,
so an extra external-rule boundary is not currently justified.
The completed linear-algebra tournament covers 12 product or interface
workloads. Analytical implementations win all 12 by complete-workload wall
clock; the closest mechanism contest is the direct-solver VJP, where analytical
is 1.53× faster than reverse Enzyme.

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

For an accelerator build, add exactly one explicit backend selection. For
example, an OpenACC-capable compiler uses
`-DFORTNUM_GPU_BACKEND=OPENACC`. OpenMP target builds additionally pass
`-DFORTNUM_OPENMP_TARGET_FLAGS="<compiler target flags>"`.

An `fpm.toml` is provided for native `fo` builds and fpm-based consumers:

```
fpm build
fo
```

The manifest's `[extra.fo.test-args]` table supplies each oracle test with its
own reference-data arguments. Plain `fpm test` cannot express different
arguments per automatically discovered test; use bare `fo` for the complete
native build, 92-test suite, static checks, and lint pipeline.

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
