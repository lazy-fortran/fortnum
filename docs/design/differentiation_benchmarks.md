# Differentiation benchmark evidence

Committed tables are reference measurements, not portable promises. Production
selection must be regenerated after a material change in hardware, compiler,
Enzyme, `fortsym`, primal source, or workload.

Every speedup or slowdown in this document compares two rows in the same table
for the same derivative operation and workload. A ratio stated as “B is
`r` times slower than A” means `runtime(B)/runtime(A) = r`; “B is `r` times
faster than A” means `runtime(A)/runtime(B) = r`. No comparison is against the
primal operation unless the table explicitly names a primal row.

## Dawson outer-expression JVP

The workload is

```text
f = dawson(x)
y = sin(f) + f**2
dy = J_y(x) v
```

It compares:

- `analytical`: a fused chain-rule expression
- `autodiff`: Enzyme differentiates the Dawson implementation and outer
  expression
- `hybrid`: Enzyme differentiates the outer expression and calls the analytical
  Dawson JVP at the operator boundary

The independent oracle is central finite difference of the complete outer
function. The hybrid test also counts analytical-rule calls, so agreement alone
cannot hide a failure to apply the custom rule.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, Flang/Clang/LLVM 22.1.8, Enzyme
`c965083`, 15 samples of 200,000 calls after three warmups.

| Candidate | Mechanism | Median ns/call | MAD ns/call | Relative to analytical |
|---|---|---:|---:|---:|
| analytical | `analytical` | 21.1663 | 0.1903 | 1.0000 |
| autodiff | `autodiff` | 23.3390 | 0.1078 | 1.1026 |
| hybrid | `hybrid` | 20.7559 | 0.1828 | 0.9806 |

Peak RSS for the shared tournament process was 2,764,800 bytes. These scalar
candidates allocate no candidate-specific dynamic storage, so this measurement
does not distinguish their memory use. `hybrid` is the raw latency winner and
is 11.07% faster than `autodiff`. It is only 1.94% faster than `analytical`,
inside the registry's 3% timing guard. With no measured memory distinction,
the deterministic stable-ID tie break selects `analytical` for this workload.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_dawson_outer.json`.

Run the validation and tournament with:

```bash
ctest --test-dir build-enzyme -R enzyme_dawson_hybrid --output-on-failure
taskset -c 4 \
  ./build-enzyme/test/ad/enzyme_dawson_hybrid.enzyme/enzyme_dawson_hybrid \
  --benchmark
```

Measure process peak memory separately so compilation is not included:

```bash
taskset -c 4 /usr/bin/time -v \
  ./build-enzyme/test/ad/enzyme_dawson_hybrid.enzyme/enzyme_dawson_hybrid \
  --benchmark
```

## Analytical linear-solve JVP factorization reuse

This workload applies 200,000 JVP directions to one dense 16-by-16 primal
system. It compares refactorizing the primal matrix for every analytical JVP
with reusing one compact partial-pivoted LU factorization. Two directions are
validated independently by central finite differences of the complete solve.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, 15 samples
after three warmups.

| Candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| refactor each JVP | `analytical` | 981.8472 | 7.2850 | 12,451,840 B |
| reuse primal LU | `analytical` | 348.6946 | 0.5354 | 12,378,112 B |

Reusing the factorization is 2.8158 times faster and its maximum observed
candidate-process RSS is 73,728 bytes lower. Memory is measured inside five
separately launched processes per candidate with `getrusage(RUSAGE_SELF)`;
therefore these figures exclude the `fo` runner. The counter is independently
tested by allocating and touching 16 MiB and requiring the reported high-water
mark to rise by at least 8 MiB.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_linear_solve_jvp_reuse.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_linear_solve_jvp refactor
taskset -c 4 fo exec bench_linear_solve_jvp reuse
fo exec bench_linear_solve_jvp refactor --peak-rss
fo exec bench_linear_solve_jvp reuse --peak-rss
```

## Analytical linear-solve VJP transpose-factorization reuse

This workload applies 200,000 cotangents to one dense 16-by-16 primal system.
It compares refactorizing the transposed primal matrix for every analytical VJP
with reusing one compact partial-pivoted LU factorization. Two cotangents are
validated independently by central finite differences of the contracted
complete solve.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, 15 samples
after three warmups.

| Candidate | Mechanism | Median ns/VJP | MAD ns/VJP | Peak RSS |
|---|---|---:|---:|---:|
| refactor transpose each VJP | `analytical` | 1,001.4936 | 2.9514 | 12,513,280 B |
| reuse transposed primal LU | `analytical` | 341.5936 | 2.1059 | 12,578,816 B |

Here “2.9318 times faster” compares the reuse row with the refactor row:
`1001.4936 / 341.5936 = 2.9318`. Maximum observed candidate-process RSS is
65,536 bytes higher for reuse. The target audit confirms this executable is
built from `bench_linear_solve_vjp.f90`; its derivative oracle is the
cotangent-contracted finite difference of the complete solve.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_linear_solve_vjp_reuse.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_linear_solve_vjp refactor
taskset -c 4 fo exec bench_linear_solve_vjp reuse
fo exec bench_linear_solve_vjp refactor --peak-rss
fo exec bench_linear_solve_vjp reuse --peak-rss
```

## Generic analytical implicit tangent boundary for scalar roots

The generic boundary accepts a converged scalar root and a callback that
returns the contracted residual products `f_x` and `f_p*tp`. It then applies
the analytical implicit rule `dx = -(f_p*tp)/f_x`; root-solver iterations are
not differentiated. An independent central-difference oracle re-solves
`x^3 + p1*x - p2 = 0` at `p + h*tp` and `p - h*tp`.

The benchmark applies 2,000,000 directions per sample to that same residual.
Both rows compute the same scalar-root JVP. The baseline calls the existing
low-level `root_jvp` with an assembled two-component residual gradient. The
candidate calls the new contracted callback boundary.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| assembled residual gradient | `analytical` | 23.2538 | 0.0708 | 12,271,616 B |
| contracted callback boundary | `analytical` | 27.9348 | 0.1342 | 12,279,808 B |

Here “1.2013 times slower” means
`27.9348 / 23.2538 = 1.2013`: it compares the callback-boundary row with the
assembled-gradient row for the same JVP, not with a primal root solve. The
assembled path remains the isolated-kernel selection. The generic boundary
costs 4.6810 ns/JVP while avoiding full residual-gradient materialization and
providing the operator interface needed for later `hybrid` residual products.
Its maximum observed candidate-process RSS is 8,192 bytes higher.

Peak memory is the maximum across five separately launched processes per
candidate, measured with `getrusage(RUSAGE_SELF)`. The machine-readable record
is
`benchmark/reference/ryzen9_5950x_scalar_root_tangent_boundary.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_scalar_root_tangent assembled
taskset -c 4 fo exec bench_scalar_root_tangent boundary
fo exec bench_scalar_root_tangent assembled --peak-rss
fo exec bench_scalar_root_tangent boundary --peak-rss
```

## Hybrid scalar-root JVP with Enzyme residual products

Both candidates use the same generic analytical implicit boundary and compute
the same scalar-root JVP. The `analytical` callback evaluates `f_x` and
`f_p*tp` explicitly. The `hybrid` callback invokes Enzyme forward mode twice
on the scalar residual `f = x^3 + p1*x - p2`: once for `f_x`, and once for the
parameter direction. Neither candidate differentiates root-solver iterations.

An independent oracle completely bisects the root at `p + h*tp` and
`p - h*tp`. Successful LLVM linking also requires Enzyme to eliminate the raw
`__enzyme_fwddiff` intrinsics, so the test cannot silently fall back to the
analytical callback.

The benchmark applies 2,000,000 directions per sample while varying the root
and preserving the residual equation. This is a residual-product comparison;
autodiff through the iteration trace remains a later tournament candidate.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, Flang/LLVM 22.1.8, Enzyme
`c96508349d9f`, Release `-O2`, 15 samples after three warmups.

| Candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| explicit residual products | `analytical` | 6.7477 | 0.0879 | 2,957,312 B |
| Enzyme forward residual products | `hybrid` | 5.1970 | 0.2172 | 2,928,640 B |

Here “1.2984 times faster” means
`6.7477 / 5.1970 = 1.2984`: it compares the `hybrid` row with the
`analytical` row for the same implicit JVP, not with a primal root solve.
`hybrid` saves 1.5507 ns/JVP and is the measured selection. Its maximum
observed self-process RSS is 28,672 bytes lower.

Peak memory is the maximum across five separately launched candidates,
measured with `getrusage(RUSAGE_SELF)`. The machine-readable record is
`benchmark/reference/ryzen9_5950x_scalar_root_hybrid_jvp.json`.

Run validation and timing with:

```bash
cmake --build build-enzyme --target enzyme_scalar_root_hybrid_build
ctest --test-dir build-enzyme -R '^enzyme_scalar_root_hybrid$' \
    --output-on-failure
taskset -c 4 \
    build-enzyme/test/ad/enzyme_scalar_root_hybrid.enzyme/enzyme_scalar_root_hybrid \
    --benchmark
```

## Hybrid scalar-root VJP

This comparison computes the same two-component VJP of the converged root for
the cubic residual `f(x,p) = x^3 + p1*x - p2`. Both candidates use the
analytical implicit adjoint boundary. The `analytical` candidate supplies
explicit residual partials; the `hybrid` candidate obtains `f_x`, `f_p1`, and
`f_p2` in one Enzyme reverse sweep and contracts the parameter partials with
the supplied root cotangent. Neither candidate differentiates root-solver
iterations.

The independent oracle centrally differences the complete scalar objective
`L(p) = 1.3*x*(p)`. Every perturbed `x*(p)` is obtained with a fresh 100-step
bisection solve. The maximum absolute error was `1.1891e-11`.

The benchmark applies 2,000,000 cotangents per sample while varying the root
and preserving the residual equation. This is a residual-product comparison;
autodiff through the iteration trace remains a later tournament candidate.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, Flang/LLVM 22.1.8, Enzyme
`c96508349d9f`, Release `-O2`, 15 samples after three warmups.

| Candidate | Mechanism | Median ns/VJP | MAD ns/VJP | Peak RSS |
|---|---|---:|---:|---:|
| explicit residual products | `analytical` | 27.3548 | 0.2323 | 2,973,696 B |
| Enzyme reverse residual products | `hybrid` | 16.3433 | 0.1202 | 2,965,504 B |

Here “1.6738 times faster” means
`27.3548 / 16.3433 = 1.6738`: it compares the `hybrid` row with the
`analytical` row for the same implicit VJP, not with a primal root solve.
`hybrid` saves 11.0115 ns/VJP and is the measured selection. Its maximum
observed self-process RSS is 8,192 bytes lower.

Peak memory is the maximum across five separately launched candidates,
measured with `getrusage(RUSAGE_SELF)`. The machine-readable record is
`benchmark/reference/ryzen9_5950x_scalar_root_hybrid_vjp.json`.

Run validation and timing with:

```bash
cmake --build build-enzyme --target enzyme_scalar_root_vjp_hybrid_build
ctest --test-dir build-enzyme -R '^enzyme_scalar_root_vjp_hybrid$' \
    --output-on-failure
taskset -c 4 \
    build-enzyme/test/ad/enzyme_scalar_root_vjp_hybrid.enzyme/enzyme_scalar_root_vjp_hybrid \
    --benchmark
```

## Hybrid fixed-quadrature integrand JVP

This workload applies a reusable 32-point Gauss-Legendre rule on `[0,1]` to

```text
f(x,p) = exp(p1*x) + sin(p2*x) + p3*x^2 + p4*x^3.
```

The `analytical` candidate evaluates the explicit integrand directional
derivative and calls `gauss_legendre_jvp`. The `hybrid` candidate obtains each
integrand directional derivative with Enzyme forward mode, then calls the same
analytical quadrature contraction. The diagnostic centrally differences two
complete fixed-quadrature evaluations. Rule construction is outside all three
timed paths because the fixed nodes and weights are reusable.

An independently integrated closed form validates four parameter directions.
The `analytical` and `hybrid` errors must be below `2e-14`; the diagnostic error
must be below `2e-10`.

The timed workload requests one, two, or four JVP directions for one scalar
integral output. Each row is a fresh process pinned to CPU 4, with 31 samples
of 10,000 complete workloads after three warmups. Reference: AMD Ryzen 9
5950X, Flang/LLVM 22.1.8, Enzyme `c96508349d9f`, Release `-O2`.

| Directions | Candidate | Mechanism | Median ns/workload | MAD ns | Peak RSS |
|---:|---|---|---:|---:|---:|
| 1 | explicit integrand JVP + fixed quadrature | `analytical` | 312.6325 | 4.0857 | 3,035,136 B |
| 1 | Enzyme integrand JVP + fixed quadrature | `hybrid` | 313.4861 | 1.0931 | 3,047,424 B |
| 1 | Enzyme through complete fixed quadrature | `autodiff` | 429.6297 | 3.4785 | 2,797,568 B |
| 1 | two complete quadrature evaluations | diagnostic | 572.1548 | 2.9014 | 3,051,520 B |
| 2 | explicit integrand JVP + fixed quadrature | `analytical` | 629.2495 | 3.1329 | 2,797,568 B |
| 2 | Enzyme integrand JVP + fixed quadrature | `hybrid` | 633.8662 | 2.8634 | 2,838,528 B |
| 2 | Enzyme through complete fixed quadrature | `autodiff` | 867.5331 | 4.9113 | 3,010,560 B |
| 2 | two complete quadrature evaluations per direction | diagnostic | 1,148.1500 | 4.9874 | 2,789,376 B |
| 4 | explicit integrand JVP + fixed quadrature | `analytical` | 1,247.7167 | 5.8860 | 2,797,568 B |
| 4 | Enzyme integrand JVP + fixed quadrature | `hybrid` | 1,263.4433 | 10.5408 | 2,797,568 B |
| 4 | Enzyme through complete fixed quadrature | `autodiff` | 1,718.8515 | 7.7717 | 3,043,328 B |
| 4 | two complete quadrature evaluations per direction | diagnostic | 2,316.5791 | 20.7341 | 2,801,664 B |

Complete-workload wall clock selects `analytical` for this integrand, but only
by 1.0027, 1.0073, and 1.0126 times at one, two, and four directions. The
`hybrid` path remains 1.8113 to 1.8335 times faster than the diagnostic.
Whole-operator `autodiff` is 1.3742 times slower than `analytical` at one
direction, but 1.3317 times faster than the diagnostic. All candidates scale
approximately linearly because each requested forward direction repeats the
complete JVP: the four/one ratios are 3.9910 (`analytical`), 4.0303
(`hybrid`), 4.0008 (`autodiff`), and 4.0489 (diagnostic).

Linux `perf stat -r 5` for 20,000 four-direction workloads gives:

| Candidate | Cycles | Instructions | Cache references | Cache misses | Miss rate |
|---|---:|---:|---:|---:|---:|
| `analytical` | 113,267,878 | 401,594,792 | 122,880 | 16,090 | 13.0941% |
| `hybrid` | 114,276,069 | 405,914,614 | 89,225 | 15,535 | 17.4110% |
| `autodiff` | 156,936,769 | 559,455,522 | 93,053 | 15,705 | 16.8775% |
| diagnostic | 210,217,773 | 727,036,238 | 174,804 | 18,374 | 10.5112% |

The cycle and instruction counts corroborate the wall-clock verdict.
Hardware cache-reference counts varied by 32% to 40% across repetitions on
this host, so they are diagnostic rather than selection-grade evidence;
absolute cache misses were close and do not change the winner. Peak RSS also
does not provide a meaningful discriminator at this small fixed rule size.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_fixed_quadrature_hybrid_jvp.json`.

Run validation and isolated candidate timing with:

```bash
cmake --build build-enzyme \
    --target enzyme_fixed_quadrature_jvp_hybrid_build
ctest --test-dir build-enzyme \
    -R '^enzyme_fixed_quadrature_jvp_hybrid$' --output-on-failure
taskset -c 4 \
    build-enzyme/test/ad/enzyme_fixed_quadrature_jvp_hybrid.enzyme/enzyme_fixed_quadrature_jvp_hybrid \
    --benchmark analytical 4
taskset -c 4 \
    build-enzyme/test/ad/enzyme_fixed_quadrature_jvp_hybrid.enzyme/enzyme_fixed_quadrature_jvp_hybrid \
    --benchmark hybrid 4
```

## Hybrid fixed-quadrature integrand VJP

This reverse-mode experiment uses the same 32-point rule, four-parameter
integrand, and independently integrated closed-form oracle as the preceding
JVP experiment. The analytical quadrature VJP first maps scalar output
cotangent `u` to node cotangents `u*w(i)`. The `analytical` candidate applies
explicit integrand partials at each node; the `hybrid` candidate uses one
Enzyme reverse sweep per node and accumulates all four parameter cotangents.
The diagnostic uses eight complete quadrature evaluations to centrally
difference all four parameters.

The timed workload requests one, two, or four scalar-output cotangents. Each
VJP returns all four parameter sensitivities. Rows are separate processes
pinned to CPU 4, with 31 samples of 10,000 complete workloads after three
warmups and a cooldown between sustained runs. Reference: AMD Ryzen 9 5950X,
Flang/LLVM 22.1.8, Enzyme `c96508349d9f`, Release `-O2`.

| Cotangents | Candidate | Mechanism | Median ns/workload | MAD ns | Peak RSS |
|---:|---|---|---:|---:|---:|
| 1 | explicit integrand VJP + fixed quadrature | `analytical` | 299.4435 | 1.1752 | 3,051,520 B |
| 1 | Enzyme integrand VJP + fixed quadrature | `hybrid` | 317.5576 | 5.2910 | 2,797,568 B |
| 1 | Enzyme through complete fixed quadrature | `autodiff` | 443.6331 | 3.1199 | 2,785,280 B |
| 1 | eight complete quadrature evaluations | diagnostic | 2,263.7966 | 12.9915 | 2,834,432 B |
| 2 | explicit integrand VJP + fixed quadrature | `analytical` | 604.4575 | 2.9396 | 2,789,376 B |
| 2 | Enzyme integrand VJP + fixed quadrature | `hybrid` | 646.4897 | 4.8511 | 2,813,952 B |
| 2 | Enzyme through complete fixed quadrature | `autodiff` | 897.8271 | 4.7469 | 2,969,600 B |
| 2 | eight complete evaluations per cotangent | diagnostic | 4,547.2562 | 45.2453 | 2,813,952 B |
| 4 | explicit integrand VJP + fixed quadrature | `analytical` | 1,213.0276 | 2.8915 | 2,793,472 B |
| 4 | Enzyme integrand VJP + fixed quadrature | `hybrid` | 1,275.7279 | 4.1718 | 2,826,240 B |
| 4 | Enzyme through complete fixed quadrature | `autodiff` | 1,804.7173 | 30.2730 | 3,026,944 B |
| 4 | eight complete evaluations per cotangent | diagnostic | 9,187.4429 | 190.7125 | 2,830,336 B |

Complete-workload wall clock selects `analytical`; it is 1.0517 to 1.0695
times faster than `hybrid` and 1.4815 times faster than whole-operator
`autodiff` for one cotangent. `Autodiff` is still 5.1029 times faster than
finite differences. Batching independent cotangents scales linearly:
four/one ratios are 4.0509 (`analytical`), 4.0173 (`hybrid`), 4.0680
(`autodiff`), and 4.0584 (diagnostic). Peak RSS does not select a candidate.

The important forward/reverse result uses complete derivative workloads, not
isolated scalar partials:

| Required derivative | `analytical` | `hybrid` | `autodiff` | Diagnostic |
|---|---:|---:|---:|---:|
| full gradient via four forward JVPs | 1,247.7167 ns | 1,263.4433 ns | 1,718.8515 ns | 2,316.5791 ns |
| full gradient via one reverse VJP | 299.4435 ns | 317.5576 ns | 443.6331 ns | 2,263.7966 ns |
| forward/reverse ratio | 4.1668 | 3.9786 | 3.8745 | 1.0233 |

For this one-output, four-input map, reverse mode wins by approximately four
times because a single VJP returns all input sensitivities, whereas forward
mode needs one JVP per input basis direction. This is a workload-shape verdict,
not a universal reverse-mode rule; multi-output crossover measurements remain
in the roadmap.

Linux `perf stat -r 5` for 20,000 single-cotangent workloads gives:

| Candidate | Cycles | Instructions | Cache references | Cache misses | Miss rate |
|---|---:|---:|---:|---:|---:|
| `analytical` | 28,978,599 | 101,314,668 | 248,437 | 16,004 | 6.4419% |
| `hybrid` | 29,425,584 | 103,934,596 | 75,711 | 16,438 | 21.7115% |
| `autodiff` | 41,127,520 | 146,615,438 | 54,478 | 15,764 | 28.9365% |
| diagnostic | 211,788,713 | 705,696,347 | 510,342 | 27,868 | 5.4607% |

Cycles and instructions corroborate the wall-clock ranking. Cache-reference
counts varied by 12% to 76% across repetitions and are therefore diagnostic,
not selection-grade; wall clock remains decisive.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_fixed_quadrature_hybrid_vjp.json`.

Run validation and isolated timing with:

```bash
cmake --build build-enzyme \
    --target enzyme_fixed_quadrature_vjp_hybrid_build
ctest --test-dir build-enzyme \
    -R '^enzyme_fixed_quadrature_vjp_hybrid$' --output-on-failure
taskset -c 4 \
    build-enzyme/test/ad/enzyme_fixed_quadrature_vjp_hybrid.enzyme/enzyme_fixed_quadrature_vjp_hybrid \
    --benchmark analytical 1
taskset -c 4 \
    build-enzyme/test/ad/enzyme_fixed_quadrature_vjp_hybrid.enzyme/enzyme_fixed_quadrature_vjp_hybrid \
    --benchmark hybrid 1
```

## Analytical frozen-trace adaptive integration

`integrate_qag_jvp` replays the Gauss-Kronrod rule on the exact subdivision
accepted by the primal and applies an explicit tangent integrand. The same
analytical trace walk underlies the QAGS, QAGP, and supported one-sided QAGIU
products. The existing behavioral test uses

```text
I(p) = integral_0^1 exp(p*x) dx,  p = 12,
```

where the low-order rule accepts three panels. It checks the JVP against the
closed-form derivative and complete adaptive central differences, and verifies
that the perturbations retain the same subdivision. Failed and nonsmooth
primal traces are separately checked to reject the derivative.

The performance workload includes everything required for a value and
derivative call: it constructs the adaptive primal trace, then either performs
one analytical tangent trace walk or two primal trace walks at `p+h` and
`p-h`. Thus the table reports combined primal-plus-derivative wall clock, not
only isolated trace replay. It evaluates 10,000 varying parameters per sample,
with 15 samples after three warmups. Reference: AMD Ryzen 9 5950X, CPU 4
pinned, GNU Fortran 16.1.1, Release build.

| Candidate | Mechanism | Median ns/value+JVP | MAD ns | Peak RSS |
|---|---|---:|---:|---:|
| adaptive primal + explicit tangent trace walk | `analytical` | 1,113.8565 | 4.4133 | 2,912,256 B |
| adaptive primal + two frozen primal trace walks | diagnostic | 1,596.3970 | 6.2468 | 2,596,864 B |

Complete-workload wall clock selects `analytical`, which is 1.4332 times
faster and saves 482.5405 ns per value-plus-derivative call. The small peak-RSS
difference is process-layout noise: both candidates reuse the same bounded
workspace and allocate no storage proportional to the number of derivative
directions.

Linux `perf stat -r 5` over the same 10,000-workload process gives:

| Candidate | Cycles | Instructions | Cache references | Cache misses | Miss rate |
|---|---:|---:|---:|---:|---:|
| `analytical` | 80,269,951 | 134,644,416 | 713,244 | 19,150 | 2.6849% |
| diagnostic | 97,032,487 | 173,874,284 | 944,545 | 19,272 | 2.0403% |

The diagnostic uses 1.2088 times as many cycles, 1.2914 times as many
instructions, and 1.3243 times as many cache references. Absolute cache misses
are nearly equal, so wall clock remains the selection criterion.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_integrate_frozen_jvp.json`.

Run the benchmark with:

```bash
cmake -S benchmark -B benchmark/build -G Ninja \
    -DCMAKE_BUILD_TYPE=Release
cmake --build benchmark/build --target bench_integrate_frozen_jvp
taskset -c 4 benchmark/build/bin/bench_integrate_frozen_jvp analytical
taskset -c 4 benchmark/build/bin/bench_integrate_frozen_jvp diagnostic
```

## Forward autodiff through a frozen adaptive trace

The first `autodiff` frozen-trace candidate records the actual three-panel
QAG trace for the preceding exponential workload, then evaluates its fixed
G7K15 weighted sum in an Enzyme-safe scalar kernel. Enzyme forward mode
differentiates that complete frozen trace. The fixture verifies that the
adaptive primal still selected `[0,0.5]`, `[0.5,0.75]`, and `[0.75,1]`;
otherwise it rejects the derivative rather than silently using a stale trace.

All timed candidates include adaptive primal trace construction:

- `analytical` replays the public generic trace JVP with an explicit tangent
  integrand;
- `autodiff` applies Enzyme forward mode to the fixed weighted-sum kernel;
- the diagnostic centrally differences two evaluations of that same frozen
  weighted sum.

The independent oracle is the closed-form derivative of the exactly integrated
exponential. The workload evaluates 10,000 varying parameters per sample, with
15 samples after three warmups. Reference: AMD Ryzen 9 5950X, CPU 4 pinned,
Flang/LLVM 22.1.8, Enzyme `c96508349d9f`, Release `-O2`.

| Candidate | Mechanism | Median ns/value+JVP | MAD ns | Peak RSS |
|---|---|---:|---:|---:|
| generic explicit-tangent trace replay | `analytical` | 2,564.1478 | 39.2610 | 2,715,648 B |
| Enzyme through fixed weighted sum | `autodiff` | 1,306.9523 | 2.4366 | 2,547,712 B |
| frozen weighted-sum central difference | diagnostic | 1,423.0742 | 6.0043 | 2,580,480 B |

Complete-workload wall clock selects `autodiff`: it is 1.9619 times faster
than the generic `analytical` trace replay and 1.0889 times faster than the
diagnostic. The mechanism name is not the explanation. The generic analytical
path calls the full Gauss-Kronrod panel evaluator, including error-estimate
bookkeeping that is unnecessary once the trace is fixed; Enzyme differentiates
a compact fixed weighted sum. A later tournament may add an equally compact
analytical replay rather than treating this workload-specific result as a
universal autodiff verdict.

Linux `perf stat -r 5` over the same 10,000-workload process gives:

| Candidate | Cycles | Instructions | Cache references | Cache misses | Miss rate |
|---|---:|---:|---:|---:|---:|
| `analytical` | 112,671,786 | 221,842,471 | 1,047,534 | 20,034 | 1.9125% |
| `autodiff` | 59,946,420 | 183,208,458 | 445,143 | 17,227 | 3.8700% |
| diagnostic | 70,736,217 | 203,712,253 | 260,635 | 21,188 | 8.1294% |

The cycle and instruction counts corroborate the wall-clock winner. Cache
references varied by 4% for `analytical`, 57% for `autodiff`, and 30% for the
diagnostic, so only the analytical cache count is stable enough for detailed
interpretation; cache counters do not override wall clock.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_integrate_frozen_autodiff_jvp.json`.

Run validation and timing with:

```bash
cmake --build build-enzyme \
    --target enzyme_adaptive_frozen_trace_jvp_build
ctest --test-dir build-enzyme \
    -R '^enzyme_adaptive_frozen_trace_jvp$' --output-on-failure
taskset -c 4 \
    build-enzyme/test/ad/enzyme_adaptive_frozen_trace_jvp.enzyme/enzyme_adaptive_frozen_trace_jvp \
    --benchmark autodiff
```

## Hybrid integrand JVP inside the frozen trace

The adaptive `hybrid` candidate keeps the accepted subdivision and
Gauss-Kronrod trace walk analytical while Enzyme forward mode supplies the
local derivative of `exp(p*x)` at each frozen node. It uses the same adaptive
primal, three-panel identity check, closed-form oracle, compiler, and
10,000-workload timing scope as the preceding frozen-trace comparison.

| Candidate | Mechanism | Median ns/value+JVP | MAD ns | Peak RSS |
|---|---|---:|---:|---:|
| generic explicit-tangent trace replay | `analytical` | 2,564.1478 | 39.2610 | 2,715,648 B |
| Enzyme integrand JVP + generic trace replay | `hybrid` | 2,604.6500 | 25.3158 | 2,957,312 B |
| Enzyme through compact fixed trace | `autodiff` | 1,306.9523 | 2.4366 | 2,547,712 B |
| frozen weighted-sum central difference | diagnostic | 1,423.0742 | 6.0043 | 2,580,480 B |

The `hybrid` is 1.0158 times slower than explicit `analytical` integrand
products and 1.9929 times slower than compact trace `autodiff`. This isolates
the interface cost: replacing the explicit local derivative with Enzyme barely
changes complete wall clock, while both generic trace-walk candidates retain
the expensive panel error bookkeeping. The later smooth-adaptive tournament
should test a compact analytical/hybrid replay before selecting a production
winner.

Linux `perf stat -r 5` for the hybrid process records 98,243,999 cycles,
222,301,272 instructions, 721,774 cache references, and 23,962 cache misses
(3.3199%). Cache references varied by 29% across repetitions and remain
diagnostic. Complete-workload wall clock selects the existing compact
`autodiff` candidate for this workload class.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_integrate_frozen_hybrid_jvp.json`.

Run the hybrid benchmark with:

```bash
taskset -c 4 \
    build-enzyme/test/ad/enzyme_adaptive_frozen_trace_jvp.enzyme/enzyme_adaptive_frozen_trace_jvp \
    --benchmark hybrid
```

## Smooth adaptive-integration tournament

The complete smooth tournament removes the structural mismatch identified
above. In addition to the generic trace walkers it includes:

- compact `analytical`: explicit integrand tangents in the fixed G7K15 sum;
- compact `hybrid`: Enzyme integrand JVPs in that same fixed sum;
- `autodiff`: Enzyme through the complete fixed sum.

Every candidate first constructs and verifies the same three-panel adaptive
primal trace. All derivative kernels are checked against the independently
integrated closed form. To reduce host frequency and scheduling bias, each
tournament run rotates candidate order over 31 interleaved samples of 2,000
complete value-plus-JVP workloads. The table reports the median of three
tournament-run medians and the median within-run MAD. Reference: AMD Ryzen 9
5950X, CPU 4 pinned, Flang/LLVM 22.1.8, Enzyme `c96508349d9f`, Release `-O2`.

| Candidate | Mechanism | Median ns/value+JVP | MAD ns | Peak RSS |
|---|---|---:|---:|---:|
| generic explicit-tangent trace replay | `analytical` | 1,668.4775 | 21.3305 | 2,613,248 B |
| compact explicit-tangent trace replay | `analytical` | 1,308.4690 | 11.0055 | 2,822,144 B |
| Enzyme through compact fixed trace | `autodiff` | 1,329.2185 | 14.5675 | 2,838,528 B |
| Enzyme integrand JVP + generic trace | `hybrid` | 1,684.0765 | 18.6100 | 2,813,952 B |
| Enzyme integrand JVP + compact trace | `hybrid` | 1,345.1885 | 15.1135 | 2,785,280 B |
| frozen-trace central difference | diagnostic | 1,481.6405 | 18.0085 | 2,801,664 B |

Complete-workload wall clock selects compact `analytical`. It is 1.0159 times
faster than whole-trace `autodiff`, 1.0281 times faster than compact `hybrid`,
and 1.1323 times faster than finite differences. The generic analytical and
hybrid walkers are 1.2751 and 1.2871 times slower because they retain panel
error-estimate work. Peak RSS is effectively flat and does not select a
candidate.

Linux `perf stat -r 5` over 10,000 workloads corroborates the compact ranking:

| Candidate | Cycles | Instructions | Cache references | Cache misses | Miss rate |
|---|---:|---:|---:|---:|---:|
| compact `analytical` | 58,931,246 | 185,458,898 | 211,356 | 17,433 | 8.2482% |
| `autodiff` | 59,309,240 | 184,308,654 | 122,758 | 16,418 | 13.3743% |
| compact `hybrid` | 60,159,688 | 187,224,891 | 197,489 | 16,263 | 8.2349% |
| diagnostic | 65,283,523 | 206,303,406 | 278,105 | 19,196 | 6.9024% |

Cache-reference dispersion reached 27% for compact analytical and is
diagnostic only. The wall-clock and cycle results agree on the winner.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_integrate_smooth_tournament_jvp.json`.

Run the interleaved tournament with:

```bash
taskset -c 4 \
    build-enzyme/test/ad/enzyme_adaptive_frozen_trace_jvp.enzyme/enzyme_adaptive_frozen_trace_jvp \
    --tournament
```

## Singular adaptive-integration tournament

The singular workload is
\(I(p)=\int_0^1 e^p/\sqrt{x}\,dx\) at \(p=0.7\). Its exact continuous
derivative is \(2e^p\), while the production derivative contract replays the
six-panel frozen QAGS trace. Candidate validation therefore uses central
differences of an independently evaluated fixed G10K21 trace as the behavioral
oracle and separately bounds the frozen-trace discretization bias.

Every timed workload constructs and verifies the same QAGS primal trace and
then evaluates one JVP. The interleaved protocol is the same as the smooth
tournament: three processes, 31 rotating-order samples, and 2,000 workloads per
sample. The table reports the median of process medians and the median
within-process MAD.

| Candidate | Mechanism | Median ns/value+JVP | MAD ns | Peak RSS |
|---|---|---:|---:|---:|
| generic explicit-tangent trace replay | `analytical` | 6,753.1975 | 297.9260 | 2,605,056 B |
| compact explicit-tangent trace replay | `analytical` | 2,862.8335 | 81.6190 | 2,592,768 B |
| Enzyme through compact fixed trace | `autodiff` | 2,839.8200 | 98.9965 | 2,818,048 B |
| Enzyme integrand JVP + generic trace | `hybrid` | 7,514.2765 | 278.0390 | 2,605,056 B |
| Enzyme integrand JVP + compact trace | `hybrid` | 3,106.1170 | 165.6670 | 2,789,376 B |
| frozen-trace central difference | diagnostic | 2,924.5145 | 127.9210 | 2,748,416 B |

Complete-workload wall clock selects compact whole-trace `autodiff`. It is
1.0081 times faster than compact `analytical`, 1.0938 times faster than compact
`hybrid`, and 1.0298 times faster than the finite-difference diagnostic.
The 0.8% lead over compact analytical is smaller than the measured dispersion,
so those two candidates are a practical tie. The generic paths are materially
slower because callback, status, and panel-error traversal remain in the
derivative workload. Peak RSS differs by only about 0.23 MB and does not
override the wall-clock result.

Single `perf stat` runs over the longer per-candidate benchmark are diagnostic:

| Candidate | Cycles | Instructions | Cache references | Cache misses | Miss rate |
|---|---:|---:|---:|---:|---:|
| compact `analytical` | 980,341,483 | 2,865,255,545 | 1,291,893 | 45,041 | 3.4864% |
| `autodiff` | 982,011,124 | 2,853,703,532 | 3,797,935 | 35,971 | 0.9471% |
| compact `hybrid` | 1,055,071,626 | 2,909,701,442 | 2,192,260 | 71,017 | 3.2395% |
| diagnostic | 993,845,109 | 2,891,188,591 | 2,321,352 | 61,982 | 2.6709% |

The cycle counts confirm the near tie between compact `analytical` and
`autodiff`; cache events are not used for selection because only one counter
run was collected. The machine-readable record is
`benchmark/reference/ryzen9_5950x_integrate_singular_tournament_jvp.json`.

Run it with:

```bash
build-enzyme/test/ad/enzyme_adaptive_frozen_trace_jvp.enzyme/enzyme_adaptive_frozen_trace_jvp \
    --singular-tournament
```

## Batched fixed-integration full-Jacobian tournament

This workload batches \(B\) independent 32-point Gauss-Legendre integrals.
Each scalar output has four private active parameters, so the complete
block-diagonal Jacobian has shape \(B\times4B\). Forward construction performs
four local JVP sweeps per integral; reverse construction performs one VJP per
scalar output. Both produce every nonzero entry of the same Jacobian.

The batch sizes are 1, 4, and 16. Every candidate is checked for all 16
parameter shifts against closed-form derivatives of the exactly integrated
terms. Timings cover complete full-Jacobian construction, use 31 samples of
2,000 batches per process, and report the median of three process medians.
Reference: AMD Ryzen 9 5950X, Flang/LLVM 22.1.8, Enzyme `c96508349d9f`,
Release `-O2`.

Forward full-Jacobian wall clock:

| Mechanism | B=1 | B=4 | B=16 | B=16/B=1 |
|---|---:|---:|---:|---:|
| `analytical` | 1,280.6130 ns | 5,117.0365 ns | 20,318.8290 ns | 15.8665 |
| `autodiff` | 1,788.8005 ns | 6,995.6110 ns | 28,073.2040 ns | 15.6940 |
| `hybrid` | 1,272.2220 ns | 5,095.4910 ns | 20,674.0080 ns | 16.2504 |
| diagnostic | 2,325.1010 ns | 9,263.2680 ns | 37,533.0330 ns | 16.1425 |

Reverse full-Jacobian wall clock:

| Mechanism | B=1 | B=4 | B=16 | B=16/B=1 |
|---|---:|---:|---:|---:|
| `analytical` | 319.2215 ns | 1,199.5500 ns | 4,884.7135 ns | 15.3019 |
| `autodiff` | 450.3835 ns | 1,871.6915 ns | 6,994.6890 ns | 15.5305 |
| `hybrid` | 340.7270 ns | 1,315.5990 ns | 5,011.5825 ns | 14.7085 |
| diagnostic | 2,323.7430 ns | 9,115.6250 ns | 36,628.8195 ns | 15.7629 |

Complete-workload wall clock selects reverse `analytical`. At \(B=16\), it is
4.1597 times faster than forward `analytical`, 1.4319 times faster than reverse
`autodiff`, 1.0260 times faster than reverse `hybrid`, and 7.4987 times faster
than finite differences. The reason is dimensional rather than terminological:
one scalar output with four active inputs needs one reverse sweep but four
forward sweeps. `analytical` and `hybrid` forward timings are a practical tie
at small batches, but reverse analytical remains the lowest full-Jacobian wall
clock at every measured batch size.

Peak RSS at \(B=16\) ranges from 2,859,008 to 3,031,040 bytes across all
candidates. For selected analytical reverse it rises only from 2,850,816 bytes
at \(B=1\) to 3,002,368 bytes at \(B=16\); the implementation streams batch
items and stores no batch-sized derivative tape.

Linux `perf stat -r 3` over the CTest-launched \(B=16\) processes gives:

| Candidate | Cycles | Instructions | Cache references | Cache misses |
|---|---:|---:|---:|---:|
| forward `analytical` | 5,685,837,952 | 20,081,441,392 | 7,736,492 | 256,097 |
| reverse `analytical` | 1,394,670,950 | 5,034,745,727 | 1,887,203 | 198,429 |
| reverse `autodiff` | 2,037,786,872 | 7,288,509,995 | 2,299,377 | 227,362 |
| reverse `hybrid` | 1,463,346,989 | 5,198,218,860 | 1,501,620 | 208,226 |
| reverse diagnostic | 10,231,175,773 | 35,197,916,516 | 40,258,188 | 278,969 |

Forward analytical uses 4.0767 times as many cycles and 3.9886 times as many
instructions as reverse analytical, matching the wall-clock explanation.
Cache-reference dispersion is high for forward analytical and the diagnostic,
so cache counters remain supporting evidence and do not override wall clock.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_integrate_batched_full_jacobian.json`.
The benchmark is launched through CTest by setting `FORTNUM_BATCH_ACTION`,
`FORTNUM_BATCH_CANDIDATE`, and `FORTNUM_BATCH_SIZE`; the ordinary test path
continues to run only independent validation.

## Hybrid ODE forward sensitivity

The first ODE `hybrid` candidate uses Enzyme forward mode only at the local RHS
boundary. `ode_integrate_jvp` remains the analytical composition layer and
advances

\[
\dot s = J_f(t,y)s + f_k(t,y)\dot k
\]

over the Cash–Karp step schedule recorded by the adaptive primal. The workload
uses \(y'=-ky\), \(y(0)=1.3\), \(t_1=2\), and \(k\) near 0.7. The independent
oracles are the closed-form sensitivity
\(-t_1y_0\exp(-kt_1)\) and central differences of complete primal solves.

Every timed analytical or hybrid workload includes adaptive primal integration
and one frozen-trace parameter JVP. The diagnostic includes a base primal plus
the two perturbed primal solves. Results are medians of three processes, each
with 31 samples of 2,000 complete workloads. Reference: AMD Ryzen 9 5950X,
Flang/LLVM 22.1.8, Enzyme `c96508349d9f`, Release `-O2`.

| Candidate | Mechanism | Median ns/value+JVP | MAD ns | Peak RSS |
|---|---|---:|---:|---:|
| explicit RHS tangent + frozen trace | `analytical` | 9,322.0850 | 177.4240 | 3,145,728 B |
| Enzyme RHS JVP + frozen trace | `hybrid` | 9,690.5995 | 364.6070 | 3,035,136 B |
| complete-solve central difference | diagnostic | 14,965.7690 | 345.1100 | 3,076,096 B |

Complete-workload wall clock selects `analytical`. It is 1.0395 times faster
than `hybrid` and 1.6054 times faster than finite differences. The hybrid is
1.5444 times faster than finite differences, showing that the local Enzyme
boundary removes most derivative-maintenance work without differentiating
adaptive control flow. This scalar one-input, one-output item does not establish
a forward-versus-reverse scaling verdict; the later many-parameter trajectory
tournament must do that.

Linux `perf stat -r 3` over CTest-launched benchmark processes records:

| Candidate | Cycles | Instructions | Cache references | Cache misses |
|---|---:|---:|---:|---:|
| `analytical` | 2,689,394,152 | 7,819,202,627 | 14,776,466 | 383,372 |
| `hybrid` | 2,680,804,009 | 7,775,468,150 | 11,247,531 | 332,953 |
| diagnostic | 4,282,025,762 | 10,540,903,123 | 11,526,373 | 468,168 |

Analytical and hybrid cycle counts are effectively tied even though analytical
has the lower wall-clock median. Cache-reference dispersion is 35% for
analytical and 7% for hybrid, so counters are diagnostic only and wall clock
selects the candidate.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_ode_hybrid_forward_sensitivity.json`.

## Analytical Cash–Karp discrete adjoint

`ode_integrate_vjp` is the exact transpose of the analytical Cash–Karp tangent
step, walked backward over the accepted-step trace selected by the primal. The
direct behavioral oracle uses a two-state linear system
\(\dot y=Ay\) and verifies

\[
J^Tu = \exp(A^Tt_1)u.
\]

This is independent of the existing tangent/adjoint dot-product identity.
Complete-solve central differences provide a second oracle.

The requested derivative object is one scalar-terminal-objective VJP. Each
timed candidate includes the adaptive primal trajectory. Reverse analytical
performs one discrete-adjoint trace walk; forward reconstruction performs two
tangent walks, one per initial-state input; the diagnostic evaluates the base
trajectory and four perturbed trajectories. Results are medians of three
processes, each with 31 samples of 1,000 workloads. Reference: AMD Ryzen 9
5950X, GNU Fortran 16.1.1, Release.

| Candidate | Mechanism | Median ns/primal+VJP | MAD ns | Peak RSS |
|---|---|---:|---:|---:|
| one reverse discrete adjoint | `analytical` | 34,270.5560 | 1,023.0970 | 5,345,280 B |
| two forward tangent sweeps | `analytical` | 58,768.7470 | 1,074.3040 | 5,382,144 B |
| complete-solve central difference | diagnostic | 62,939.3870 | 1,150.4480 | 5,292,032 B |

Complete-workload wall clock selects the reverse discrete adjoint. It is
1.7149 times faster than forward reconstruction and 1.8365 times faster than
finite differences. Peak RSS is equal within 91 KB and does not affect the
selection.

Linux `perf stat -r 3` over CTest-launched benchmark processes records:

| Candidate | Cycles | Instructions | Cache references | Cache misses |
|---|---:|---:|---:|---:|
| reverse `analytical` | 4,926,950,485 | 18,091,883,235 | 22,343,502 | 382,233 |
| forward reconstruction | 8,317,476,053 | 28,965,614,332 | 23,260,473 | 445,412 |
| diagnostic | 8,850,465,175 | 30,614,064,380 | 18,595,675 | 411,590 |

Cycles and instructions corroborate the wall-clock winner. Cache counters are
supporting evidence only. This two-input, one-objective workload demonstrates
the expected reverse advantage for a VJP; state-, parameter-, and objective-
count scaling remains the responsibility of the later many-parameter
trajectory tournament.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_ode_discrete_adjoint.json`.

## Cash–Karp parameter-adjoint scaling

`ode_integrate_parameter_vjp` extends the analytical discrete adjoint with a
stage callback that accumulates

\[
\left(\frac{\partial f}{\partial p}\right)^T\bar f_i
\]

for every Runge–Kutta stage cotangent \(\bar f_i\). The callback owns the
parameter layout, so the reverse trace does not form a state-by-parameter
Jacobian.

The scaling workload is
\(y'=-\operatorname{mean}(p)y\), \(y(0)=1.3\), \(t_1=2\), with one scalar
terminal objective and 1, 4, or 16 active parameters. Every gradient component
has the closed form

\[
-\bar y_1 t_1 y_0
\exp[-\operatorname{mean}(p)t_1]/n_p.
\]

The reverse candidate performs one primal and one adjoint trace regardless of
parameter count. Forward performs one primal and one tangent trace per
parameter. The diagnostic evaluates a base primal and two complete perturbed
primal solves per parameter. Results are medians of three processes, each with
31 samples of 200 complete gradients. Reference: AMD Ryzen 9 5950X, GNU
Fortran 16.1.1, Release.

| Candidate | Mechanism | 1 parameter | 4 parameters | 16 parameters |
|---|---|---:|---:|---:|
| one parameter-accumulating adjoint | `analytical` reverse | 21,043.5950 ns | 20,967.5000 ns | 22,712.5350 ns |
| one tangent per parameter | `analytical` forward | 18,128.4000 ns | 53,328.5150 ns | 193,192.9500 ns |
| complete-solve central differences | diagnostic | 18,766.0550 ns | 55,831.3800 ns | 206,896.3650 ns |

Wall clock selects forward for one parameter: it is 1.1608 times faster than
reverse. Reverse wins at four parameters by 2.5434 times and at sixteen by
8.5060 times. From one to sixteen parameters, reverse grows only 1.0793 times,
while forward grows 10.6569 times. This is the measured forward/reverse
crossover the selector must preserve; neither mode is a universal default.

Peak RSS stays between 5.32 and 5.50 MB for all candidates and parameter
counts. The parameter gradient itself is the only storage proportional to
\(n_p\); there is no parameter-sized trajectory tape.

Selected `perf stat -r 3` results:

| Candidate | Parameters | Cycles | Instructions | Cache references | Cache misses |
|---|---:|---:|---:|---:|---:|
| reverse | 1 | 606,161,833 | 1,889,467,959 | 5,010,051 | 250,626 |
| forward | 1 | 549,756,725 | 1,600,964,051 | 5,290,049 | 252,998 |
| reverse | 16 | 666,099,573 | 2,371,565,182 | 5,328,395 | 271,439 |
| forward | 16 | 5,492,501,896 | 16,738,933,055 | 6,656,145 | 433,791 |
| diagnostic | 16 | 5,899,163,995 | 17,719,627,576 | 7,548,423 | 358,595 |

Cycles and instructions corroborate both the one-parameter forward winner and
the sixteen-parameter reverse winner. Cache counts remain supporting evidence;
complete-gradient wall clock determines selection.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_ode_parameter_adjoint_scaling.json`.

## Vector-root candidate tournament

The coupled vector tournament uses

```text
F1 = x1^2 + x2 - p1
F2 = x1 + x2^2 - p2
```

and compares the same four mechanisms as the scalar tournament. The
`analytical` and `hybrid` candidates differentiate the residual equation and
solve the implicit tangent or adjoint system. The `autodiff` candidate uses
Enzyme on a fixed 12-step coupled Newton trace. The finite-difference
diagnostic differences complete fixed-step solves.

Every candidate is checked against fresh complete solves at independently
perturbed parameters. Maximum JVP errors were `9.7185e-12` for the first three
candidates and zero for the diagnostic. Maximum VJP errors were `3.7337e-11`,
`3.7337e-11`, `3.7337e-11`, and `1.6653e-12`, respectively.

The benchmark applies 100,000 varying directions or cotangents per sample.
Reference: AMD Ryzen 9 5950X, CPU 4 pinned, Flang/LLVM 22.1.8, Enzyme
`c96508349d9f`, Release `-O2`, 15 samples after three warmups.

| JVP candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| implicit with explicit residual products | `analytical` | 80.7424 | 0.4353 | 2,424,832 B |
| implicit with Enzyme residual products | `hybrid` | 80.4377 | 0.3774 | 2,224,128 B |
| Enzyme through 12 Newton steps | `autodiff` | 230.2933 | 1.7658 | 2,465,792 B |
| complete-solve central difference | diagnostic | 386.5402 | 1.9719 | 2,273,280 B |

The `hybrid` JVP has the lowest median runtime and uses 200,704 fewer bytes of
maximum observed RSS than `analytical`. It is 1.0038 times faster than
`analytical`, 2.8630 times faster than `autodiff` through the iterations, and
4.8055 times faster than the finite-difference diagnostic.

| VJP candidate | Mechanism | Median ns/VJP | MAD ns/VJP | Peak RSS |
|---|---|---:|---:|---:|
| implicit with explicit residual products | `analytical` | 83.3377 | 0.8708 | 2,535,424 B |
| implicit with Enzyme residual products | `hybrid` | 82.5682 | 0.5024 | 2,510,848 B |
| Enzyme through 12 Newton steps | `autodiff` | 228.8764 | 1.0845 | 2,457,600 B |
| componentwise complete-solve differences | diagnostic | 392.2317 | 3.4115 | 2,437,120 B |

The reverse-mode `hybrid` VJP has the lowest median runtime and uses 24,576
fewer bytes of maximum observed RSS than `analytical`. It is 1.0093 times
faster than `analytical`, 2.7720 times faster than reverse `autodiff` through
the iterations, and 4.7504 times faster than the diagnostic.

The machine-readable records are
`benchmark/reference/ryzen9_5950x_vector_root_tournament_jvp.json` and
`benchmark/reference/ryzen9_5950x_vector_root_tournament_vjp.json`.

Run validation and timing with:

```bash
cmake --build build-enzyme --target enzyme_vector_root_hybrid_build \
    enzyme_vector_root_vjp_hybrid_build
ctest --test-dir build-enzyme -R '^enzyme_vector_root(_vjp)?_hybrid$' \
    --output-on-failure
taskset -c 4 \
    build-enzyme/test/ad/enzyme_vector_root_hybrid.enzyme/enzyme_vector_root_hybrid \
    --benchmark
taskset -c 4 \
    build-enzyme/test/ad/enzyme_vector_root_vjp_hybrid.enzyme/enzyme_vector_root_vjp_hybrid \
    --benchmark
```

## Analytical differentiation under a fixed-bound integral

For inactive bounds, `integrate_fixed_jvp` applies

```text
d/dp integral_a^b f(x,p) dx = integral_a^b (df/dp)(x,p) dx.
```

The independent test uses `f(x,p)=exp(p*x)` on `[0,1]`. It compares the
integrated tangent with both the closed-form derivative and a central
difference of two complete primal integrations.

The benchmark varies `p` over 101 values and computes 10,000 derivatives per
sample. `analytical` performs one adaptive integration of `df/dp`.
`diagnostic` performs two complete adaptive primal integrations at `p+h` and
`p-h`. Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1,
Release build, 15 samples after three warmups.

| Candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| integrate the analytical tangent | `analytical` | 430.3008 | 3.5697 | 12,550,144 B |
| central difference of complete integrals | diagnostic | 820.7775 | 14.6786 | 12,177,408 B |

The `analytical` JVP is 1.9075 times faster in complete-workload wall clock,
saving 390.4767 ns per derivative. It uses 372,736 more bytes of maximum
observed process RSS, so this workload has a runtime-memory tradeoff. Wall
clock is the primary metric here, and selects `analytical`.

Linux `perf stat -r 5` over the same 10,000-call process gives:

| Candidate | Cycles | Instructions | Cache references | Cache misses | Miss rate |
|---|---:|---:|---:|---:|---:|
| `analytical` | 649,863,590 | 1,499,721,281 | 25,971,138 | 4,477,666 | 17.2409% |
| diagnostic | 1,111,056,138 | 2,366,450,831 | 33,800,218 | 4,597,997 | 13.6035% |

Although the diagnostic has a lower miss rate, it performs 1.3015 times as
many cache references and 1.0269 times as many absolute cache misses. It also
executes 1.5779 times as many instructions and consumes 1.7097 times as many
cycles. Cache counters include the identical `fo exec` launcher overhead and
are therefore supporting evidence; measured complete-workload wall clock
remains decisive.

This case has one active input and one scalar output, so it does not establish
a forward/reverse crossover. The roadmap now requires input/output-count
scaling and cache evidence for tournaments where those dimensions vary.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_integrate_fixed_jvp.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_integrate_fixed_jvp analytical
taskset -c 4 fo exec bench_integrate_fixed_jvp diagnostic
fo exec bench_integrate_fixed_jvp analytical --peak-rss
fo exec bench_integrate_fixed_jvp diagnostic --peak-rss
perf stat -r 5 -e cycles,instructions,cache-references,cache-misses \
    taskset -c 4 fo exec bench_integrate_fixed_jvp analytical
```

## Analytical moving-lower-bound term

`integrate_moving_lower_jvp` implements the Leibniz rule

```text
dI/dp = integral_a^b df/dp dx - f(a,p) da/dp
```

for an active lower bound and inactive upper bound. The independent test uses
`f(x,p)=exp(p*x)`, `a(p)=0.1*p`, and `b=1`. It compares against both a
closed-form derivative and complete primal integrations at `p+h` and `p-h`,
including the perturbed lower limit.

The benchmark varies `p` over 101 values and computes 10,000 derivatives per
sample. `analytical` integrates the tangent once and adds the endpoint term.
The finite-difference diagnostic performs two complete adaptive integrations.
Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| tangent integral plus lower endpoint | `analytical` | 428.1954 | 2.6470 | 12,304,384 B |
| complete moving-bound differences | diagnostic | 868.7511 | 3.4445 | 12,304,384 B |

The `analytical` candidate is 2.0289 times faster in complete-workload wall
clock, saving 440.5557 ns per derivative, with identical maximum observed
process RSS.

Linux `perf stat -r 5` over the same process gives:

| Candidate | Cycles | Instructions | Cache references | Cache misses | Miss rate |
|---|---:|---:|---:|---:|---:|
| `analytical` | 658,611,224 | 1,538,209,367 | 25,879,115 | 4,602,608 | 17.7850% |
| diagnostic | 946,269,077 | 2,391,248,616 | 29,936,298 | 4,577,791 | 15.2918% |

The diagnostic uses 1.4368 times as many cycles, 1.5546 times as many
instructions, and 1.1568 times as many cache references. The analytical path
has 0.54% more absolute cache misses despite doing less total work;
complete-workload wall clock still selects it.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_integrate_moving_lower_jvp.json`.

## Analytical moving-upper-bound term

`integrate_moving_upper_jvp` implements

```text
dI/dp = integral_a^b df/dp dx + f(b,p) db/dp
```

for an inactive lower bound and active upper bound. The independent test uses
`f(x,p)=exp(p*x)`, `a=0`, and `b(p)=0.5+0.1*p`. It compares with a closed-form
derivative and complete primal integrations whose upper limits are perturbed.

The benchmark computes 10,000 derivatives per sample over 101 parameter
values. Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1,
Release build, 15 samples after three warmups.

| Candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| tangent integral plus upper endpoint | `analytical` | 458.6955 | 6.3670 | 12,312,576 B |
| complete moving-bound differences | diagnostic | 775.7548 | 13.7940 | 12,537,856 B |

The `analytical` candidate is 1.6912 times faster in complete-workload wall
clock, saves 317.0593 ns per derivative, and uses 225,280 fewer bytes of
maximum observed process RSS.

Linux `perf stat -r 5` gives:

| Candidate | Cycles | Instructions | Cache references | Cache misses | Miss rate |
|---|---:|---:|---:|---:|---:|
| `analytical` | 670,655,351 | 1,561,690,214 | 26,954,668 | 4,801,576 | 17.8135% |
| diagnostic | 912,613,852 | 2,407,807,332 | 27,799,974 | 4,796,019 | 17.2519% |

The diagnostic uses 1.3608 times as many cycles, 1.5418 times as many
instructions, and 1.0314 times as many cache references. Absolute cache misses
differ by only 0.12%; wall clock and peak memory both select `analytical`.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_integrate_moving_upper_jvp.json`.

## Hybrid vector-root JVP

This comparison computes the same two-component JVP of the converged root for

```text
F1(x,p) = x1^2 + x2 - p1
F2(x,p) = x1 + x2^2 - p2.
```

Both candidates use the analytical implicit tangent boundary. The
`analytical` candidate supplies explicit `F_x` and `F_p*tp`. The `hybrid`
candidate uses six statically resolved Enzyme forward sweeps: four produce the
two-by-two state Jacobian and two produce the contracted parameter residual
JVP. Neither candidate differentiates nonlinear-solver iterations.

The independent oracle centrally differences complete 12-step Newton solves
at `p+h*tp` and `p-h*tp`. The maximum absolute error for both candidates was
`9.7185e-12`.

The benchmark applies 2,000,000 directions per sample while varying both root
components and preserving the residual equation. Reference: AMD Ryzen 9 5950X,
CPU 4 pinned, Flang/LLVM 22.1.8, Enzyme `c96508349d9f`, Release `-O2`, 15
samples after three warmups.

| Candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| explicit residual products | `analytical` | 83.8493 | 1.1293 | 2,490,368 B |
| Enzyme forward residual products | `hybrid` | 81.8199 | 0.6058 | 2,658,304 B |

Here “1.0248 times faster” means
`83.8493 / 81.8199 = 1.0248`: the isolated `hybrid` JVP saves 2.0294 ns
relative to the `analytical` JVP. The comparison is not against a primal root
solve. Hybrid has the lower isolated runtime, while analytical uses 167,936
fewer bytes of maximum observed self-process RSS. There is therefore no
production selection from this microbenchmark alone; the roots tournament
must resolve the tradeoff for a complete application workload.

Peak memory is the maximum across five separately launched candidates,
measured with `getrusage(RUSAGE_SELF)`. The machine-readable record is
`benchmark/reference/ryzen9_5950x_vector_root_hybrid_jvp.json`.

Run validation and timing with:

```bash
cmake --build build-enzyme --target enzyme_vector_root_hybrid_build
ctest --test-dir build-enzyme -R '^enzyme_vector_root_hybrid$' \
    --output-on-failure
taskset -c 4 \
    build-enzyme/test/ad/enzyme_vector_root_hybrid.enzyme/enzyme_vector_root_hybrid \
    --benchmark
```

## Hybrid vector-root VJP

This comparison computes the same two-parameter VJP for the two-state residual
used by the vector-root JVP benchmark. Both candidates use the analytical
implicit adjoint boundary: assemble `F_x`, solve `F_x^T*lambda=u`, then compute
`-F_p^T*lambda`. The `analytical` candidate supplies explicit residual
partials. The `hybrid` candidate obtains a four-scalar gradient from each
residual component with Enzyme reverse mode, both for state-Jacobian assembly
and for the contracted parameter VJP. Neither differentiates solver
iterations.

The independent oracle centrally differences the complete objective
`L(p)=u^T*x*(p)` component by component. Every perturbed root comes from a
fresh 12-step Newton solve. The maximum absolute error for both candidates was
`3.7337e-11`.

The benchmark applies 2,000,000 cotangents per sample while varying both root
components and preserving the residual equation. Reference: AMD Ryzen 9 5950X,
CPU 4 pinned, Flang/LLVM 22.1.8, Enzyme `c96508349d9f`, Release `-O2`, 15
samples after three warmups.

| Candidate | Mechanism | Median ns/VJP | MAD ns/VJP | Peak RSS |
|---|---|---:|---:|---:|
| explicit residual products | `analytical` | 73.8151 | 0.5271 | 2,441,216 B |
| Enzyme reverse residual products | `hybrid` | 74.1869 | 0.9943 | 2,420,736 B |

Here “1.0050 times faster” means
`74.1869 / 73.8151 = 1.0050`: the isolated `analytical` VJP saves 0.3718 ns
relative to `hybrid`. The comparison is not against a primal root solve.
Analytical has the lower isolated runtime, while hybrid uses 20,480 fewer
bytes of maximum observed self-process RSS. There is no production selection
from this microbenchmark alone; the roots tournament must resolve the tradeoff
for a complete application workload.

Peak memory is the maximum across five separately launched candidates,
measured with `getrusage(RUSAGE_SELF)`. The machine-readable record is
`benchmark/reference/ryzen9_5950x_vector_root_hybrid_vjp.json`.

Run validation and timing with:

```bash
cmake --build build-enzyme --target enzyme_vector_root_vjp_hybrid_build
ctest --test-dir build-enzyme -R '^enzyme_vector_root_vjp_hybrid$' \
    --output-on-failure
taskset -c 4 \
    build-enzyme/test/ad/enzyme_vector_root_vjp_hybrid.enzyme/enzyme_vector_root_vjp_hybrid \
    --benchmark
```

## Generic analytical implicit adjoint boundary for scalar roots

The generic adjoint boundary accepts a converged scalar root and a callback
that returns `f_x` and the contracted residual VJP `f_p^T*u`. It applies the
analytical implicit factor `-1/f_x`, without differentiating root-solver
iterations. The callback output uses an explicit extent tied to `size(p)` so
the boundary does not require an assumed-shape result descriptor.

The independent oracle differentiates the scalar objective
`L(p) = 1.3*x*(p)`, where `x*(p)` is obtained by completely re-solving
`x^3 + p1*x - p2 = 0`. Each parameter component is checked with a central
difference of complete objective evaluations.

The benchmark applies 2,000,000 root cotangents per sample. Both rows compute
the same two-component scalar-root VJP. The baseline supplies an assembled
residual gradient to the existing low-level `root_vjp`; the candidate uses the
new contracted callback boundary.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Mechanism | Median ns/VJP | MAD ns/VJP | Peak RSS |
|---|---|---:|---:|---:|
| assembled residual gradient | `analytical` | 15.7552 | 0.0287 | 12,304,384 B |
| contracted callback boundary | `analytical` | 35.6404 | 0.2652 | 12,521,472 B |

Here “2.2621 times slower” means
`35.6404 / 15.7552 = 2.2621`: it compares the callback-boundary row with the
assembled-gradient row for the same VJP, not with a primal root solve. The
assembled path remains the isolated-kernel selection. The generic boundary
costs 19.8852 ns/VJP and 217,088 bytes more maximum observed self-process RSS,
while providing the reverse residual-product interface needed for later
`hybrid` composition.

Peak memory is the maximum across five separately launched processes per
candidate, measured with `getrusage(RUSAGE_SELF)`. The machine-readable record
is
`benchmark/reference/ryzen9_5950x_scalar_root_adjoint_boundary.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_scalar_root_adjoint assembled
taskset -c 4 fo exec bench_scalar_root_adjoint boundary
fo exec bench_scalar_root_adjoint assembled --peak-rss
fo exec bench_scalar_root_adjoint boundary --peak-rss
```

## Generic analytical implicit tangent boundary for vector roots

The vector-root tangent boundary accepts a converged root and a callback that
returns the state Jacobian `F_x` and contracted residual tangent `F_p*tp`. It
then solves `F_x*dx = -(F_p*tp)` analytically, without differentiating the
root-solver iterations. Factorization reuse and additional reliability options
remain separate ROADMAP items.

The independent oracle uses the nonlinear two-state system
`F1 = x1^2 + x2 - p1`, `F2 = x1 + x2^2 - p2`. It completely re-solves at
`p + h*tp` and `p - h*tp` and compares the resulting central-difference
direction with the boundary JVP.

The benchmark applies 2,000,000 directions per sample to a two-state,
two-parameter residual. Both rows compute the same vector-root JVP. The
baseline uses preassembled `F_x` and `F_p`; the candidate evaluates `F_x` and
the contracted `F_p*tp` through the new callback boundary.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| preassembled residual Jacobians | `analytical` | 75.6729 | 0.7939 | 12,537,856 B |
| contracted callback boundary | `analytical` | 107.5941 | 0.8444 | 12,513,280 B |

Here “1.4218 times slower” means
`107.5941 / 75.6729 = 1.4218`: it compares the callback-boundary row with the
preassembled row for the same JVP, not with a primal vector-root solve. The
preassembled path remains the isolated-kernel selection. The boundary costs
31.9212 ns/JVP while providing the operator interface required for later
`hybrid` residual products. Its maximum observed self-process RSS is 24,576
bytes lower.

Peak memory is the maximum across five separately launched processes per
candidate, measured with `getrusage(RUSAGE_SELF)`. The machine-readable record
is
`benchmark/reference/ryzen9_5950x_vector_root_tangent_boundary.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_vector_root_tangent assembled
taskset -c 4 fo exec bench_vector_root_tangent boundary
fo exec bench_vector_root_tangent assembled --peak-rss
fo exec bench_vector_root_tangent boundary --peak-rss
```

## Generic analytical implicit adjoint boundary for vector roots

The vector-root adjoint boundary accepts separate callbacks for the state
Jacobian `F_x` and the contracted parameter VJP. It solves
`F_x^T*lambda = u`, evaluates `F_p^T*lambda`, and negates that result. Solver
iterations remain inactive. The separation follows the mathematical operator
boundary and permits a later autodiff parameter-VJP callback without requiring
a full `F_p`.

The independent oracle uses the nonlinear two-state system
`F1 = x1^2 + x2 - p1`, `F2 = x1 + x2^2 - p2`. For the scalar objective
`L(p) = u^T*x*(p)`, it completely re-solves each positive and negative
parameter perturbation and checks both components of the boundary VJP.

The benchmark applies 2,000,000 cotangents per sample to a two-state,
two-parameter residual. Both rows compute the same vector-root VJP. The
baseline uses preassembled `F_x` and full `F_p`; the candidate evaluates `F_x`
and the contracted `F_p^T*lambda` through the new boundary.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Mechanism | Median ns/VJP | MAD ns/VJP | Peak RSS |
|---|---|---:|---:|---:|
| preassembled residual Jacobians | `analytical` | 105.9837 | 1.4521 | 12,316,672 B |
| contracted callback boundary | `analytical` | 94.1294 | 0.6195 | 12,472,320 B |

Here “1.1259 times faster” means
`105.9837 / 94.1294 = 1.1259`: it compares the callback-boundary row with the
preassembled row for the same VJP, not with a primal vector-root solve. The
callback boundary is the isolated-kernel runtime selection, saving
11.8543 ns/VJP by avoiding a full parameter-Jacobian product after the adjoint
solve. Its maximum observed self-process RSS is 155,648 bytes higher.

Peak memory is the maximum across five separately launched processes per
candidate, measured with `getrusage(RUSAGE_SELF)`. The machine-readable record
is
`benchmark/reference/ryzen9_5950x_vector_root_adjoint_boundary.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_vector_root_adjoint assembled
taskset -c 4 fo exec bench_vector_root_adjoint boundary
fo exec bench_vector_root_adjoint assembled --peak-rss
fo exec bench_vector_root_adjoint boundary --peak-rss
```

## Converged vector-root Jacobian reuse

The analytical-Jacobian vector-root solver can optionally return the Jacobian
that its callback evaluated with the last accepted iterate. On convergence,
this is the converged `F_x` required by analytical and hybrid implicit
derivative products. The implementation copies it only at solver exit, not
after every Newton step.

The independent behavioral oracle checks the returned matrices for the
Rosenbrock-gradient, Powell-singular, and circle-line systems against central
finite differences of residual values at the converged roots. The maximum
absolute matrix-entry error was `1.0702e-8`.

The complete benchmark workload solves an eight-state dense nonlinear system
with trigonometric residual terms and an analytical Jacobian. The baseline
solves and then reevaluates the residual/Jacobian callback at the converged
root. The reuse candidate requests the already evaluated converged Jacobian
from the same solve. Each sample contains 5,000 complete solve-plus-Jacobian
workloads.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Mechanism | Median ns/workload | MAD ns/workload | Peak RSS |
|---|---|---:|---:|---:|
| solve then reevaluate Jacobian | `analytical` | 2,042.0032 | 18.8718 | 12,521,472 B |
| return and reuse converged Jacobian | `analytical` | 1,963.0000 | 10.2972 | 12,316,672 B |

Here “1.0402 times faster” means
`2042.0032 / 1963.0000 = 1.0402`: reuse saves 79.0032 ns for the same complete
solve-plus-Jacobian workload. It is not a comparison with an isolated
Jacobian callback. Reuse also lowers maximum observed self-process RSS by
204,800 bytes and is the measured selection for this workload.

Peak memory is the maximum across five separately launched candidates,
measured with `getrusage(RUSAGE_SELF)`. The machine-readable record is
`benchmark/reference/ryzen9_5950x_multiroot_jacobian_reuse.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_multiroot_jacobian_reuse reevaluate
taskset -c 4 fo exec bench_multiroot_jacobian_reuse reuse
fo exec bench_multiroot_jacobian_reuse reevaluate --peak-rss
fo exec bench_multiroot_jacobian_reuse reuse --peak-rss
```

## Vector-root JVP factorization reuse

`multiroot_jvp_factored` accepts the compact LU arrays and pivots obtained by
factoring the converged state Jacobian once. It forms the contracted residual
tangent and applies the analytical implicit solve without refactorizing
`F_x` for every parameter direction.

The independent oracle compares the factored JVP with a central difference of
complete nonlinear root solves at `p+h*tp` and `p-h*tp`. Its maximum absolute
error was `6.0430e-12`.

The benchmark uses a dense 16-state Jacobian and parameter Jacobian. The
baseline `multiroot_jvp` factors the same converged Jacobian for every
direction. The reuse candidate factors once before timing and applies those
factors to 200,000 varying directions per sample.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| refactor converged Jacobian per JVP | `analytical` | 1,052.7729 | 17.6672 | 12,484,608 B |
| reuse converged Jacobian LU | `analytical` | 356.1775 | 1.3245 | 12,312,576 B |

Here “2.9558 times faster” means
`1052.7729 / 356.1775 = 2.9558`: factor reuse saves 696.5954 ns for the same
JVP. Reuse also lowers maximum observed self-process RSS by 172,032 bytes and
is the measured selection for this workload.

Peak memory is the maximum across five separately launched candidates,
measured with `getrusage(RUSAGE_SELF)`. The machine-readable record is
`benchmark/reference/ryzen9_5950x_multiroot_jvp_factorization.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_multiroot_jvp_factorization refactor
taskset -c 4 fo exec bench_multiroot_jvp_factorization reuse
fo exec bench_multiroot_jvp_factorization refactor --peak-rss
fo exec bench_multiroot_jvp_factorization reuse --peak-rss
```

## Vector-root VJP transpose-factorization reuse

`multiroot_vjp_factored` accepts a compact LU of the transposed converged state
Jacobian. It solves the analytical adjoint equation and contracts
`-F_p^T*lambda` without refactorizing `F_x^T` for every root cotangent.

The independent oracle compares the factored VJP with componentwise central
differences of the complete scalar objective `L(p)=u^T*x*(p)`. Every perturbed
root is obtained from a complete nonlinear solve. The maximum absolute error
was `3.7637e-11`.

The benchmark uses a dense 16-state Jacobian and parameter Jacobian. The
baseline `multiroot_vjp` transposes and factors the same converged Jacobian for
every cotangent. The reuse candidate factors the transpose once before timing
and applies those factors to 200,000 varying cotangents per sample.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Mechanism | Median ns/VJP | MAD ns/VJP | Peak RSS |
|---|---|---:|---:|---:|
| refactor transposed Jacobian per VJP | `analytical` | 1,187.7362 | 17.0501 | 12,320,768 B |
| reuse transposed Jacobian LU | `analytical` | 398.2992 | 3.2079 | 12,505,088 B |

Here “2.9820 times faster” means
`1187.7362 / 398.2992 = 2.9820`: factor reuse saves 789.4370 ns for the same
VJP. Reuse is the isolated runtime winner, while refactorization uses 184,320
fewer bytes of maximum observed self-process RSS. This microbenchmark
therefore leaves production selection to the complete roots workload.

Peak memory is the maximum across five separately launched candidates,
measured with `getrusage(RUSAGE_SELF)`. The machine-readable record is
`benchmark/reference/ryzen9_5950x_multiroot_vjp_factorization.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_multiroot_vjp_factorization refactor
taskset -c 4 fo exec bench_multiroot_vjp_factorization reuse
fo exec bench_multiroot_vjp_factorization refactor --peak-rss
fo exec bench_multiroot_vjp_factorization reuse --peak-rss
```

## Analytical implicit tangent product for fixed points

For a converged fixed point `x = G(x,p)`, the new analytical JVP solves
`(I - G_x)*dx = G_p*tp`. It treats the fixed-point iterations as inactive and
reuses the converged state plus local map derivatives supplied by the caller.

The independent oracle iterates the nonlinear map

```text
G1 = tanh(0.2*x1 + 0.1*x2 + p1)
G2 = tanh(0.05*x1 + 0.25*x2 + p2)
```

to convergence at `p + h*tp` and `p - h*tp`. Its central difference is
compared with the implicit JVP.

The benchmark compares that same two-state, two-parameter JVP. The analytical
row uses a reusable converged fixed point and map derivatives. The reference
row performs two complete fixed-point solves for every direction. Analytical
uses 2,000,000 directions per sample; the slower reference uses 50,000.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Kind | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| analytical implicit solve | `analytical` | 97.8032 | 2.4065 | 12,554,240 B |
| two complete re-solves | reference oracle | 1,224.4016 | 10.1642 | 12,550,144 B |

Here “12.5190 times faster” means
`1224.4016 / 97.8032 = 12.5190`: it compares the analytical row with the
complete-re-solve reference for the same JVP, not with a standalone primal
iteration. Analytical saves 1,126.5984 ns/JVP and is the runtime selection; its
maximum observed self-process RSS is 4,096 bytes higher.

Peak memory is the maximum across five separately launched processes per
candidate, measured with `getrusage(RUSAGE_SELF)`. The machine-readable record
is `benchmark/reference/ryzen9_5950x_fixed_point_tangent.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_fixed_point_tangent analytical
taskset -c 4 fo exec bench_fixed_point_tangent reference
fo exec bench_fixed_point_tangent analytical --peak-rss
fo exec bench_fixed_point_tangent reference --peak-rss
```

## Analytical implicit adjoint product for fixed points

For a converged fixed point `x = G(x,p)`, the analytical VJP solves
`(I - G_x)^T*lambda = u` and returns `G_p^T*lambda`. It treats the fixed-point
iterations as inactive and reuses the converged state plus local map
derivatives supplied by the caller.

The independent oracle uses the same nonlinear two-state map as the tangent
test. For `L(p) = u^T*x*(p)`, it completely re-solves positive and negative
perturbations of each parameter and compares the componentwise central
difference with the implicit VJP.

The analytical row uses a reusable converged fixed point and map derivatives.
The reference row performs four complete fixed-point solves for each
two-component VJP. Analytical uses 2,000,000 cotangents per sample; the slower
reference uses 25,000.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Kind | Median ns/VJP | MAD ns/VJP | Peak RSS |
|---|---|---:|---:|---:|
| analytical implicit solve | `analytical` | 132.3983 | 0.9601 | 12,550,144 B |
| four complete re-solves | reference oracle | 2,495.1440 | 8.1037 | 12,529,664 B |

Here “18.8457 times faster” means
`2495.1440 / 132.3983 = 18.8457`: it compares the analytical row with the
complete-re-solve reference for the same VJP, not with a standalone primal
iteration. Analytical saves 2,362.7457 ns/VJP and is the runtime selection; its
maximum observed self-process RSS is 20,480 bytes higher.

Peak memory is the maximum across five separately launched processes per
candidate, measured with `getrusage(RUSAGE_SELF)`. The machine-readable record
is `benchmark/reference/ryzen9_5950x_fixed_point_adjoint.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_fixed_point_adjoint analytical
taskset -c 4 fo exec bench_fixed_point_adjoint reference
fo exec bench_fixed_point_adjoint analytical --peak-rss
fo exec bench_fixed_point_adjoint reference --peak-rss
```

## Reusable preconditioner hook for implicit products

The multiroot analytical JVP, VJP, and scalar gradient accept an optional
preconditioned-solve callback plus caller-owned reusable context. The validation
applies one diagonal preconditioner across a JVP and VJP, checks that the same
context is invoked twice, and compares both products with finite differences of
the complete root solve.

The benchmark applies 200,000 JVP directions to a dense 16-by-16 residual
Jacobian. The hook candidate applies a reusable diagonal left preconditioner
before the direct solve; this is deliberately a simple reference implementation,
not a claim that diagonal preconditioning suits dense direct solves.

| Candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| default dense solve | `analytical` | 2,967.4902 | 38.2232 | 23,982,080 B |
| diagonal preconditioned hook | `analytical` | 5,338.7971 | 166.7145 | 23,887,872 B |

The preconditioned hook is 1.7991 times slower on this workload, so evidence
retains the default dense solve. The reusable hook is intended for later
iterative and application-specific candidates where iteration reduction can
repay its cost. Peak RSS again measures the complete `fo exec` process tree.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_multiroot_preconditioner.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_multiroot_preconditioner default
taskset -c 4 fo exec bench_multiroot_preconditioner preconditioned
```

## Build-time consumption of committed selections

The production build reads the selected Dawson outer-JVP candidate from its
committed tournament record with CMake's native JSON parser. CMake then emits a
build-local Fortran wrapper, so the selected call is resolved before compilation
and no registry dispatch remains in the numerical kernel.

An independent CMake fixture verifies parsing with a synthetic `hybrid`
selection. The existing complete-expression finite-difference test exercises
the generated production wrapper selected by the real record.

Clean configure measurements disable tests and examples and use a fresh Ninja
build directory for every sample:

| Configuration | Median ms | MAD ms | Peak RSS |
|---|---:|---:|---:|
| committed benchmark record | 680.8280 | 14.0784 | 36,429,824 B |
| hardcoded fallback | 678.9285 | 4.3068 | 36,442,112 B |

Reading and generating from the committed record adds 1.8995 ms, or 0.28%, to
median clean-configure time in this sample and does not increase measured peak
RSS. The committed-record path is selected because it connects benchmark
evidence to production without adding runtime dispatch.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_build_selection.json`.

## Static per-workload selection registry

CMake now generates a fixed registry from the four committed derivative
tournament records. Lookup keys are exact operator, derivative product, and
workload strings; values are the committed selected candidate IDs. An unknown
key returns no selection. Numerical wrappers continue to use generated direct
calls, so registry lookup is a pre-loop configuration operation.

The benchmark cycles over all four keys for 500,000 calls per sample:

| Candidate | Median ns/lookup | MAD ns/lookup | Peak RSS |
|---|---:|---:|---:|
| direct hardcoded ID | 24.1286 | 0.1279 | 23,859,200 B |
| static registry lookup | 211.5917 | 1.7550 | 23,547,904 B |

Static lookup adds 187.4631 ns. This cost is intentionally not paid in hot
loops; the registry exists to resolve or inspect a workload selection before
entering a kernel. Peak RSS measures the complete `fo exec` process tree.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_static_selection_registry.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_static_registry direct
taskset -c 4 fo exec bench_static_registry registry
```

## Conditioning diagnostics for implicit products

The multiroot analytical JVP, VJP, and scalar gradient expose an optional
reciprocal 1-norm condition estimate for the residual Jacobian. The estimate
solves against every coordinate vector and therefore remains opt-in. Validation
uses matrices whose exact reciprocal condition numbers are `1/3` and `1e-8`.

The 16-by-16 benchmark compares the unchanged product path with requesting the
diagnostic:

| Candidate | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---:|---:|---:|
| plain analytical product | 2,936.5041 | 9.2533 | 23,650,304 B |
| product plus reciprocal condition | 46,947.1965 | 928.5938 | 23,453,696 B |

Exact inverse-norm estimation is 15.9874 times as expensive for this workload.
The production default therefore remains the plain product; callers request the
diagnostic where derivative reliability matters. Peak RSS measures the complete
`fo exec` process tree.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_multiroot_condition.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_multiroot_condition plain
taskset -c 4 fo exec bench_multiroot_condition diagnostic
```

## Reliability status for ill-conditioned implicit products

The multiroot analytical JVP, VJP, and scalar gradient accept an optional
minimum reciprocal condition. If the measured condition falls below that
threshold, the product is zeroed and returns `FORTNUM_DOMAIN_ERROR` with an
unreliable-derivative message. A matrix with exact reciprocal condition `1e-8`
is independently verified to fail a `1e-6` threshold and pass `1e-9`.

Checking reliability performs the opt-in exact condition estimate:

| Candidate | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---:|---:|---:|
| plain analytical product | 2,927.8098 | 9.5219 | 23,977,984 B |
| product plus reliability threshold | 52,479.1808 | 242.1733 | 24,236,032 B |

Per-call reliability checking is 17.9244 times as expensive and raises measured
process-tree peak RSS by 258,048 bytes. It therefore remains explicit; callers
should request or cache it at solver and optimization boundaries rather than
inside repeated direction loops.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_multiroot_reliability.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_multiroot_condition plain
taskset -c 4 fo exec bench_multiroot_condition reliability
```

## Reliability status at the implicit JVP boundary

The callback-based `multiroot_implicit_jvp` now forwards an optional reciprocal
condition output and minimum acceptable reciprocal condition to its analytical
implicit solve. An independently known diagonal Jacobian with reciprocal
condition `1e-8` is reported within `1e-22`, rejected by a `1e-6` threshold,
and accepted by a `1e-9` threshold. Rejection zeroes the JVP and reports
`FORTNUM_DOMAIN_ERROR`.

The benchmark applies 10,000 directions per sample to a dense 16-state
callback residual. Both rows evaluate the same residual products and implicit
JVP. The reliability row additionally requests the exact inverse-norm-based
condition estimate and applies a passing threshold.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| plain implicit JVP boundary | `analytical` | 1,336.6043 | 8.3176 | 12,488,704 B |
| boundary plus reliability check | `analytical` | 16,699.6221 | 185.5896 | 12,537,856 B |

Here “12.4941 times slower” means
`16699.6221 / 1336.6043 = 12.4941`: reliability checking adds 15,363.0178 ns
to the same JVP. It also adds 49,152 bytes of maximum observed self-process
RSS. The plain product remains the repeated-direction selection; reliability
is an explicit boundary diagnostic that callers should evaluate or cache
outside hot loops.

Peak memory is the maximum across five separately launched candidates,
measured with `getrusage(RUSAGE_SELF)`. The machine-readable record is
`benchmark/reference/ryzen9_5950x_multiroot_implicit_jvp_reliability.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_multiroot_implicit_jvp_reliability plain
taskset -c 4 fo exec bench_multiroot_implicit_jvp_reliability reliability
fo exec bench_multiroot_implicit_jvp_reliability plain --peak-rss
fo exec bench_multiroot_implicit_jvp_reliability reliability --peak-rss
```

## Reliability status at the implicit VJP boundary

The callback-based `multiroot_implicit_vjp` now reports the same opt-in
reliability information as its JVP counterpart. An independently known
diagonal Jacobian with reciprocal condition `1e-8` is reported within `1e-22`,
rejected by a `1e-6` threshold, and accepted by a `1e-9` threshold. Rejection
zeroes the VJP and reports `FORTNUM_DOMAIN_ERROR`.

The benchmark applies 10,000 cotangents per sample to a dense 16-state
callback residual. Both rows evaluate the same state Jacobian, transposed
implicit solve, and parameter VJP. The reliability row additionally requests
the exact inverse-norm-based condition estimate and applies a passing
threshold.

Reference: AMD Ryzen 9 5950X, CPU 4 pinned, GNU Fortran 16.1.1, Release build,
15 samples after three warmups.

| Candidate | Mechanism | Median ns/VJP | MAD ns/VJP | Peak RSS |
|---|---|---:|---:|---:|
| plain implicit VJP boundary | `analytical` | 1,348.4467 | 2.4196 | 12,197,888 B |
| boundary plus reliability check | `analytical` | 16,751.9471 | 261.8053 | 12,275,712 B |

Here “12.4231 times slower” means
`16751.9471 / 1348.4467 = 12.4231`: reliability checking adds 15,403.5004 ns
to the same VJP. It also adds 77,824 bytes of maximum observed self-process
RSS. The plain product remains the repeated-cotangent selection; reliability
is an explicit boundary diagnostic that callers should evaluate or cache
outside hot loops.

Peak memory is the maximum across five separately launched candidates,
measured with `getrusage(RUSAGE_SELF)`. The machine-readable record is
`benchmark/reference/ryzen9_5950x_multiroot_implicit_vjp_reliability.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_multiroot_implicit_vjp_reliability plain
taskset -c 4 fo exec bench_multiroot_implicit_vjp_reliability reliability
fo exec bench_multiroot_implicit_vjp_reliability plain --peak-rss
fo exec bench_multiroot_implicit_vjp_reliability reliability --peak-rss
```

## Scalar-root candidate tournament

The scalar tournament compares both derivative directions for the same cubic
root `x*(p)` defined by `x^3 + p1*x - p2 = 0`:

- `analytical`: explicit residual partials plus the implicit root rule;
- `hybrid`: Enzyme residual products plus the analytical implicit root rule;
- `autodiff`: Enzyme differentiates a fixed 12-step Newton solve;
- finite-difference diagnostic: central differences of complete 12-step
  Newton solves.

All candidates are checked against central differences of independent
100-step bisection solves. Maximum absolute errors for the JVP were
`2.2072e-12`, `2.2072e-12`, `2.2071e-12`, and zero in table order. Maximum
VJP errors were `1.1891e-11`, `1.1891e-11`, `1.1891e-11`, and `1.4433e-11`.
The fixed Newton trace is intentionally smooth; these results do not make
autodiff through branch-dependent solver termination a general default.

The benchmark applies 100,000 varying directions or cotangents per sample.
Reference: AMD Ryzen 9 5950X, CPU 4 pinned, Flang/LLVM 22.1.8, Enzyme
`c96508349d9f`, Release `-O2`, 15 samples after three warmups.

| JVP candidate | Mechanism | Median ns/JVP | MAD ns/JVP | Peak RSS |
|---|---|---:|---:|---:|
| implicit with explicit residual products | `analytical` | 8.6395 | 0.2406 | 2,818,048 B |
| implicit with Enzyme residual products | `hybrid` | 8.0483 | 0.1378 | 2,797,568 B |
| Enzyme through 12 Newton steps | `autodiff` | 127.8436 | 0.4421 | 2,793,472 B |
| complete-solve central difference | diagnostic | 109.3232 | 0.5895 | 2,830,336 B |

The `hybrid` JVP is the measured winner. It is 1.0735 times faster than
`analytical`, 15.8845 times faster than `autodiff` through the iterations, and
13.5834 times faster than the finite-difference diagnostic. Its peak RSS is
20,480 bytes below `analytical`.

| VJP candidate | Mechanism | Median ns/VJP | MAD ns/VJP | Peak RSS |
|---|---|---:|---:|---:|
| implicit with explicit residual products | `analytical` | 17.3524 | 0.3155 | 2,830,336 B |
| implicit with Enzyme residual products | `hybrid` | 16.7265 | 0.0793 | 2,703,360 B |
| Enzyme through 12 Newton steps | `autodiff` | 177.6881 | 0.9734 | 2,711,552 B |
| componentwise complete-solve differences | diagnostic | 186.4905 | 0.9139 | 2,834,432 B |

The reverse-mode `hybrid` VJP is also the measured winner. It is 1.0374 times
faster than `analytical`, 10.6231 times faster than reverse `autodiff` through
the iterations, and 11.1494 times faster than the finite-difference
diagnostic. Its peak RSS is 126,976 bytes below `analytical`.

The machine-readable records are
`benchmark/reference/ryzen9_5950x_scalar_root_tournament_jvp.json` and
`benchmark/reference/ryzen9_5950x_scalar_root_tournament_vjp.json`.

Run validation and timing with:

```bash
cmake --build build-enzyme --target enzyme_scalar_root_hybrid_build \
    enzyme_scalar_root_vjp_hybrid_build
ctest --test-dir build-enzyme -R '^enzyme_scalar_root(_vjp)?_hybrid$' \
    --output-on-failure
taskset -c 4 \
    build-enzyme/test/ad/enzyme_scalar_root_hybrid.enzyme/enzyme_scalar_root_hybrid \
    --benchmark
taskset -c 4 \
    build-enzyme/test/ad/enzyme_scalar_root_vjp_hybrid.enzyme/enzyme_scalar_root_vjp_hybrid \
    --benchmark
```
