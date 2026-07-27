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
