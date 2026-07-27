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
