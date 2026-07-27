# Differentiation benchmark evidence

Committed tables are reference measurements, not portable promises. Production
selection must be regenerated after a material change in hardware, compiler,
Enzyme, `fortsym`, primal source, or workload.

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
| refactor each JVP | `analytical` | 4,863.9215 | 82.1883 | 23,453,696 B |
| reuse primal LU | `analytical` | 813.4458 | 4.2406 | 23,977,984 B |

Reusing the factorization is 5.9794 times faster. Its measured peak RSS is
524,288 bytes higher. Peak RSS was measured with `/usr/bin/time -v` around the
`fo exec` process tree, so it includes the runner and is a conservative
end-to-end number, not candidate-private allocation.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_linear_solve_jvp_reuse.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_linear_solve_jvp refactor
taskset -c 4 fo exec bench_linear_solve_jvp reuse
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
| refactor transpose each VJP | `analytical` | 5,401.6115 | 55.3969 | 24,080,384 B |
| reuse transposed primal LU | `analytical` | 1,207.9997 | 3.6217 | 23,900,160 B |

Reusing the transpose factorization is 4.4715 times faster and its measured
peak RSS is 180,224 bytes lower. As for the JVP benchmark, peak RSS includes the
`fo exec` runner and is therefore an end-to-end process-tree measurement.

The machine-readable record is
`benchmark/reference/ryzen9_5950x_linear_solve_vjp_reuse.json`.

Run the benchmark with:

```bash
cd benchmark
fo build
taskset -c 4 fo exec bench_linear_solve_vjp refactor
taskset -c 4 fo exec bench_linear_solve_vjp reuse
```
