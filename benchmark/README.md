# fortnum benchmarks

Micro-benchmark harness for primal call overhead. Each benchmark is a
`(name, callable)` pair; the harness times the callable over many
repetitions and reports nanoseconds per call.

## Run

The benchmark build is standalone and never touches `src/` or `test/`:

```bash
cmake -S benchmark -B build-bench -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-bench -j
./build-bench/bin/bench_main
```

Use `Release`. At `-O0` the numbers measure the compiler's missing
optimizer, not the routine.

Output is one row per case: name, repetitions, ns/call, and the
derivative columns (empty until M6 fills them).

## Add a benchmark

A kernel is a function returning `real(dp)`. Return the result so the
optimizer cannot delete the work:

```fortran
function kernel_my_routine() result(sink)
    use my_module, only: my_routine
    real(dp) :: sink
    sink = my_routine(...)
end function kernel_my_routine
```

Register it in `bench_main.f90`:

```fortran
call registry%add("my_routine", kernel_my_routine)
```

If the routine pulls in more `src/` modules, add their source paths to
`add_executable(bench_main ...)` in `benchmark/CMakeLists.txt`.

## Derivative overhead (M6, issue #36)

`registry%add` takes an optional `deriv=` kernel of the same shape. Pass
the JVP/VJP variant of a routine and the harness times both, then reports
the primal/derivative ratio:

```fortran
call registry%add("rosenbrock", primal_rosenbrock, deriv=vjp_rosenbrock)
```

The result type and the report table already carry the derivative
columns, so M6 adds autodiff kernels without changing the harness.

## Regression gate

`bench_main --json` emits a machine-readable run that `gate.py` compares
against `baseline.json`. A primal benchmark slower than `baseline * factor`
(default 2.0) fails the gate with a nonzero exit. CI runs it on every push.

```bash
./build-bench/bin/bench_main --json | python3 benchmark/gate.py
# or against a saved run, with a custom factor:
python3 benchmark/gate.py --run run.json --factor 1.5
```

The JSON schema is one object per benchmark under `benchmarks`:
`name`, `reps`, `ns_per_call`, a `backend` tag
(`analytic | implicit | trace | generated | primal`), and the
derivative-product fields `deriv_ns_per_call`, `jvp_primal`, `vjp_primal`,
`grad_primal`, `hvp_primal`. The derivative fields are `null` until M6
supplies autodiff kernels; the schema reserves them now so the baseline
format does not change when they arrive.

Derivative-overhead ratios are checked but non-blocking by default (runner
noise on those ratios is not yet characterized). Pass `--gate-derivative`
to make them fail the gate once a noise budget is set.

### Refresh the baseline

Regenerate after an intended performance change, on a quiet machine, in a
Release build:

```bash
cmake -S benchmark -B build-bench -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-bench -j
./build-bench/bin/bench_main --json > benchmark/baseline.json
```

Commit the new `baseline.json` with the change that motivated it. The same
command refreshes derivative timings: once `bench_main` registers JVP/VJP
kernels, their `deriv_ns_per_call` and `*_primal` ratios populate
automatically in the regenerated baseline.
