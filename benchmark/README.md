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
