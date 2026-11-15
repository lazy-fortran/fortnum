# Building fortnum

## Requirements

- a Fortran 2018 compiler (gfortran, ifx; `nvfortran` for OpenMP GPU builds)
- CMake 3.22 or newer, or fpm

CMake fetches the [fortgen](https://github.com/lazy-fortran/fortgen) runtime;
fpm takes it from `fpm.toml`.

## CMake

```bash
cmake -S . -B build -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
cmake --install build
```

| Option | Default | Effect |
| --- | --- | --- |
| `FORTNUM_BUILD_TESTING` | on for a top-level build | Build the test suite |
| `FORTNUM_BUILD_EXAMPLES` | on for a top-level build | Build the examples under `examples/` |
| `FORTNUM_CHECK_GENERATED` | off | Regenerate and compare the committed `fortsym` kernels |
| `FORTNUM_USE_BENCHMARK_SELECTION` | on | Read the committed derivative selections |
| `FORTNUM_ENABLE_ENZYME` | off | Build the supported CPU Enzyme tests |
| `FORTNUM_ENZYME_REQUIRED` | off | Fail if the requested Enzyme toolchain is unavailable |
| `FORTNUM_ENABLE_OPENMP` | on | Link host OpenMP for the Krylov vector kernels when available |
| `FORTNUM_GPU_BACKEND` | `NONE` | `NONE`, `OPENACC` or `OPENMP` |

Selecting an unavailable GPU backend fails configuration; GPU execution never
falls back silently to the host. OpenMP GPU builds currently require the NVIDIA
HPC SDK `nvfortran` compiler.

### Quad precision

CMake probes `selected_real_kind(33, 4931)` (`FORTNUM_HAVE_REAL128`). Without
it (for example with `nvfortran`) configuration warns, the modules
`fortnum_interval_qp`, `fortnum_cfft_rigorous` and `fortnum_fseries` are not
built, and tests with `real128` reference data are reported by CTest as
disabled. gfortran and ifx build the full library. fpm has no per-compiler source
selection and always builds everything.

## WebAssembly

fortnum and its test suite build for WebAssembly (wasm32) with LFortran and
emscripten through the CMake toolchain file `tools/wasm/emscripten-lfortran.cmake`;
`tools/wasm/test_wasm.sh [BUILD_DIR] [CTEST_ARGUMENT ...]` configures, builds and
runs the tests under node. It needs LFortran built with `WITH_TARGET_WASM`
(`LFORTRAN`, default `lfortran` on `PATH`), the emscripten SDK (`EMSDK_PATH`) and
node (`NODE`). LFortran's `main` lacks fixes this needs; the
[lazy-fortran fork](https://github.com/lazy-fortran/lfortran) carries them and
lists them in its `FORK.md`.
With `-DFORTNUM_WASM_ENVIRONMENT=shell` or `web` the test programs are single
files with the oracle data embedded, for a JavaScript shell such as
JavaScriptCore's `jsc` or for a browser page. WebAssembly has no OpenMP runtime
and runs on one thread; the OpenMP tests link the serial runtime functions of
`tools/wasm/omp_serial.c`.

## fpm

```bash
fpm build --profile release
fpm test --profile release --runner "python3 scripts/fpm_test_runner.py"
```

The runner (Python 3.11 or newer) passes the reference-data files listed under
`[extra.fo.test-args]` in `fpm.toml` to the tests that need them.

The contributor workflow, including the `fo` tool that runs static checks,
builds and tests, is described in [CONTRIBUTING.md](../CONTRIBUTING.md).

## Generated kernels

Some kernels are generated with [fortsym](https://github.com/lazy-fortran/fortsym)
and committed under `src/generated`. The tested fortsym revision is recorded in
`tools/codegen/fortsym.lock`. Set `FORTNUM_CODEGEN_VERBOSE=1` for diagnostic
generator output.
