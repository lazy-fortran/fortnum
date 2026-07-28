# CPU Flang and Enzyme toolchain

Status: supported for tested CPU candidate shapes.

Enzyme is optional. Normal library builds contain no runtime Enzyme dependency.

## Pipeline

The CMake helper performs:

```text
Fortran source
  -> Flang LLVM IR
  -> llvm-link for multi-source fixtures
  -> opt with the Enzyme pass plugin
  -> Flang link with the Fortran runtime
```

Flang, `opt`, `llvm-link`, and the Enzyme plugin must target the same LLVM
major release. The repository tests the new pass-manager form:

```bash
opt -load-pass-plugin=/path/to/LLVMEnzyme-MAJOR.so \
    -passes=enzyme input.ll -S -o output.ll
```

## Configure

```bash
cmake -S . -B build-enzyme -G Ninja \
  -DFORTNUM_ENABLE_ENZYME=ON \
  -DFORTNUM_ENZYME_PLUGIN=/path/to/LLVMEnzyme-MAJOR.so
cmake --build build-enzyme
ctest --test-dir build-enzyme -R enzyme --output-on-failure
```

Discovery sets the Flang, LLVM utility, version, plugin, and availability cache
variables. With `FORTNUM_ENZYME_REQUIRED=OFF`, an unavailable plugin skips the
optional path. With `FORTNUM_ENZYME_REQUIRED=ON`, configuration fails.

## Tested ABI profiles

Generated wrappers cover:

- one to five active scalar `real(real64)` inputs
- scalar output with forward JVP and reverse VJP
- one scalar plus a fixed-size vector
- fixed-size arrays with inactive integer controls
- optional analytical forward custom rules

Active arrays are contiguous and explicit-shape. Size and mode arguments pass
by value as inactive controls.

Assumed-shape descriptors, allocatable active values, polymorphic active
values, optional active arguments, procedure pointers, and arbitrary derived
types are unsupported until a dedicated real-Enzyme test establishes their
ABI.

## Generated boundaries

Raw `__enzyme_fwddiff` and `__enzyme_autodiff` declarations live only in
temporary generated wrappers or isolated ABI smoke probes. Hybrid fixtures use
the shared wrapper generator and shared timing, statistics, memory, and
custom-rule counter support.

The repository guard rejects new raw intrinsics or duplicate fixture helpers
under `cmake/enzyme/hybrid`.

## Operator rules

Enzyme differentiates local smooth kernels. Solvers and adaptive algorithms
remain outside the pass:

- root fixtures differentiate residual components or fixed-iteration
  comparators, then use analytical implicit boundaries
- integration fixtures differentiate integrands or a frozen trace
- ODE fixtures differentiate a local RHS inside an analytical sensitivity
  recurrence
- direct-solver fixtures compare Enzyme with factor-reusing implicit products

CPU Enzyme does not launch GPU kernels and does not differentiate OpenACC or
OpenMP target regions.

## Mixed-mode second order

`cmake/enzyme/hvp/forward_over_reverse.f90` proves nested mixed mode through
the released pipeline. Its inner reverse call generates a nonlinear
least-squares gradient; an outer forward call differentiates that gradient
along one direction. The transformed HVP is checked against a closed-form
Hessian contraction.

`cmake/enzyme/hvp/reverse_over_forward.f90` applies the complementary order:
an inner forward call forms a scalar directional derivative and an outer
reverse call differentiates it with respect to the primal inputs. It is
checked against the same contraction.
