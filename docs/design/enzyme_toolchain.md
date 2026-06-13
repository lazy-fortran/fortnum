# Flang + Enzyme autodiff toolchain

Status: accepted (issue #38, M6.2). Describes how fortnum builds and checks
Enzyme-produced derivatives, and how to activate the path on a machine that
has the Enzyme LLVM pass plugin.

fortnum is primal-first and stays fully buildable without Enzyme. The Enzyme
path is opt-in: `FORTNUM_ENABLE_ENZYME=OFF` by default, and even when ON a
missing plugin only skips the derivative tests. Nothing in `src/` depends on
Enzyme.

## Pipeline

Enzyme differentiates LLVM IR, so the route runs the Flang front end down to
IR, applies the Enzyme pass with `opt`, then compiles and links the result
with Flang again to pull in the Fortran runtime.

For a single source `kernel.f90` at optimization level `O2`:

```
flang-new -O2 -S -emit-llvm kernel.f90 -o kernel.ll
opt -load-pass-plugin=<ENZYME_PLUGIN>.so -passes=enzyme kernel.ll -S -o kernel.enzyme.ll
flang-new -O2 kernel.enzyme.ll -o kernel
```

Multiple sources each emit one IR file; `llvm-link -S a.ll b.ll -o linked.ll`
merges them before the `opt` stage so Enzyme sees the whole module. The kernel
modules come first, the source holding the `program` last.

`FortnumAddEnzymeTest.cmake` runs these stages from a generated CMake `-P`
driver. Any stage that returns nonzero fails the test build.

## Version assumptions

- Flang and LLVM `opt` must come from the same LLVM release. The plugin is
  built against one LLVM version and its name encodes the major:
  `LLVMEnzyme-<major>.so`.
- Developed and verified against flang-new, opt, and clang 22.1.6 with the
  Enzyme plugin for LLVM 22 (`LLVMEnzyme-22.so`).
- `FortnumEnzyme.cmake` reads the major version from `opt --version` and uses
  it to name the plugin it searches for, so the discovered plugin matches the
  `opt` that loads it.

## New vs legacy pass manager

LLVM 22 uses the new pass manager. The working invocation is
`-load-pass-plugin=<plugin> -passes=enzyme`.

Older Enzyme or LLVM builds expose the pass only through the legacy pass
manager. The fallback there is `-load=<plugin>` with `-enzyme` instead of
`-passes=`:

```
opt -load=<ENZYME_PLUGIN>.so -enzyme legacy_input.ll -S -o out.ll
```

fortnum targets the new-PM form. If a site pins an LLVM old enough to need the
legacy form, adjust the `opt` invocation in `FortnumAddEnzymeTest.cmake`; the
rest of the pipeline is unchanged.

## Conservative interface subset

The Enzyme path is tested against the simplest memory layouts first
(ad.md sec. 3). The smoke tests and the first generated entry points stay
inside this subset:

- scalar `real(real64)` arguments, passed by value;
- contiguous explicit-shape `real(real64)` arrays, with an integer length
  passed by value as an inactive size argument.

Descriptors, allocatables, assumed-shape arrays, polymorphism, optional active
arguments, and derived-type components are out of subset until a dedicated
Enzyme test covers that shape. No test, no support; the primal may still use
those features freely.

The differentiable kernel is an ordinary `bind(c)` Fortran function so its IR
symbol name is predictable and unmangled. The raw Enzyme entry points
(`__enzyme_autodiff` for reverse mode, `__enzyme_fwddiff` for forward mode) are
declared as `bind(c)` interfaces inside the test wrapper only. They never
appear in the public API; public callers use the `foo_jvp` / `foo_vjp` /
`foo_grad` / `foo_hvp` names (ad.md sec. 2).

### Argument activity in the wrapper

- Reverse mode: an active array `x` is paired with a shadow array `dx` that
  Enzyme accumulates the gradient into. The shadow is passed right after the
  primal array in the `__enzyme_autodiff` call. The output seed is the implicit
  unit, so the shadow holds the VJP `d out / d x`.
- Forward mode: an active array `x` is paired with a tangent array `dx`. The
  return value is the directional derivative (JVP) along the supplied tangents.
- A size or other inactive argument is forwarded by value with no shadow.

## Discovery and gating

`FortnumEnzyme.cmake` runs only when `FORTNUM_ENABLE_ENZYME=ON` and sets:

- `FORTNUM_FLANG_EXECUTABLE` (default search: `flang-new`, then `flang`),
- `FORTNUM_LLVM_OPT_EXECUTABLE` (`opt` and versioned names),
- `FORTNUM_LLVM_VERSION` (parsed from `opt --version`),
- `FORTNUM_ENZYME_PLUGIN` (searched under standard install prefixes),
- `FORTNUM_ENZYME_AVAILABLE` (TRUE only when all of the above resolve).

When the plugin is absent:

- `FORTNUM_ENZYME_REQUIRED=OFF` (default): the derivative tests register as
  skipped and the suite stays green.
- `FORTNUM_ENZYME_REQUIRED=ON`: configure fails with a hard error.

## Installing the plugin to activate the path

The plugin is not packaged on most distributions; build Enzyme against the
local LLVM 22 and either install it to a standard prefix or point fortnum at it.

Either install `LLVMEnzyme-22.so` to a searched prefix
(`/usr/lib`, `/usr/local/lib`, `/opt/enzyme/lib`), or set the environment
variable `ENZYME_PLUGIN_DIR` to its directory, or pass the path explicitly:

```
cmake -S . -B build -G Ninja \
    -DFORTNUM_ENABLE_ENZYME=ON \
    -DFORTNUM_ENZYME_PLUGIN=/path/to/LLVMEnzyme-22.so
```

Reconfigure and the derivative tests build and run instead of skipping. Verify
with `ctest -R enzyme`.
