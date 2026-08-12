# CPU microarchitecture variants and runtime dispatch

Status: investigated; mechanism deferred pending a distributed-binary consumer.

## Question

For CPU targets, "one Fortran kernel" is not automatically the end of the
story: the same source compiled for different microarchitectures differs in
vector width, FMA availability, and instruction selection (`-march`).  This
document records (1) whether any supported Fortran compiler exposes function
multiversioning and (2) whether compiling generated kernels at
`-march=native` leaves a measured gap that a runtime-dispatch mechanism would
recover.

## Compiler multiversioning availability

GCC's function multiversioning (`__attribute__((target_clones(...)))` and the
resulting ifunc resolver) is a C/C++ frontend feature.  The GNU Fortran
frontend exposes no equivalent:

- `target_clones` and `target` are not accepted in `!GCC$ attributes`
  statements; gfortran rejects them with `Unknown attribute`.
- gfortran's `-march`/`-mtune` options build one fixed target per object;
  there is no per-function multiversioning directive.

Neither ifx nor flang exposes a Fortran function-multiversioning directive
that would let one source emit per-microarchitecture clones and an ifunc
dispatcher.  The observation that drives the generated-kernel approach stands:
because `fortsym` emits the same expression N times under different module
names, compiling each with a different `-march` and dispatching at runtime is
the only route on the supported compilers, and it is nearly free relative to
hand-written per-architecture assembly.

## The cheap measurement before building

Before building per-microarchitecture variants and a runtime selector, the
relevant question is whether a baseline `-march=x86-64-v2` build leaves a gap
that `-march=native` recovers for the existing generated kernels.  The two
existing kernels bracket the expected behavior:

- `fortnum_dawson_outer` fused JVP: transcendental-heavy (sin/cos), scalar.
- `fortnum_det3_jvp` batch: pure arithmetic, FMA-able, autovectorizable.

The evidence is collected by `benchmark/collect_march_gap.py`, which builds
each kernel at `-O3 -march=x86-64-v2` and `-O3 -march=native`, interleaves
runs to cancel run-order bias on shared hosts, and records pooled medians
plus a verdict against a 5% reproducible-noise margin.  The committed record
is `benchmark/reference/xeon_e5_2630v4_march_gap.json`.

On the Intel Xeon E5-2630 v4 (Broadwell, AVX2/FMA) reference host the native
to baseline ratio stayed within 0.97--1.04 for both kernels across repeated
runs: `-march=native` never led the baseline by more than the noise margin.
The scalar transcendental kernel is dominated by `sin`/`cos` and does not
expose a microarchitecture gap; the batched arithmetic kernel is
memory/compute bound at the batch sizes exercised and shows the same result.

## Decision

Defer the runtime microarchitecture-dispatch mechanism.  A single-machine
`-march=native` build already matches the baseline `x86-64-v2` build for the
measured generated kernels, so there is no benefit to emit per-microarchitecture
variants and dispatch at runtime within this repository.

Revisit only when a concrete consumer needs a distributed binary that must
match native performance on heterogeneous machines.  At that point, because
the kernel is generated, the OpenBLAS-style path is cheap: emit the same
expression under N names, compile each with its own `-march`, and add a
runtime selector keyed on detected CPU features.  Prefer a small number of
microarchitecture levels (for example `x86-64-v2` and `x86-64-v3`) over a
long per-model list, and keep blocking and layout parameters derived from an
analytical model (cache sizes, register count, vector width) rather than an
ATLAS-style empirical search, consistent with the BLIS route.
