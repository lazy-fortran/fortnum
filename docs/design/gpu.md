# CPU and GPU differentiation contract

Status: analytical GPU pilots implemented on released NVHPC.

## Support matrix

| Target | `autodiff` | `analytical` | `hybrid` |
| --- | --- | --- | --- |
| CPU | supported for tested Enzyme shapes | supported | supported for tested compositions |
| GPU | unsupported | supported for listed device leaves | unsupported |

The GPU path does not claim parity with the CPU derivative surface. A requested
unavailable backend fails. Device tests prove execution on the accelerator, so
host fallback cannot satisfy them.

## Architecture

GPU support separates three concerns:

1. `fortsym` emits a pure numerical leaf and optional procedure annotations.
2. `fortnum` owns a small batch loop and backend scheduling directives.
3. the application owns data lifetime and persistent residency.

OpenACC and OpenMP target call the same numerical leaf. The mathematical body
contains no backend branch.

This is the useful part of a Kokkos-like design without a container library,
execution-space framework, task graph, memory manager, or runtime scheduler.

## Device-leaf contract

A device leaf has:

- a pure procedure
- explicit scalar or contiguous fixed-shape data
- no allocation or deallocation
- no I/O
- no mutable global state
- no polymorphic arguments
- no procedure pointers
- no unsupported Fortran descriptors

`fortsym` may add OpenACC `routine seq` and OpenMP `declare target`
annotations. It does not generate launch schedules or data regions.

## Configure

`FORTNUM_GPU_BACKEND` accepts:

- `NONE`
- `OPENACC`
- `OPENMP`

The default is `NONE`. OpenMP target also requires explicit
`FORTNUM_OPENMP_TARGET_FLAGS`, because host-only OpenMP is not a GPU backend.

Released NVHPC supplies both validated NVIDIA paths. OpenMP target is accepted
only with `nvfortran`; other compilers fail configuration until they pass the
same gates on released software.

## Acceptance gates

Every supported leaf and backend requires:

1. byte-stable regeneration from the lock file
2. an independent CPU mathematical oracle
3. mandatory non-host device execution
4. validation at representative inputs and boundaries
5. transfer-inclusive complete-workload wall clock
6. resident-data wall clock
7. peak host and device memory
8. compiler, flags, driver, hardware, source, and generator provenance

Profiling records bandwidth, occupancy, register use, and spills when released
tools expose reliable values.

## Implemented coverage

| Family | Analytical GPU products | Independent oracle |
| --- | --- | --- |
| Dawson | value, JVP, VJP, fused products | CPU formula and adjoint identity |
| multi-input scalar | JVP/VJP scaling at 2, 4, 8, 16 inputs | generated CPU products |
| fixed-cell Lagrange | value, JVP, VJP | cubic formula and adjoint identity |
| fixed quadrature | arbitrary-order JVP/VJP batch loops | exact polynomials and adjoint identity |
| fixed-size linear algebra | 3 by 3 determinant and inverse JVP/VJP | finite difference, matrix identity, adjoint identity |
| FFT | length-8 complex JVP/VJP | direct DFT and complex adjoint identity |
| scalar root residual | local products plus implicit JVP | closed-form root |
| ODE trace | two-state tangent and adjoint maps | closed-form matrix composition |
| terminal ODE objective | value plus reverse gradient | closed-form terminal state and gradient |

These are individually validated leaves. Unsupported dimensions and algorithms
do not silently route to a different GPU implementation.

## Data layout and residency

The measured multi-input layout is `x(batch, active)`. The batch index is the
leftmost Fortran index, so adjacent GPU threads access adjacent elements.

Resident and transfer-inclusive measurements answer different questions:

- resident time selects repeated kernels inside an existing device data region
- transfer-inclusive time selects a complete host-call workload

Applications should keep reusable inputs, traces, and outputs resident across
composed kernels. The persistent terminal-objective benchmark measures forward
state propagation, objective construction, reverse propagation, and gradient
return as one application.

## Fusion and algebraic variants

Generated value/JVP and value/VJP fusion avoids duplicate primal work and extra
launches. `fortsym` can emit raw, simplified, and factored variants from one
symbolic DAG. Symbolic equivalence and numerical boundary tests run before
native timing.

Post-CSE operation count filters candidates. Native complete-workload wall
clock decides the production form. Smaller source or fewer operations do not
override slower execution.

Only the selected production leaf is committed. Losing variants are
reproducible temporary outputs.

## Backend selection

OpenACC and OpenMP target are separate compiled executables. A workload lookup
selects the backend before launch loops from:

- product
- batch size
- residency
- validated timing
- peak device allocation as a deterministic tie-breaker

Unknown keys return unavailable. Current measured backend differences for the
multi-input pilot are practical ties, with the deterministic memory tie-break
selecting OpenACC.

## Performance interpretation

GPU speedup depends on batch size, arithmetic intensity, returned arrays, and
residency. Small launch-bound workloads often remain on CPU. Large resident
workloads can favor GPU strongly. Returning many VJP arrays can make transfer
cost dominate even when the resident kernel is fast.

Exact measurements live in `benchmark/reference/rtx5060ti_*.json`. Aggregate
plot data and `fortplot` source live under `benchmark/report`.

Generate the GPU figures:

```bash
output_dir=$(mktemp -d)
cd benchmark/report
fo exec plot_gpu_report \
  data/gpu_batch_scaling.csv \
  data/gpu_active_scaling.csv \
  data/gpu_profile.csv \
  "${output_dir}"
```

Generated PNG files stay outside the repository.

## Exclusions

The supported design excludes:

- Enzyme differentiation of GPU kernels or offload regions
- Enzyme integration into `nvfortran`
- custom MLIR or LLVM lowering
- CUDA-specific derivative source
- a replacement for Kokkos
- a CPU Enzyme callback that launches an analytical GPU rule

Reconsider an exclusion only when a released toolchain supports the complete
path and a measured application needs it.
