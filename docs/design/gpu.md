# CPU and GPU differentiation contract

This document defines the supported differentiation paths and the minimum
contract for numerical leaves that may execute on a device. Public mechanism
names follow [ad.md](ad.md): `autodiff`, `analytical`, and `hybrid`.

## Support matrix

“Supported” means the path has a released-toolchain build, an independent
numerical oracle, proof that execution occurred on the requested target, and
committed complete-workload wall-clock and peak-memory evidence. A planned path
must not silently behave as a supported one.

| Target | `autodiff` | `analytical` | `hybrid` |
| --- | --- | --- | --- |
| Serial and parallel CPU | Supported selectively through the pinned Flang/Enzyme pipeline | Supported | Supported selectively as Enzyme composition with analytical custom rules |
| OpenACC GPU | Not supported | Planned, one generated leaf at a time | Not supported |
| OpenMP-target GPU | Not supported | Planned, one generated leaf at a time | Not supported |

CPU remains the complete feature path. Initial GPU work supports only explicit
`analytical` leaves that independently pass the device gates below. A missing
GPU derivative is reported as unavailable; it is never executed on the host
under a GPU label.

Direct Enzyme differentiation of GPU code, Enzyme inside OpenACC or OpenMP
offload, custom MLIR lowering, and a CUDA-specific derivative implementation
are outside the supported design. They may be reconsidered only after a
released upstream toolchain proves the complete path and a real consumer
justifies it.

## Device-leaf contract

A device leaf is the smallest mathematical value, JVP, or VJP kernel. The same
Fortran procedure body must be usable by serial CPU, parallel CPU, OpenACC, and
OpenMP target. A conforming leaf has:

- a `pure` procedure;
- scalar arguments or explicit-shape arrays with explicit, contiguous storage;
- fixed bounds known from arguments or compile-time parameters;
- no allocation or deallocation;
- no I/O, `stop`, or `error stop`;
- no mutable module or global state;
- no polymorphism, type-bound dispatch, or procedure pointers;
- no assumed-shape, allocatable, pointer, or other compiler-sensitive
  descriptors at the device boundary;
- status returned through a scalar argument when failure is possible; and
- stable floating-point semantics unless a separately validated relaxed
  contract is selected.

`fortsym` owns the explicit mathematical body and may emit only the leaf-level
annotations `!$omp declare target` and `!$acc routine seq`. It does not own
parallel scheduling, data movement, memory management, or runtime backend
selection.

`fortnum` owns thin batch wrappers. Those wrappers provide the parallel loop,
choose OpenACC or OpenMP target at configuration time, and manage application
data residency. Both backends call the identical numerical leaf and use the
same logical layout. Backend and derivative-candidate selection occur before
the launch loop.

## Acceptance gates

Each GPU leaf is unsupported until all applicable gates pass:

1. Byte-stable regeneration from the pinned `fortsym` revision.
2. Independent CPU numerical validation of value and derivative products.
3. GPU agreement with that oracle over normal, boundary, and failure cases.
4. Positive device-execution proof: OpenMP requires
   `omp_is_initial_device() == .false.`; OpenACC requires a released-runtime
   device check. Host fallback is failure.
5. Transfer-inclusive complete-workload wall clock and resident-data kernel
   wall clock measured separately.
6. Peak device allocation recorded, plus bandwidth, occupancy, registers, and
   spills where released tools expose reliable values.
7. Workload, layout, batch size, compiler, flags, device, driver, source
   revision, generator revision, and residency policy recorded.

The production candidate minimizes validated complete-workload wall clock,
subject to peak-memory constraints. Kernel-only timing, operation count, or a
mechanism name cannot select the winner.

## Immediate implementation order

The first pilot is the generated Dawson fused value/JVP leaf:

1. add annotation-only device emission to `fortsym`;
2. regenerate the same mathematical leaf with both optional annotations;
3. add `FORTNUM_GPU_BACKEND=NONE|OPENACC|OPENMP`, with unavailable selections
   failing configuration;
4. add one shared batch-wrapper template;
5. validate and measure OpenACC and OpenMP-target execution independently; and
6. retain only paths that pass the acceptance gates.

Coverage then expands by measured kernel family. No general portability
framework is introduced until at least two real kernels require the same
abstraction.
