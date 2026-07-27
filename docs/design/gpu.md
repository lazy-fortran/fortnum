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
| OpenACC GPU | Not supported | Dawson fused value/JVP pilot | Not supported |
| OpenMP-target GPU | Not supported | Dawson fused value/JVP pilot | Not supported |

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

## Shared batch-wrapper pilot

`benchmark/gpu/fortnum_gpu_batch_wrappers.f90` contains one Dawson fused
value/JVP loop body. OpenACC `parallel loop` and OpenMP
`target teams distribute parallel do` annotate that same loop; both call the
same generated numerical leaf. The wrapper owns scheduling only and contains
no derivative expression.

The CPU check establishes that this boundary does not add overhead before GPU
work begins. On an AMD Ryzen 9 5950X, GNU Fortran 16.1.1 `-O3`, CPU 4 pinned,
31 samples of 100 batches with 65,536 elements each:

| Candidate | Median ns/element | MAD ns | Peak RSS |
| --- | ---: | ---: | ---: |
| Shared wrapper | 9.1550 | 0.0573 | 12,578,816 B |
| Direct reference loop | 9.6499 | 0.0621 | 11,771,904 B |

The wrapper was 1.0541 times faster in this run; the evidence supports “no
measured CPU overhead,” not a general speedup claim. The independent oracle
evaluates the value and JVP formulas separately at 17 points. Machine-readable
evidence is in
`benchmark/reference/ryzen9_5950x_gpu_batch_wrapper_cpu.json`.

## OpenACC Dawson pilot

The generated analytical Dawson fused value/JVP is validated on a real NVIDIA
device by `test_openacc_dawson_offload`. A device-side
`acc_on_device(acc_device_nvidia)` check makes host fallback a test failure.
The independent CPU oracle evaluates
`sin(f) + f*f` and
`v*(1 - 2*x*f)*(cos(f) + 2*f)` at 4,096 points.

On an NVIDIA GeForce RTX 5060 Ti with driver 610.43.03, NVHPC 26.5
`nvfortran -fast -O3 -acc`, CPU 4 pinned, 31 samples of 1,048,576 elements:

| Complete workload | Median ms | MAD ms | Host peak RSS |
| --- | ---: | ---: | ---: |
| OpenACC, transfers included | 4.0772 | 0.1452 | 158,973,952 B |
| CPU generated-leaf loop | 8.1771 | 0.0571 | 47,403,008 B |

The transfer-inclusive OpenACC workload was 2.006 times faster for this batch.
This result includes three host-to-device arrays and two device-to-host arrays
on every sample. It does not measure resident-data kernel time or peak device
allocation; those belong to later checklist items. Machine-readable evidence
is in `benchmark/reference/rtx5060ti_openacc_dawson_transfer.json`.

## OpenMP-target Dawson pilot

The same generated analytical leaf, shared batch loop, flat array layout, and
1,048,576-element workload are compiled by NVHPC for OpenMP target. The
behavioral test runs with `OMP_TARGET_OFFLOAD=MANDATORY` and requires
`omp_is_initial_device()` to be false inside a target region. It then checks
the GPU value and JVP against the same independently evaluated CPU formulas at
4,096 points.

The released NVHPC compiler exposed that its device linker did not resolve the
previous module-list form of `declare target`. `fortsym` now emits the
procedure-local Fortran form documented by NVIDIA, at revision
`bba0d57c3278dce36abb991739e8128d1f3cfc53`. No numerical expression or
execution schedule was duplicated to fix it.

On the same NVIDIA GeForce RTX 5060 Ti with driver 610.43.03, NVHPC 26.5
`nvfortran -fast -O3 -mp=gpu -mp`, CPU 4 pinned, 31 samples:

| Complete workload | Median ms | MAD ms | Host peak RSS |
| --- | ---: | ---: | ---: |
| OpenMP target, transfers included | 4.0091 | 0.0741 | 165,863,424 B |
| CPU generated-leaf loop | 6.6490 | 0.0370 | 54,300,672 B |

The transfer-inclusive OpenMP-target workload was 1.6585 times faster for this
batch. As with the OpenACC result, resident-data kernel time and peak device
allocation are not part of this item. Machine-readable evidence is in
`benchmark/reference/rtx5060ti_openmp_dawson_transfer.json`.

## Persistent-data Dawson pilot

Both validated backends now keep the three input arrays and two output arrays
on the device across timed calls. Initial copies and the final output update
are outside the resident-data timing. Each synchronous timed call still
launches the shared 1,048,576-element batch wrapper. The behavioral tests run
both transfer-inclusive and resident paths and check both results against the
independent analytical CPU oracle.

On the same RTX 5060 Ti, NVHPC 26.5, CPU 4 pinned, with 31 samples:

| Backend and residency | Median ms | MAD ms | Host peak RSS |
| --- | ---: | ---: | ---: |
| OpenACC, transfers included | 3.9251 | 0.0839 | 159,031,296 B |
| OpenACC, resident | 0.2098 | 0.0003 | 144,572,416 B |
| OpenMP target, transfers included | 3.9800 | 0.0990 | 163,774,464 B |
| OpenMP target, resident | 0.2101 | 0.0003 | 149,323,776 B |

Persistent data reduces complete-call wall clock by 18.7088 times for OpenACC
and 18.9434 times for OpenMP target in this workload. The two resident paths
are within 0.14 percent of one another, so this measurement does not select a
backend. Peak device allocation is deliberately deferred to its dedicated
profiling checklist item. Machine-readable evidence is in
`benchmark/reference/rtx5060ti_dawson_residency.json`.

## Analytical VJP offload

`dawson_vjp_batch` adds scheduling around the generated contracted VJP leaf;
it contains no derivative expression. The OpenACC and OpenMP target directives
annotate the same loop and call the same `fortsym`-generated procedure. Both
transfer-inclusive and persistent-data paths run on the GPU.

The independent batched adjoint oracle checks

```text
sum(u * Jv) = sum(v * transpose(J)u)
```

at 4,096 nonuniform points after both residency paths. The JVP is also checked
against its independently evaluated analytical formula, so agreement is not
accepted merely because both products came from the same generator.

On the RTX 5060 Ti, NVHPC 26.5, CPU 4 pinned, with 31 samples of 1,048,576
elements:

| VJP workload | Median ms | MAD ms | Host peak RSS |
| --- | ---: | ---: | ---: |
| OpenACC, transfers included | 2.9421 | 0.0341 | 159,084,544 B |
| OpenACC, resident | 0.1690 | 0.0001 | 152,985,600 B |
| CPU leaf in OpenACC build | 6.7489 | 0.0313 | 47,685,632 B |
| OpenMP target, transfers included | 2.8081 | 0.0069 | 165,961,728 B |
| OpenMP target, resident | 0.1690 | 0.0001 | 157,728,768 B |
| CPU leaf in OpenMP build | 5.6961 | 0.0192 | 54,521,856 B |

Transfer-inclusive GPU execution is 2.2939 times faster than its CPU comparator
for OpenACC and 2.0285 times faster for OpenMP target. Persistent data improves
the GPU paths by another 17.4089 and 16.6160 times, respectively.
Machine-readable evidence is in
`benchmark/reference/rtx5060ti_dawson_vjp.json`.

## GPU fusion tournament

The value, contracted JVP, contracted VJP, fused value/JVP, and fused value/VJP
leaves all come from the same symbolic DAG. The tournament compares one fused
launch with two separate launches; it does not duplicate any expression.
Independent analytical and adjoint oracles validate every wrapper on both
devices.

Three batch regimes are measured:

- 256 elements: launch dominated;
- 65,536 elements: throughput transition; and
- 1,048,576 elements: sustained throughput when resident and
  transfer/memory dominated when copied for every call.

The tables show fused/separate median milliseconds, followed by the speedup
from fusion. Each cell is based on 31 samples after warmup.

| Backend | Elements | Residency | Value/JVP ms | Speedup | Value/VJP ms | Speedup |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| OpenACC | 256 | transfer | 0.0410 / 0.0480 | 1.1707x | 0.0441 / 0.0510 | 1.1565x |
| OpenACC | 256 | resident | 0.0081 / 0.0151 | 1.8642x | 0.0081 / 0.0150 | 1.8519x |
| OpenACC | 65,536 | transfer | 0.3419 / 0.3672 | 1.0740x | 0.3500 / 0.3579 | 1.0226x |
| OpenACC | 65,536 | resident | 0.0210 / 0.0329 | 1.5667x | 0.0210 / 0.0331 | 1.5762x |
| OpenACC | 1,048,576 | transfer | 3.9251 / 4.0579 | 1.0338x | 4.0021 / 4.0850 | 1.0207x |
| OpenACC | 1,048,576 | resident | 0.2100 / 0.3052 | 1.4533x | 0.2100 / 0.3052 | 1.4533x |
| OpenMP target | 256 | transfer | 0.0410 / 0.0480 | 1.1707x | 0.0410 / 0.0489 | 1.1927x |
| OpenMP target | 256 | resident | 0.0081 / 0.0159 | 1.9630x | 0.0081 / 0.0159 | 1.9630x |
| OpenMP target | 65,536 | transfer | 0.3391 / 0.3550 | 1.0469x | 0.3419 / 0.3519 | 1.0292x |
| OpenMP target | 65,536 | resident | 0.0210 / 0.0329 | 1.5667x | 0.0209 / 0.0332 | 1.5885x |
| OpenMP target | 1,048,576 | transfer | 4.0281 / 4.9010 | 1.2167x | 4.4570 / 4.7629 | 1.0686x |
| OpenMP target | 1,048,576 | resident | 0.2100 / 0.3119 | 1.4852x | 0.2101 / 0.3071 | 1.4617x |

At the largest batch, host peak RSS ranges from 151,797,760 to 174,792,704 B
across the 16 backend, residency, and fusion candidates; there is no
fusion-memory winner independent of backend and residency. Fused execution is
the wall-clock winner in every measured workload. Machine-readable medians,
MADs, and candidate-specific peak RSS are in
`benchmark/reference/rtx5060ti_dawson_fusion.json`.

## Immediate implementation order

The first pilot is the generated Dawson fused value/JVP leaf:

1. add annotation-only device emission to `fortsym`; complete at
   `349bc6257a22b416093624bd04dd1ed8a83852d0`;
2. regenerate the same mathematical leaf with both optional annotations;
   complete, with byte-stable regeneration;
3. add `FORTNUM_GPU_BACKEND=NONE|OPENACC|OPENMP`, with unavailable selections
   failing configuration; complete;
4. add one shared batch-wrapper template; complete for fused Dawson value/JVP;
5. validate and measure OpenACC and OpenMP-target execution independently;
   both are complete for the Dawson fused value/JVP pilot; and
6. retain only paths that pass the acceptance gates.

Coverage then expands by measured kernel family. No general portability
framework is introduced until at least two real kernels require the same
abstraction.

### Configuration

`FORTNUM_GPU_BACKEND` defaults to `NONE` and adds no accelerator flags.
`OPENACC` requires CMake to detect `OpenACC::OpenACC_Fortran` for the selected
released compiler. `OPENMP` additionally requires explicit
`FORTNUM_OPENMP_TARGET_FLAGS`; host-only OpenMP is rejected because it cannot
prove that a target region has a real offload target. The selected backend is
propagated through the `fortnum` target. This configuration gate establishes
compiler availability only; a backend remains unsupported until the
device-execution and performance gates above pass.
