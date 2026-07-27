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
| OpenACC GPU | Not supported | Validated generated explicit-kernel pilots | Not supported |
| OpenMP-target GPU | Not supported | Validated generated explicit-kernel pilots | Not supported |

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

## Multi-input scalar-output scaling

The next generated family uses one symbolic definition,

```text
f(x) = sum(sin(x_k)) + 0.5 * sum(x_k)^2,
```

with fixed 2-, 4-, 8-, and 16-input variants. `fortsym` emits fused
analytical value/JVP and value/VJP leaves; wrappers contain scheduling only.
The corrected DAG product generation preserves linear structure: post-CSE
operation counts are 16, 30, 58, and 114 for JVP and 15, 27, 51, and 99 for
VJP. The earlier repeated-partial construction would have produced 353 and
338 operations at 16 inputs.

Both real-device tests cover all four input counts. Independent CPU formulas
check the value, every JVP, and every VJP component at 1,024 nonuniform points;
the batched adjoint identity is a second oracle. OpenACC requires an NVIDIA
device, and OpenMP target requires `omp_is_initial_device()` to be false.

At 65,536 points with one product and resident data, 31-sample medians are:

| Product | Inputs | CPU ms | OpenACC ms | OpenMP ms | Fastest GPU/CPU |
| --- | ---: | ---: | ---: | ---: | ---: |
| JVP | 2 | 0.6559 | 0.0319 | 0.0310 | 21.16x |
| JVP | 4 | 1.3759 | 0.0551 | 0.0551 | 24.97x |
| JVP | 8 | 2.9981 | 0.1040 | 0.1049 | 28.83x |
| JVP | 16 | 6.8441 | 0.2150 | 0.2160 | 31.83x |
| VJP | 2 | 0.6759 | 0.0310 | 0.0310 | 21.80x |
| VJP | 4 | 1.3740 | 0.0542 | 0.0549 | 25.35x |
| VJP | 8 | 3.0999 | 0.1011 | 0.1020 | 30.66x |
| VJP | 16 | 7.7820 | 0.2072 | 0.2010 | 38.72x |

For eight inputs, increasing directions or cotangents from 1 to 4 to 16
increases wall clock approximately 4x and 16x on CPU and both GPUs. This
benchmark executes one contracted product per direction or cotangent; it
therefore measures the expected linear repeated-product cost, not a blocked
multi-product optimization.

Batch size determines the useful execution path:

| Product | Points | CPU ms | Best resident GPU ms | Best transfer GPU ms |
| --- | ---: | ---: | ---: | ---: |
| JVP | 256 | 0.0129 | 0.0141 | 0.0448 |
| JVP | 65,536 | 3.0201 | 0.1040 | 0.9389 |
| JVP | 1,048,576 | 48.5580 | 1.4620 | 14.3049 |
| VJP | 256 | 0.0129 | 0.0141 | 0.0449 |
| VJP | 65,536 | 3.0561 | 0.1011 | 1.0409 |
| VJP | 1,048,576 | 49.3989 | 1.4169 | 14.5468 |

Thus the CPU wins the launch-bound 256-point workload. With 65,536 or more
points, GPU wins even when transfers are included; resident execution wins by
29.0 to 34.9 times. OpenACC and OpenMP target differ by less than one percent
in every resident comparison, so this item does not select a backend.

Product-specific allocation avoids counting unused derivative buffers. Live
array storage is `8*n*(p + p*q + 2*q)` bytes for either product: 17,825,792 B
at `(p,n,q)=(16,65536,1)`, 88,080,384 B at `(8,65536,16)`, and 150,994,944 B
at `(8,1048576,1)`. Measured pinned-CPU peak RSS is respectively about 23.4,
93.4, and 156.3 MB. GPU-process host RSS includes compiler runtime and mapping
overhead and ranges from 129.7 to 260.9 MB in those regimes. Peak device
allocation and cache/bandwidth counters remain deliberately assigned to the
profiling checklist item.

Machine-readable medians, MADs, memory, exact toolchain, and validation
evidence are in
`benchmark/reference/rtx5060ti_multi_input_products.json`.

## Multi-input layout selection

The same generated eight-input leaves were benchmarked with the only two
admissible fixed layouts:

- `x(batch,active)`, where adjacent GPU threads access contiguous elements;
- `x(active,batch)`, where one thread's eight inputs are contiguous but
  adjacent threads access with stride eight.

Both layouts use the same leaf, arithmetic, launch schedule, allocation size,
and transfer volume. Real-device tests compare their value, JVP, and every VJP
component at 1,024 points. No layout abstraction or runtime dispatch was
introduced.

Resident 31-sample medians are:

| Product | Points | CPU batch/active ms | OpenACC batch/active ms | OpenMP batch/active ms |
| --- | ---: | ---: | ---: | ---: |
| JVP | 256 | 0.0129 / 0.0131 | 0.0141 / 0.0141 | 0.0141 / 0.0141 |
| JVP | 65,536 | 3.0458 / 3.0949 | 0.1040 / 0.1040 | 0.1049 / 0.1040 |
| JVP | 1,048,576 | 49.0649 / 49.7191 | 1.4620 / 1.4620 | 1.4620 / 1.4620 |
| VJP | 256 | 0.0131 / 0.0131 | 0.0141 / 0.0148 | 0.0140 / 0.0150 |
| VJP | 65,536 | 3.0690 / 3.0871 | 0.1011 / 0.1061 | 0.1011 / 0.1070 |
| VJP | 1,048,576 | 49.3751 / 49.2291 | 1.4181 / 1.4222 | 1.4172 / 1.4191 |

The layouts tie for resident GPU JVP within measurement resolution.
Batch-first is 4.9 percent faster for the 65,536-point OpenACC VJP and
5.8 percent faster for OpenMP target, and is 1.3 to 1.6 percent faster for the
larger CPU JVPs. Large transfer-inclusive measurements varied by more than the
layout effect and did not consistently favor either ordering even over 155
samples. Since transfer byte counts are identical, that noise is not used to
select the numerical layout.

At 1,048,576 points, both layouts require the same 150,994,944 live array
bytes. CPU peak RSS is 156.2 to 156.3 MB. GPU-process host peak RSS ranges from
196.8 to 260.6 MB with no consistent layout advantage.

`x(batch,active)` is selected: it provides the expected coalesced mapping, is
never materially slower while resident, wins the mid-size VJP, and is slightly
better for CPU JVP. The alternate wrappers remain benchmark-only evidence;
they are not a public layout framework. Full medians, repeated transfer
measurements, memory, validation, and toolchain provenance are in
`benchmark/reference/rtx5060ti_multi_input_layout.json`.

## Released-tool GPU profile

The selected batch-first layout was profiled at 1,048,576 points with released
Nsight Compute 2026.2.1 and Nsight Systems 2026.1.3. Nsight Compute collected
the `SpeedOfLight`, `LaunchStats`, `Occupancy`, and
`MemoryWorkloadAnalysis` sections for one launch after three skipped launches.
Nsight Systems traced CUDA allocation events; peak device allocation is the
maximum running allocation-minus-free total. Performance counters were read
through passwordless `sudo` because this host disables unprivileged access.

| Backend | Product | Duration ms | GB/s | Compute % | DRAM % | Registers/thread | Achieved occupancy | Spills | Peak device allocation |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| OpenACC | JVP | 1.46 | 102.42 | 84.66 | 23.21 | 147 | 24.41% | 0 | 153,051,422 B |
| OpenACC | VJP | 1.41 | 99.37 | 84.68 | 22.52 | 150 | 24.31% | 0 | 153,051,422 B |
| OpenMP target | JVP | 1.46 | 102.47 | 84.67 | 23.23 | 147 | 24.52% | 0 | 153,125,080 B |
| OpenMP target | VJP | 1.41 | 99.65 | 84.66 | 22.58 | 150 | 24.55% | 0 | 153,125,080 B |

Theoretical occupancy is 25 percent and is register-limited. No local-memory
spilling request occurs. L1 hit rate is 52.45 to 52.72 percent and L2 hit rate
is 43.66 to 43.95 percent. Compute throughput is about 84.7 percent while DRAM
throughput is only 22.5 to 23.2 percent, so reducing bytes alone is not the
leading optimization for this transcendental-heavy kernel. The near-identical
profiles also corroborate that the two scheduling backends reach equivalent
generated device code.

The live numerical arrays account for 150,994,944 B; the roughly 2.1 MB
difference to peak device allocation is runtime allocation overhead. Complete
profile metadata and the extraction definition are in
`benchmark/reference/rtx5060ti_multi_input_profile.json`.

## Reproducible GPU figures

`benchmark/report/app/plot_gpu_report.f90` generates six focused figures from
the normalized CSV copies of the committed benchmark evidence:

- JVP and VJP complete wall clock versus batch size on log-log axes, including
  CPU, resident GPU, and transfer-inclusive GPU paths;
- JVP and VJP wall clock versus active-input count on log axes;
- achieved device-memory bandwidth; and
- peak device allocation.

The runtime figures use an Okabe-Ito-derived color-safe palette plus distinct
line styles and markers. Backend and product names label every profile bar, so
color is not required to decode them. Axes state quantities, units, and log
scales explicitly. The rendered color figures and a grayscale conversion were
inspected for labels, overlap, and redundant series distinction.

Reproduce them from the repository root:

```bash
mkdir -p /tmp/fortnum-gpu-report
cd benchmark/report
fo exec plot_gpu_report \
  data/gpu_batch_scaling.csv \
  data/gpu_active_scaling.csv \
  data/gpu_profile.csv \
  /tmp/fortnum-gpu-report
```

The repository commits the generator and data only. `benchmark/report` ignores
its output directory and all PNG files.

## Benchmark provenance

Every committed RTX 5060 Ti benchmark record now carries the same auditable
provenance fields: exact fortnum source revision, residency policy,
device-execution proof, GPU and driver, compiler and flags, and exact `fortsym`
revision. Historical records point to the commit that introduced the measured
executable and evidence; the profile record points to the selected-layout
source it measured. These identifiers describe the measured state rather than
the later documentation-only commit that added normalized provenance keys.

## Backend selection

OpenACC and OpenMP target remain separately compiled candidates; mixing their
compiler modes into one runtime binary would add complexity without improving
the numerical kernel. The benchmark-level
`select_multi_input_gpu_backend` lookup therefore resolves the backend before
the application chooses an executable and before any launch loop.

The lookup contains only the 12 measured combinations of JVP/VJP,
transfer/resident data, and 256/65,536/1,048,576 points. An unknown key returns
unavailable; it is never rounded to a nearby size and never silently defaults.
Both candidates have passed the same independent formulas and real-device
tests.

All measured backend wall clocks are tied within the existing three-percent
timing threshold after the repeated large-transfer sample. The deterministic
secondary comparison selects OpenACC because its profiled peak device
allocation is 153,051,422 B versus 153,125,080 B for OpenMP target. Thus the
winner follows measured time and memory, not the backend name or the
`analytical` mechanism shared by both candidates.

The exact workload table and selection rule are in
`benchmark/reference/rtx5060ti_gpu_backend_selection.json`. Real-device tests
exercise every measured lookup and prove that an unmeasured size is rejected.

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
