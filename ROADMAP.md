# Differentiation implementation plan

This file tracks the remaining work for differentiable programming in
`fortnum`. Public terminology is defined in `docs/design/ad.md`: `autodiff`,
`analytical`, and `hybrid`.

## Hard execution rules

Work through this file one checkbox at a time.

1. Select exactly one unchecked item.
2. Keep the implementation as small as possible while fully satisfying that
   item.
3. Avoid overengineering. Prefer simple, direct, effective code over general
   frameworks, speculative abstractions, and unused flexibility.
4. Add an independent behavioral oracle and run the relevant test matrices.
5. Measure and report performance for every item. Report at least runtime and
   peak memory, with before and after numbers or competing-candidate numbers.
   Record the workload, compiler, hardware, and measurement method. Wall-clock
   time for the complete workload is the primary selection metric.
6. Measure scaling with active-input count, output count, and derivative
   directions wherever those dimensions can change the forward-versus-reverse
   verdict. Report CPU work, peak memory, and cache behavior in addition to
   wall-clock time when the host exposes reliable counters.
7. Update committed benchmark evidence when the item affects differentiation
   performance or candidate selection.
8. Check off the completed item in this file.
9. Commit the implementation, tests, evidence, and checklist update together.
10. Push the commit to `main`.
11. Only after the push succeeds, select the next unchecked item.
12. Keep one mathematical source of truth for value, JVP, VJP, and fused
    products. Generate mechanical derivatives and wrappers instead of manually
    duplicating expressions.
13. Keep backend selection and derivative-candidate selection outside hot
    loops. OpenMP target and OpenACC must execute the same numerical leaf.
14. Treat silent GPU-to-host fallback as a failure. Transfer-inclusive
    complete-workload wall clock is the primary GPU metric; also report
    resident-kernel time, peak device memory, bandwidth, occupancy, registers,
    and spills where reliable tools expose them.
15. Preserve stable floating-point semantics by default. Reassociation,
    approximate math, and relaxed contraction require an explicit numerical
    contract, an independent accuracy oracle, and a complete-workload win.
16. Do not build speculative portability frameworks. Add a shared abstraction
    only after at least two real kernels need the same behavior.

Do not combine several checklist items into one implementation cycle. If an item
turns out to contain several independent changes, split it into smaller
checkboxes before writing code.

A mechanism name never selects a winner. Production selection requires
independent validation plus measured application runtime and peak memory.

The initial target matrix is deliberately asymmetric:

| Target | `autodiff` | `analytical` | `hybrid` |
| --- | --- | --- | --- |
| CPU | Enzyme forward and reverse | generated or hand-written | Enzyme composition with analytical custom rules |
| GPU | not initially supported | supported incrementally | not initially supported |

An analytical residual product followed by an analytical implicit tangent or
adjoint solve remains `analytical`. The public term `hybrid` requires autodiff
composition with at least one analytical rule.

## Completed foundation

- [x] Define `autodiff`, `analytical`, and `hybrid` as candidate mechanisms.
- [x] Implement derivative provenance and quality merging.
- [x] Implement deterministic candidate metadata and pre-loop selection.
- [x] Build a real Flang and Enzyme custom-rule path.
- [x] Generate and consume a fused Dawson value/JVP kernel with `fortsym`.
- [x] Benchmark analytical, autodiff, and hybrid Dawson candidates.
- [x] Implement analytical implicit linear-solve JVP and VJP products.
- [x] Implement combined-active Lagrange and B-spline JVP/VJP products.

## Core infrastructure

- [x] Add factorization reuse to the analytical linear-solve JVP.
- [x] Add transpose-factorization reuse to the analytical linear-solve VJP.
- [x] Add reusable preconditioner hooks for implicit products.
- [x] Connect committed benchmark records to build-time candidate selection.
- [x] Generate a static per-workload selection registry.
- [x] Add conditioning diagnostics to implicit derivative products.
- [x] Add derivative-reliability status for ill-conditioned implicit products.
- [x] Measure candidate-specific peak memory instead of shared-process RSS.
- [x] Audit the linear-solve VJP benchmark target and refresh its committed
  Release and candidate-memory evidence.

## Roots and residual equations

- [x] Add a generic analytical implicit tangent boundary for scalar roots.
- [x] Add a generic analytical implicit adjoint boundary for scalar roots.
- [x] Add a generic analytical implicit tangent boundary for vector roots.
- [x] Add a generic analytical implicit adjoint boundary for vector roots.
- [x] Add analytical implicit tangent products for fixed points.
- [x] Add analytical implicit adjoint products for fixed points.
- [x] Add a hybrid scalar-root JVP using forward-mode autodiff residual products.
- [x] Add a hybrid scalar-root VJP using reverse-mode autodiff residual products.
- [x] Add a hybrid vector-root JVP using forward-mode autodiff residual products.
- [x] Add a hybrid vector-root VJP using reverse-mode autodiff residual products.
- [x] Return and reuse converged Jacobians from analytical vector-root solves.
- [x] Reuse converged root Jacobian factorizations for JVPs.
- [x] Reuse converged transposed root Jacobian factorizations for VJPs.
- [x] Report JVP reliability near singular root Jacobians.
- [x] Report VJP reliability near singular root Jacobians.
- [x] Benchmark analytical implicit, hybrid, autodiff-through-iterations, and
  finite-difference scalar-root candidates.
- [x] Benchmark analytical implicit, hybrid, autodiff-through-iterations, and
  finite-difference vector-root candidates.

## Integration and quadrature

- [x] Add analytical differentiation under a fixed-bound integral.
- [x] Add analytical moving-lower-bound terms.
- [x] Add analytical moving-upper-bound terms.
- [x] Compose autodiff integrand JVPs with analytical fixed quadrature.
- [x] Compose autodiff integrand VJPs with analytical fixed quadrature.
- [x] Add an analytical frozen-trace adaptive-integration candidate.
- [x] Add an autodiff frozen-trace adaptive-integration candidate.
- [x] Add a hybrid adaptive-integration candidate.
- [x] Benchmark fixed quadrature candidates.
- [x] Benchmark smooth adaptive-integration candidates.
- [x] Benchmark singular adaptive-integration candidates.
- [x] Benchmark batched integration candidates.

## ODE solvers and events

- [x] Add hybrid forward sensitivities using autodiff RHS JVPs.
- [x] Add an analytical discrete adjoint for one explicit ODE method.
- [x] Accumulate parameter VJPs in the discrete adjoint.
- [x] Add checkpointed reverse differentiation as a candidate.
- [x] Add recomputation-based reverse differentiation as a candidate.
- [x] Add analytical implicit-stage tangent products.
- [x] Add analytical implicit-stage adjoint products.
- [x] Differentiate event times through the event residual equation.
- [x] Document and test the continuous sensitivity contract.
- [x] Document and test the discrete sensitivity contract.
- [x] Benchmark a short nonstiff trajectory.
- [x] Benchmark a long nonstiff trajectory.
- [x] Benchmark a stiff trajectory.
- [x] Benchmark a many-parameter trajectory.
- [x] Benchmark an event-driven trajectory.

## Linear algebra

- [x] Add an analytical determinant JVP.
- [x] Add an analytical determinant VJP.
- [x] Add an analytical inverse JVP for callers that require an inverse.
- [x] Add an analytical inverse VJP for callers that require an inverse.
- [x] Add a reusable LU factorization object.
- [x] Add multiple-right-hand-side tangent solves.
- [x] Add multiple-right-hand-side adjoint solves.
- [x] Add a forward-mode autodiff direct-solver JVP comparator.
- [x] Add a reverse-mode autodiff direct-solver VJP comparator.
- [x] Add an autodiff iterative-solver comparator where applicable.
- [x] Add a hybrid BLAS or LAPACK custom-rule candidate where measurement
  justifies it.
- [x] Benchmark all implemented linear-algebra candidates.

## Evidence report

- [x] Add cumulative `analytical`, finite-difference diagnostic, `autodiff`,
  and `hybrid` statistics with reproducible `fortplot` figure generators.

## Interpolation and splines

- [x] Add analytical products for active support-node locations.
- [x] Add analytical products for active knot locations where defined.
- [x] Add analytical implicit differentiation of fitted spline coefficients.
- [x] Define derivative status at interpolation-cell crossings.
- [x] Define derivative status at knot crossings.
- [x] Benchmark separate and fused combined-active Lagrange products.
- [x] Benchmark separate and fused combined-active B-spline products.
- [x] Compare analytical basis products with fixed-span autodiff candidates.

## Generated-kernel policy and reproducibility

- [x] Document which symbolic, generated numerical, generated Enzyme, and
  hand-written artifacts are committed, and record the current duplication
  baseline.
- [x] Extend and test `fortsym` provenance support in its own repository when
  the emitter needs capabilities beyond its current public interface.
- [x] Add a `fortsym` revision lock and emit the exact revision in every
  generated numerical-kernel banner.
- [x] Regenerate Dawson, determinant, and inverse kernels with current module,
  purity, line-wrapping, operation-count, and `fo` regeneration metadata.
- [x] Add a temporary-directory regeneration check that byte-compares every
  committed generated numerical kernel.
- [x] Run generated-source drift checking in CTest and CI.
- [x] Record compiler flags, native symbol size, `fortsym` revision, and
  structural operation counts in relevant benchmark records.
- [x] Generate fused value/JVP, fused value/VJP, separate JVP/VJP, and
  contracted products from one symbolic DAG.

## Reliable performance-portable GPU execution

This tranche starts immediately after the preceding single-DAG item. It adopts
the useful Kokkos-like separation of mathematical kernels, execution schedule,
and data layout without implementing containers, execution spaces, a task
graph, a memory manager, or a runtime dispatch framework.

One pure Fortran numerical leaf must serve serial CPU, CPU parallel, OpenMP
target, and OpenACC builds wherever the compiler supports it. `fortsym` owns
explicit mathematical generation; `fortnum` owns small execution wrappers and
application data residency. GPU work initially supplies `analytical`
derivatives. CPU Enzyme remains the `autodiff` backend and supplies CPU
`hybrid` composition.

- [x] Document the supported CPU/GPU differentiation matrix and the minimal
  device-leaf contract: pure procedures, explicit contiguous data, and no
  allocation, I/O, mutable global state, polymorphism, procedure pointers, or
  unsupported descriptors.
- [x] Extend and test `fortsym` device-leaf emission in its own repository,
  limited to optional OpenMP `declare target` and OpenACC `routine seq`
  annotations; add no scheduling or GPU runtime abstraction.
- [x] Pin the tested `fortsym` revision in `fortnum`, regenerate the Dawson
  leaf with device annotations, and prove byte-stable regeneration.
- [x] Add the build-time choice
  `FORTNUM_GPU_BACKEND=NONE|OPENACC|OPENMP`, defaulting to `NONE`; selecting an
  unavailable backend must fail configuration instead of falling back.
- [x] Add one shared benchmark-only batch-wrapper template whose numerical
  loop body is identical for OpenACC and OpenMP target.
- [x] Offload the generated analytical Dawson fused value/JVP through OpenACC
  with a released `nvfortran`; prove NVIDIA-device execution and validate
  against an independent CPU analytical oracle.
- [x] Offload the identical Dawson leaf, layout, and workload through OpenMP
  target with the same released compiler; require
  `omp_is_initial_device()==false`.
- [x] Add persistent-data execution for the Dawson pilot and report
  transfer-inclusive and resident-data wall-clock measurements separately.
- [x] Offload a generated analytical VJP leaf and validate it with an
  independent adjoint identity.
- [x] Benchmark fused versus separate value/JVP and value/VJP GPU kernels
  across launch-bound, throughput, and memory-bound batch sizes.
- [x] Offload a generated multi-input scalar-output kernel and measure
  analytical JVP and VJP scaling over active inputs, derivative directions,
  output cotangents where applicable, and batch sizes.
- [x] Compare admissible array layouts by complete-workload wall clock and
  select one coalesced layout without exposing a general layout framework.
- [x] Record peak device allocation, achieved bandwidth, occupancy, register
  use, and spills with released profiling tools.
- [x] Add reproducible GPU tables and `fortplot` source generators for runtime,
  scaling, memory, bandwidth, and backend crossovers; commit no PNG files.
- [x] Record GPU compiler, device, flags, driver, `fortsym` revision, source
  revision, residency policy, and device-execution proof in benchmark records.
- [x] Select OpenACC or OpenMP target per measured workload before entering
  the launch loop; neither backend name nor derivative mechanism selects a
  winner.
- [ ] Add released Flang/OpenMP-target as a supported compiler combination only
  after it passes the same device-execution, oracle, memory, and wall-clock
  gates.
- [x] Add one analytical GPU implicit-composition pilot using generated local
  residual products and the existing analytical tangent or adjoint boundary;
  classify the complete candidate as `analytical`.
- [x] Expand GPU coverage to special functions with generated Dawson
  value/JVP/VJP leaves, real-device validation, resident and
  transfer-inclusive wall clock, and peak-memory evidence.
- [x] Expand GPU coverage to interpolation with generated fixed-cell Lagrange
  value/JVP/VJP leaves, independent cubic and adjoint oracles, and resident and
  transfer-inclusive wall-clock and peak-memory evidence.
- [x] Expand GPU coverage to fixed quadrature with arbitrary-order analytical
  JVP/VJP batch loops, exact-polynomial and adjoint oracles, and resident and
  transfer-inclusive wall-clock and peak-memory evidence.
- [x] Expand GPU coverage to fixed-size linear algebra with generated 3x3
  determinant and inverse JVP/VJP leaves, finite-difference, matrix-product,
  and adjoint oracles, plus resident and transfer-inclusive wall-clock and
  peak-memory evidence.
- [x] Expand GPU coverage to FFT rules with one shared analytical length-8
  radix-2 JVP/VJP leaf, direct-DFT and complex-adjoint oracles, and resident
  and transfer-inclusive wall-clock and peak-memory evidence.
- [x] Expand GPU coverage to residual kernels through the generated scalar-root
  residual and analytical implicit JVP pilot.
- [x] Expand GPU coverage to fixed-trace ODE products with shared two-state
  forward-tangent and reverse-adjoint transition-map recurrences, closed-form
  and adjoint oracles, trace-length scaling, cache counters, and resident and
  transfer-inclusive wall-clock and peak-memory evidence.
- [x] Run an application-level GPU benchmark with persistent data and select
  candidates from end-to-end wall clock and peak memory rather than isolated
  kernel latency.

### GPU source ownership and exclusions

The symbolic input, selected generated leaf, provenance, regeneration test, and
benchmark evidence are committed. Temporary algebraic variants, generated
Enzyme wrappers, and generated batch wrappers remain reproducible build-tree
artifacts. OpenACC and OpenMP target may have distinct compiler configuration
and scheduling directives, but never distinct copies of a mathematical body.

Do not implement direct Enzyme differentiation of GPU kernels or offload
regions, Enzyme integration into `nvfortran`/OpenACC, custom MLIR, a
CUDA-specific derivative branch, or a Kokkos replacement. Do not add a
CPU-Enzyme outer derivative that launches an analytical GPU custom rule until
persistent-device semantics are proven and an application benchmark requires
that boundary. Reconsider these exclusions only when a released upstream
toolchain supports the complete path reliably and a concrete consumer
justifies the complexity.

- [x] Generate several algebraic variants, prove them equivalent, and rank
  them by post-CSE operation count before native CPU, OpenACC, and OpenMP-target
  benchmarking.
- [x] Commit only the measured production kernel and the generator inputs
  needed to reproduce temporary losing variants.

## Shared CPU Enzyme infrastructure

This section removes CPU fixture duplication. It does not imply GPU Enzyme
support.

- [x] Record a machine-readable baseline for Enzyme fixture duplication,
  build time, runtime, peak RSS, and native code size.
- [x] Add one internal support module for environment parsing, timing, warmups,
  median/MAD, standardized output, and peak-RSS access.
- [x] Generate wrappers for the proven scalar-kernel shape with one to four
  active scalar inputs, scalar output, forward JVP, reverse VJP, and an
  optional analytical forward custom rule.
- [x] Add shared custom-rule counter support that proves analytical rule
  selection without kernel-specific counter boilerplate.
- [x] Migrate the Bessel fixture as the pilot and require equivalent validation
  plus no complete-workload regression beyond 3% or combined dispersion.
- [x] Migrate Dawson and the scalar Enzyme smoke fixtures.
- [x] Migrate fixed-span B-spline Enzyme fixtures.
- [x] Migrate direct-solver JVP and VJP Enzyme fixtures.
- [x] Migrate the iterative-solver Enzyme fixture.
- [x] Migrate scalar-root JVP and VJP Enzyme fixtures.
- [x] Migrate vector-root JVP and VJP Enzyme fixtures.
- [x] Migrate fixed-quadrature JVP and VJP Enzyme fixtures.
- [x] Migrate adaptive-integration Enzyme fixtures.
- [x] Migrate the ODE forward-sensitivity Enzyme fixture.
- [x] Remove superseded duplicate helpers and reject new copies in repository
  checks.

## Analytical and hybrid generation migration

- [x] Establish a compact documentation map and rewrite the README, contributor
  guide, benchmark guide, and repository templates around current behavior.
- [x] Rewrite the API guide and module design documents as concise contracts
  checked against the production source.
- [x] Rewrite the differentiation, GPU, migration, and evidence documents;
  remove append-only implementation history while preserving authoritative
  benchmark data and decisions.
- [x] Add documentation drift checks for internal links, source ownership,
  terminology, generated revisions, and stale future-tense claims.

- [x] Inventory every derivative kernel as `fortsym`-generated, hand-written
  algorithmic, hand-written stable recurrence, implicit solve, frozen trace,
  or generated hybrid boundary.
- [x] Replace manually transcribed explicit chain-rule expressions with
  `fortsym` value/JVP/VJP generation where the generated candidate is admissible.
- [x] Generate contracted JVP expressions directly without materializing a
  Jacobian.
- [x] Generate contracted VJP expressions directly without materializing a
  Jacobian.
- [x] Generate fused and separate value/product variants from one symbolic
  definition.
- [x] Generate special-function outer kernels while retaining stable primal
  algorithms and recurrences at operator boundaries.
- [x] Generate explicit scalar local residual products while retaining
  hand-written implicit tangent and adjoint solves.
- [x] Generate explicit vector local residual products while retaining
  hand-written implicit tangent and adjoint solves.
- [x] Add stability-preserving symbolic transformations with equivalence and
  numerical boundary tests.
- [x] Add target-aware factoring only when native complete-workload benchmarks
  justify it.
- [x] Re-run candidate selection and cumulative reports after the generation
  migration.

## Special functions

- [x] Run an analytical, autodiff, and hybrid Bessel tournament.
- [x] Run an analytical, autodiff, and hybrid gamma-function tournament.
- [ ] Run an analytical, autodiff, and hybrid error-function tournament.
- [ ] Run an analytical, autodiff, and hybrid hypergeometric tournament.
- [ ] Generate region-specific analytical candidates with `fortsym`.
- [ ] Benchmark fused and separate value/JVP special-function candidates.
- [ ] Benchmark representative small-argument regions.
- [ ] Benchmark representative transition regions.
- [ ] Benchmark representative asymptotic regions.

## FFT and transforms

- [ ] Benchmark the analytical FFT JVP against autodiff of the implementation.
- [ ] Benchmark the analytical FFT VJP against autodiff of the implementation.
- [ ] Verify complex-adjoint and normalization conventions independently.
- [ ] Add a hybrid custom rule for an external FFT library.
- [ ] Benchmark scalar FFT workloads.
- [ ] Benchmark batched FFT workloads.
- [ ] Benchmark an FFT-based spectral objective.

## Module tournaments

- [ ] Complete the special-functions tournament and commit its evidence.
- [ ] Complete the interpolation tournament and commit its evidence.
- [ ] Complete the FFT tournament and commit its evidence.
- [ ] Complete the quadrature tournament and commit its evidence.
- [ ] Complete the roots tournament and commit its evidence.
- [ ] Complete the linear-algebra tournament and commit its evidence.
- [ ] Complete the ODE tournament and commit its evidence.

Every tournament record must include:

- [ ] Independent validation error.
- [ ] Median runtime and dispersion.
- [ ] Candidate-specific peak memory.
- [ ] Generated or native code size where relevant.
- [ ] Hardware and compiler identity.
- [ ] Reusable primal state.
- [ ] The selected candidate and deterministic selection rationale.
- [ ] Scaling over representative active-input and output counts.
- [ ] Forward-mode versus reverse-mode crossover evidence.
- [ ] CPU work and cache counters where supported.

## Application-level benchmarks

- [ ] Benchmark a root-constrained scalar objective.
- [ ] Benchmark a parameterized integral.
- [ ] Benchmark spline fitting followed by a downstream gradient.
- [ ] Benchmark an FFT-based spectral objective.
- [ ] Benchmark nonlinear least squares.
- [ ] Benchmark stiff ODE sensitivity.
- [ ] Benchmark nonstiff ODE sensitivity.
- [ ] Benchmark an event-driven trajectory.
- [ ] Benchmark a small PDE or residual application.
- [ ] Benchmark a representative `itpplasma` workload.
- [ ] Measure end-to-end forward/reverse scaling over active-input counts.
- [ ] Measure end-to-end forward/reverse scaling over output counts.
- [ ] Record wall-clock, peak-memory, and cache-performance crossover curves.

## Second order

- [ ] Identify a concrete consumer that requires an HVP.
- [ ] Add an autodiff forward-over-reverse HVP candidate.
- [ ] Add an autodiff reverse-over-forward HVP candidate.
- [ ] Add an analytical or hybrid implicit HVP where justified.
- [ ] Validate HVPs with an independent second-directional or high-precision
  oracle.
- [ ] Measure HVP runtime and peak memory.
- [ ] Add a full Hessian only when a measured consumer requires one.

## Repository quality

- [ ] Resolve or explicitly baseline the `fo lint` array-temporary warnings.
- [x] Repair the fpm and `fo test` oracle-argument mismatch.
- [ ] Add CI coverage for the real Flang and Enzyme hybrid pipeline.
- [ ] Add scheduled benchmark runs on the reference Ryzen 9 5950X.
- [ ] Re-run candidate selection after material compiler changes.
- [ ] Re-run candidate selection after material Enzyme changes.
- [ ] Re-run candidate selection after material `fortsym` changes.
- [ ] Re-run candidate selection after material primal-code changes.
- [ ] Re-run candidate selection on new target hardware.
