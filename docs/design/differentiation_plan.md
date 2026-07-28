# Differentiation implementation plan

Status: active. [ROADMAP.md](../../ROADMAP.md) is the authoritative checklist.

## Objective

`fortnum` is primal-first, derivative-plural, and benchmark-selected.

For each important mathematical operator and derivative product:

1. enumerate admissible `autodiff`, `analytical`, and `hybrid` candidates
2. generate mechanical local algebra and wrappers
3. validate each candidate independently
4. measure realistic complete workloads
5. select by wall clock and peak memory

The stable public surface remains:

```text
foo
foo_jvp
foo_vjp
foo_grad
foo_hvp
```

Candidate symbols, generated wrappers, and selection metadata remain internal.

## Implemented foundation

| Area | Current capability |
| --- | --- |
| contract | product-specific candidate sets and standard terminology |
| status | merged provenance, quality, and reliability |
| selection | deterministic workload metadata and pre-loop dispatch |
| generation | one-DAG values, contracted JVP/VJP, fused products, variants, provenance, device annotations |
| CPU autodiff | tested Flang/Enzyme scalar, fixed-vector, and fixed-array wrapper profiles |
| hybrid boundaries | special functions, roots, integration, direct solves, and ODE RHS products |
| implicit products | scalar/vector roots, fixed points, linear solves, spline fits, ODE stages, and events |
| traces | adaptive integration and ODE tangent/adjoint replay |
| GPU | generated analytical leaves for selected fixed kernels and one application |
| evidence | machine-readable CPU/GPU records and reproducible `fortplot` reports |

The ownership inventory at
[derivative_kernel_inventory.csv](derivative_kernel_inventory.csv) distinguishes
generated explicit algebra, hand-written algorithms, stable recurrences,
implicit solves, frozen traces, and generated Enzyme boundaries.

## Generation policy

Use one mathematical source for a local explicit operator. `fortsym` may emit:

- value-only kernels
- contracted JVP and VJP kernels
- fused value/JVP and value/VJP kernels
- raw, simplified, and factored variants
- CPU and annotation-only device leaves
- temporary Enzyme wrappers

Generate a full Jacobian only when a caller reuses it enough to justify
materialization.

Keep the following hand-written:

- special-function regime selection and stable recurrences
- adaptive trace construction
- root, equilibrium, and factorization orchestration
- implicit tangent and adjoint solves
- checkpoint and recomputation schedules
- GPU scheduling and data residency
- independent validation oracles

Generated code must pass symbolic equivalence where available, numerical
boundary tests, byte-stable regeneration, native compilation, and complete
workload benchmarks.

## Candidate tournament

Every tournament records:

| Dimension | Required evidence |
| --- | --- |
| operator | mathematical definition and primal contract |
| product | JVP, VJP, gradient, HVP, or matrix |
| activity | active and inactive arguments |
| validation | independent error, adjoint identity, or residual equation |
| workload | inputs, outputs, directions, cotangents, batch, reusable state |
| runtime | median complete-workload wall clock and dispersion |
| memory | candidate-specific peak memory |
| scaling | forward/reverse crossover dimensions |
| hardware | compiler, flags, CPU/GPU, affinity, and revisions |
| counters | cache, work, code size, transfer, or device counters when useful |
| selection | winner and deterministic rationale |

Value-plus-product and derivative-only workloads are different contracts. A
comparison must return the same values before it can support “faster”.

## CPU and GPU roles

CPU remains the complete feature path. Enzyme supplies CPU `autodiff`
candidates and CPU `hybrid` composition.

GPU uses explicit generated `analytical` leaves. OpenACC and OpenMP target
share the same mathematical body. Direct GPU autodiff is outside the supported
path.

## Remaining sequence

The current ROADMAP sequence is:

1. migrate admissible hand-written explicit chain rules to generated products
2. finish special-function and FFT tournaments
3. consolidate module-level tournament evidence
4. run complete downstream applications
5. add second-order products only for demonstrated consumers
6. keep CI, benchmark selection, and documentation synchronized with compiler,
   Enzyme, `fortsym`, primal, and hardware changes

Each ROADMAP checkbox is implemented, tested, measured, committed, pushed, and
installed before the next begins.
