# Derivative contract

Status: normative.

## Terminology

`fortnum` uses three public mechanism names:

- `autodiff`: compiler or source-transformation differentiation of program
  operations
- `analytical`: explicit expressions, stable recurrences, linear operators,
  implicit tangent or adjoint solves, and frozen-trace algorithms
- `hybrid`: autodiff composition with at least one analytical rule at a
  mathematical operator boundary

These names describe candidate implementations. They are not exclusive classes
assigned to procedures.

Finite differences and complex step are validation diagnostics by default.
They become production candidates only when their accuracy, wall clock, and
memory satisfy the same gates.

## Products

For \(f:\mathbb R^n\rightarrow\mathbb R^m\), callers request:

| Product | Definition | Typical scaling pressure |
| --- | --- | --- |
| JVP | \(Jv\) | favorable for few input directions |
| VJP | \(J^Tu\) | favorable for few output cotangents |
| gradient | \(J^T\) for a scalar output | reverse or analytical adjoint |
| HVP | \(\nabla^2 f\,v\) | mixed forward/reverse or analytical second order |
| Jacobian | all entries of \(J\) | storage and recovery cost |

Public derivative names append `_jvp`, `_vjp`, `_grad`, or `_hvp` to the
primal operator. A product may return the primal value in the same call when
fusion shares work.

Full Jacobians and Hessians are added only for consumers that need the matrices.
Contracted products are the default.

## Active arguments

Every product states which arguments are active. Values that select control
flow or representation are inactive unless a separate mathematical contract
defines their derivative:

- dimensions and integer selectors
- tolerances and iteration limits
- adaptive accept/reject decisions
- cell and knot-span indices
- random seeds and counters
- status and diagnostic outputs

Activity is a property of one product. A knot location can be active in one
B-spline product and inactive in a coefficient-only product.

## Candidate boundaries

Custom rules follow mathematical operators:

| Operator | Leading analytical candidates |
| --- | --- |
| explicit local expression | generated contracted JVP/VJP |
| stable special function | identity or recurrence around the stable primal |
| FFT or fixed quadrature | operator and adjoint operator |
| root, fixed point, equilibrium | implicit tangent or adjoint solve |
| linear solve | factor-reusing tangent or transpose solve |
| adaptive integration or ODE | frozen trace, continuous sensitivity, or discrete adjoint |
| event time | implicit derivative of the event residual |
| interpolation | fixed-cell or fixed-span basis products |

Autodiff can differentiate smooth code between these boundaries. A `hybrid`
candidate may use autodiff for local residual products and an analytical
implicit solve for the outer operator.

## Implicit products

For a state defined by

\[
R(y,p)=0,
\]

the analytical tangent solves

\[
R_y\,\dot y=-R_p\,\dot p.
\]

For a scalar objective \(L(y,p)\), the adjoint solves

\[
R_y^T\lambda=L_y^T,
\qquad
\frac{dL}{dp}=L_p-\lambda^T R_p.
\]

The implementation solves these systems. It never forms \(R_y^{-1}\).
Converged Jacobians, factorizations, and preconditioners are reusable primal
state. Near-singular systems report derivative reliability explicitly.

Differentiating solver iterations remains an admissible comparator. It is not
the unmeasured default.

## Adaptive products

An adaptive primal records its accepted schedule. A frozen-trace product
differentiates that fixed numerical map. It does not differentiate accept/reject
logic or a changing step count.

A continuous sensitivity differentiates the underlying mathematical problem.
A discrete adjoint differentiates the numerical step composition. Their
answers can differ at finite resolution. Each interface states which contract
it implements.

## Provenance and quality

`fortnum_ad_status_t` combines primal status with:

- mechanism provenance
- derivative quality

Quality is independent of mechanism. Exact, approximate, and nonsmooth results
can come from different implementations. Status merging preserves the worst
quality and every relevant provenance bit.

Legacy backend constants remain ABI aliases. New documentation and code use
the public terminology in this document.

## Selection

Each candidate is independently validated before timing. Production selection
uses the requested product and a workload class that records relevant
dimensions, including:

- active inputs and outputs
- JVP directions or VJP cotangents
- batch size
- reusable factors, traces, or plans
- CPU/GPU target and data residency

Complete-workload median wall clock is the primary metric. Candidate-specific
peak memory is mandatory. Dispersion, code size, CPU work, cache counters, data
transfer, bandwidth, occupancy, registers, and spills are recorded when they
affect the decision.

Selection occurs before hot loops. Unknown workload classes return unavailable
or use an explicitly documented fallback. They are not silently guessed.

## `fortsym` ownership

`fortsym` owns symbolic DAG construction, differentiation, simplification,
common-subexpression elimination, contracted products, fused products, and
Fortran emission. `tools/codegen/fortsym.lock` identifies the tested revision.

Commit:

- symbolic specifications and generators
- the selected production kernel
- provenance and regeneration commands
- independent tests and benchmark evidence

Keep routine Enzyme wrappers and losing algebraic candidates in the build tree.
Do not commit generated figures.

Hand-written code retains:

- stable recurrences and region selection
- adaptive and solver orchestration
- factorization and preconditioner reuse
- implicit tangent and adjoint solves
- CPU/GPU execution schedules
- independent validation oracles

The machine-checked ownership list is
[derivative_kernel_inventory.csv](derivative_kernel_inventory.csv).

## GPU scope

CPU supports selected `autodiff`, `analytical`, and `hybrid` candidates. GPU
support currently consists of validated generated `analytical` leaves.
Direct Enzyme differentiation of offload regions is unsupported.

See [gpu.md](gpu.md) for the device contract.

## Adding a product

1. define the mathematical operator and requested product
2. declare active and inactive arguments
3. enumerate admissible mechanisms
4. generate mechanical local algebra where appropriate
5. add an independent oracle
6. test branches, regimes, singularities, and adjoint identities
7. benchmark realistic workload classes
8. record selection evidence
9. update the API, inventory, and ROADMAP
