# Differentiation implementation plan

Status: implementation in progress. Normative downstream of
`docs/design/ad.md` and motivated by
`docs/performance_optimal_differentiation.md`.

## 1. Goal

`fortnum` will be primal-first, derivative-plural, and benchmark-selected.
Every important derivative product can have `autodiff`, `analytical`, and
`hybrid` candidates. Independent validation establishes correctness. Measured
application runtime and peak memory select production implementations.

The stable public surface remains:

```text
foo(...)
foo_jvp(...)
foo_vjp(...)
foo_grad(...)
foo_hvp(...)
```

Candidate generation, registration, selection, and generated symbols remain
internal.

## 2. Current baseline

The library already contains:

- backend-independent value, JVP, VJP, gradient, and HVP interfaces
- a flat active-vector layout for optimizer-facing kernels
- approximately 45 first-order derivative procedures
- analytical special-function identities and recurrences
- analytical FFT and fixed-quadrature linear products
- analytical frozen-trace products for adaptive integration and ODE integration
- analytical implicit products for scalar and vector roots
- generic analytical root tangent and adjoint callback boundaries
- analytical implicit fixed-point tangent and adjoint products
- analytical interpolation and B-spline products
- finite-difference, complex-step, smoothness, and adjoint-identity test helpers
- an optional Flang and Enzyme test pipeline with three isolated smoke kernels

The first implementation slice now also contains:

- a real Flang 22 and Enzyme pipeline, including C analytical custom rules
- a generated, public fused Dawson outer-expression value/JVP kernel
- a real `hybrid` Dawson boundary whose test proves the custom rule was called
- deterministic candidate metadata and a pre-loop selector
- combined-active Lagrange and B-spline product rules
- analytical implicit linear-solve JVP and VJP products
- reusable linear-solve JVP and transposed-VJP factorizations
- generic scalar-root and vector-root tangent and adjoint boundaries
- analytical fixed-point JVP and VJP products
- hybrid scalar-root JVP and VJP candidates using Enzyme forward and reverse
  residual products
- hybrid vector-root JVP and VJP candidates using Enzyme forward and reverse
  residual products
- optional reuse of the converged analytical vector-root Jacobian and its LU
  or transposed-LU factorization across JVP directions or VJP cotangents
- reciprocal-condition reporting and threshold rejection on the vector-root
  implicit JVP boundary
- committed runtime, dispersion, memory, validation, hardware, and toolchain
  evidence for the first tournament

The main remaining pieces are hybrid residual products, root-factorization
reuse, broad module tournaments, application-level selection, and justified
second-order products.

## 3. Terminology

- `autodiff`: compiler or source-transformation differentiation of smooth code.
- `analytical`: explicit expressions, recurrences, implicit solves, tangent or
  adjoint models, linear operators, and frozen-trace rules.
- `hybrid`: autodiff composition with analytical rules at mathematical operator
  interfaces.

Use these terms in new code, documentation, issues, and pull requests. Legacy
terms such as `transparent`, `analytic_rule`, `implicit_rule`, and `trace_rule`
may describe existing implementations during migration. They no longer imply
exclusive ownership of a procedure.

## 4. Symbolic algebra dependency

`fortnum` uses `fortsym` for build-time symbolic algebra and code generation.
The development dependency is `../lazy-fortran/fortsym`. The first generator
fixes the currently exercised contract: expression DAGs, differentiation,
SymEngine-backed simplification, fused kernel emission, operation counting, and
exact regeneration banners.

`fortsym` is still evolving. Extend and test it there when a planned candidate
needs more capability. Do not invent an API in this repository, and do not add a
runtime dependency on it.

## 4.1 Implementation status

| Plan area | Status | Evidence |
|---|---|---|
| terminology and candidate contract | complete | `docs/design/ad.md`, `AGENTS.md` |
| status/provenance composition | complete | `ad_status_merge` behavioral tests |
| generated analytical JVP | complete for first slice | `dawson_outer_jvp`, generated source and finite-difference oracle |
| real `hybrid` boundary | complete for first slice | Enzyme Dawson custom-rule call assertion |
| candidate registry | complete for static selection | deterministic validation, timing, memory, code-size and ID ordering tests |
| symbolic generation | complete for first slice | `gen_dawson_outer`, `fortsym` tests, regeneration banner |
| interpolation interface rules | complete | simultaneous-activity directional and adjoint tests |
| implicit linear solve | JVP/VJP and factorization reuse complete | finite-difference, adjoint, and reuse benchmarks |
| roots and fixed points | scalar/vector tournaments, analytical boundaries, hybrid JVP/VJP, Jacobian/factor reuse, and implicit JVP/VJP reliability complete | complete-solve finite-difference, exact-condition, Jacobian, and scalar-objective oracles |
| integration | analytical differentiation under fixed bounds complete; moving bounds, hybrid products, and tournaments pending | closed-form derivative and complete-integral finite differences |
| ODE hybridization | pending | existing analytical products remain candidates |
| module/application tournaments | pending except Dawson | first committed table in `differentiation_benchmarks.md` |
| second order | pending | implement only for demonstrated consumers |

## 5. Work plan

Difficulty is relative to this repository:

- S: one focused change
- M: several coordinated modules or tests
- L: substantial build, compiler, or numerical work
- XL: compiler-sensitive work with research risk

| Order | Terminology | Location | Deliverable | Difficulty | Completion gate |
|---:|---|---|---|---|---|
| 1 | all | `docs/`, contributor templates, repository instructions | Candidate-set contract, standard terminology, `fortsym` ownership rule | S | No instruction still requires exactly one derivative policy. |
| 2 | `analytical` | `src/ad/fortnum_ad_interfaces.f90` | Hybrid provenance and deterministic quality/status merging | S | Unit tests cover mixed backend and worst-quality propagation. |
| 3 | `autodiff` | `cmake/`, `src/ad/`, `src/CMakeLists.txt` | One generated JVP linked into the production library | L | Public product uses a real generated symbol and passes an independent directional oracle. |
| 4 | `hybrid` | `src/ad/`, special-function wrappers | Autodiff outer expression with an analytical Dawson or Bessel rule at the call boundary | M to L | Pure autodiff, analytical, and hybrid candidates agree with finite difference, and the test proves the analytical rule was selected. |
| 5 | `analytical` and `hybrid` | `src/roots/fortnum_roots.f90`, `src/roots/fortnum_multiroot.f90` | Generic implicit tangent and adjoint boundary with autodiff residual products as candidates | M to L | Residual linearization and adjoint identities pass, including ill-conditioning status. |
| 6 | `analytical` and `hybrid` | `src/quadrature/fortnum_integrate.f90` | Autodiff integrand products composed with analytical differentiation under the integral or frozen trace | M to L | Moving-boundary and fixed-boundary contracts are distinguished and independently validated. |
| 7 | `analytical` | `src/linalg/fortnum_linalg.f90` | JVP/VJP rules for determinant, inverse, and solve with factorization-reuse hooks | M | Products pass finite-difference and adjoint tests without forming explicit inverses in solve rules. |
| 8 | `analytical` | interpolation and B-spline modules | Combined-active JVP/VJP for evaluation point and coefficients or nodal values | S to M | Product-rule results pass directional and adjoint tests across fixed cells and spans. |
| 9 | `analytical` and `hybrid` | `src/ode/fortnum_ode.f90`, ODE method modules | Autodiff RHS products inside analytical forward sensitivity and discrete adjoint traces, including parameter VJP accumulation | L to XL | Discrete derivative contract passes frozen-trace finite difference and adjoint tests. |
| 10 | all | new internal registry and `benchmark/` | Candidate metadata, workload classes, benchmark evidence, and pre-loop selector | L | Registry reproduces selected winners and records runtime, peak memory, and validation error. |
| 11 | all | module tournaments | Special, interpolation, FFT, quadrature, roots, linear algebra, and ODE candidate tournaments | L | Each selected winner has committed validation and representative benchmark evidence. |
| 12 | all | downstream mini-applications | Root-constrained objective, parameterized integral, spline fit, FFT objective, and ODE sensitivity benchmarks | L | Selection minimizes complete workload cost, not only isolated derivative latency. |
| 13 | all | second-order modules | HVP candidates where consumers justify them | L to XL | Independent second-directional or high-precision oracle and memory evidence exist. |
| 14 | `analytical` and `hybrid` | `fortsym` build-time integration | Symbolic DAG, simplification variants, contracted-product generation, operation counts, and code generation | M to L | Every committed generated kernel records an exact regeneration command and has an independent behavioral oracle. |

## 6. First vertical slice

The first implementation milestone is intentionally small:

1. choose a scalar special function with an existing analytical JVP
2. expose a compiler-stable primal and analytical rule boundary
3. generate an autodiff derivative of an outer nonlinear expression
4. make the generated derivative use the analytical rule at the inner call
5. compare the `autodiff`, `analytical`, and `hybrid` candidates against an
   independent finite-difference or complex-step oracle
6. benchmark runtime and peak memory

Dawson or a real modified Bessel function is the preferred target. This proves
the hybrid mechanism before adding solver traces, tapes, or polymorphic
contexts.

## 7. Candidate tournament

Every performance-relevant product follows the same sequence:

1. state the mathematical operator and active arguments
2. state the requested product
3. enumerate admissible `autodiff`, `analytical`, and `hybrid` candidates
4. implement an independent behavioral oracle
5. validate directional, adjoint, implicit-residual, branch, and regime behavior
6. benchmark realistic workload classes
7. record runtime, peak memory, accuracy, and reusable primal state
8. select the winner outside hot loops
9. repeat after material changes to primal code, compiler, hardware, or workload

Finite differences are reference candidates by default. They become production
fallbacks only when their validation and performance justify that role.

## 8. Interface rules

- Custom-rule boundaries follow mathematical operators.
- Autodiff may run outside a boundary, inside its local residual or physics
  kernel, or both.
- Implicit candidates differentiate defining equations.
- Adaptive candidates state whether they differentiate a frozen discrete trace
  or a continuous mathematical problem.
- Contracted JVPs and VJPs are preferred to full matrices unless reuse justifies
  matrix construction.
- Candidate selection is per product. A Jacobian winner does not select the JVP
  or VJP winner.
- Dispatch resolves before hot loops.
- Simple interoperable scalar and explicit-array boundaries come first.
  Optional active arguments, descriptors, allocatables, and polymorphism require
  dedicated compiler tests.

## 9. Validation rules

Repository-state checks are not tests. Each derivative candidate needs an
independent behavioral oracle appropriate to the product:

- central finite difference with a step-size convergence check
- complex step where the primal is analytic and complex-safe
- high-precision reference
- adjoint identity $u^T(Jv)=v^T(J^Tu)$
- implicit residual identity $R_y\,dy+R_p\,dp=0$
- an independently derived analytical result

A comparison between two implementations generated from the same expression is
corroboration, not an independent oracle.

## 10. Performance evidence

Record at least:

- value time
- derivative time
- fused value-and-derivative time
- peak memory
- validation error
- workload dimensions and direction count
- reusable traces, factorizations, or preconditioners
- compiler and hardware identity

Application-level evidence has priority over isolated kernel latency when the
application changes optimizer iterations, cache behavior, batching, or memory
pressure.

The first measured table and its machine-readable record are in
`docs/design/differentiation_benchmarks.md` and `benchmark/reference/`.
