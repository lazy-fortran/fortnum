# ADR: Derivative contract for every public module

Status: accepted (issue #37, M6.1), amended by
`docs/performance_optimal_differentiation.md`. Normative for all later module
ADRs and implementation issues.

fortnum is primal-first. A routine computes a value, and that value is correct
on its own, before any derivative exists. But every public procedure is written
so a derivative can be added later without changing the primal signature or the
primal behavior. This document fixes the contract that makes "later" cheap: what
a derivative is called, which arguments carry it, and which candidate
mechanisms may produce it.

The contract is the shared vocabulary. A module ADR does not redefine these
terms; it cites the section numbers here and fills in the module-specific
choices. `CONTRIBUTING.md` summarizes the candidate terms for contributors
and names this file as the source of truth. Where the two disagree, this file
wins.

## 1. Candidate mechanisms

A public procedure may have several implementations of the same derivative
product. The public terminology is:

- `autodiff`: a compiler or source-transformation backend differentiates a
  smooth implementation.
- `analytical`: an explicit expression, recurrence, implicit equation,
  tangent or adjoint solve, fixed linear operator, or frozen-trace rule supplies
  the product.
- `hybrid`: autodiff composes a larger smooth calculation while one or more
  mathematical operators use analytical custom rules at their interfaces.

These terms describe implementations, not exclusive procedure classes. For
example, `root_jvp` can have an analytical implicit implementation, an autodiff
implementation that differentiates solver iterations, and a hybrid
implementation that obtains residual products by autodiff and applies the
analytical implicit rule at the solver boundary.

Candidate metadata may refine the public terms with mechanisms such as
`autodiff_forward`, `autodiff_reverse`, `analytical_expression`,
`analytical_recurrence`, `analytical_implicit_tangent`,
`analytical_implicit_adjoint`, `analytical_frozen_trace_forward`,
`analytical_frozen_trace_reverse`, and `hybrid_custom_rule`.

`finite_difference_reference` is a validation or fallback candidate.
`primal_only` means no derivative is mathematically meaningful for the declared
active arguments. A `primal_only` declaration carries a one-line
justification.

For each derivative product, documentation records:

1. every mathematically admissible candidate
2. the independent validation oracle
3. the workload classes that affect selection
4. benchmark evidence for the selected implementation

No candidate wins from its mechanism name. Production selection minimizes
measured application runtime and peak memory subject to the derivative contract
and validation tolerance.

## 2. Naming convention

For a primal routine `foo`, the public derivative entry points are:

- `foo_jvp(...)`: forward-mode product (Jacobian times vector).
- `foo_vjp(...)`: reverse-mode product (vector times Jacobian).
- `foo_grad(...)`: gradient, where the output is scalar and a gradient is what
  the caller wants.
- `foo_hvp(...)`: Hessian-vector product, where second order is meaningful.
- `foo_ad(...)`: optional combined convenience wrapper over the above.

Not every primal needs every name. A scalar-output routine offers `foo_grad`; a
vector-output routine offers `foo_jvp` and `foo_vjp`. Add the names the
derivative contract and output shape justify, no more.

Autodiff-generated entry points stay internal. Public callers use the `foo_*`
names above and never call a raw `__enzyme_*` symbol. Generated symbols,
analytical rules, dispatch, and registry data are implementation details.

## 3. Active argument rules

An argument is active if the derivative flows through it, inactive otherwise.
The classification is part of the contract and is declared per procedure.

- Active real and complex arrays use contiguous explicit-shape or assumed-size
  wrappers on the autodiff path first. The simplest memory layout is the one
  the compiler path is tested against. Add richer layouts only with dedicated
  tests.
- Descriptors, allocatables, polymorphism, optional active arguments, and
  derived-type components are allowed on an autodiff path only when a dedicated
  compiler test covers that shape. No test means no autodiff support. The primal
  may still use them freely.
- Integers, sizes, orders, keys, branch modes, RNG seeds, status flags, and
  workspace capacities are inactive. They select behavior or report it; they do
  not carry a derivative.

The `fortnum_status_t` object is inactive. Status reporting is a side channel,
not a differentiable output.

## 4. Per-module candidate table

The table lists leading candidates. It does not preselect a winner.

| Module | Leading candidates | Required operator boundary |
|---|---|---|
| `special` | analytical expressions and recurrences, autodiff of the numerical approximation, hybrid custom rules | Preserve argument regions, branch cuts, scaling, and recurrence direction. |
| `fft` | analytical transform and adjoint transform, autodiff of the implementation, hybrid external-library rule | Treat the transform as a linear operator. |
| `quadrature` | analytical weighted sum, autodiff of the fixed rule | Nodes and weights are inactive unless their parameter dependence is explicitly requested. |
| `levin` | analytical rational-transform products, frozen-selection variants, finite-difference reference | The selected order is inactive within a candidate and nonsmooth at a selection change. |
| `integrate` | analytical differentiation under the integral, analytical frozen trace, autodiff over the frozen trace, hybrid integrand autodiff | Do not make refinement logic the only derivative definition. |
| `ode` | analytical forward sensitivity, analytical discrete or continuous adjoint, autodiff over a frozen trace, hybrid autodiff RHS products | State whether the derivative contract is discrete or continuous. |
| `roots` | analytical implicit tangent and adjoint, hybrid autodiff residual products, autodiff of iterations as a measured comparator | Reuse the converged Jacobian, factorization, or preconditioner where possible. |
| `linalg` | analytical implicit solve rule, factorization-specific rule, autodiff solver comparator, hybrid external-library rule | Use solves and transpose solves, never an explicit inverse by default. |
| `interp` | analytical basis products, autodiff inside a fixed cell, analytical implicit coefficient solve, hybrid composition | Cell and knot-span selection is inactive. |
| `rng` | `primal_only` for seed-to-draw kernels | Estimator derivatives belong above the RNG layer. |

## 5. How module docs reference this contract

A module ADR or an implementation issue is normative downstream of this file. It
states, for each public procedure:

1. the admissible `autodiff`, `analytical`, and `hybrid` candidates from
   section 1
2. the active and inactive arguments per section 3
3. the derivative entry-point names per section 2
4. the independent validation oracle
5. benchmark workload classes and evidence for any selected winner

It cites this document by section number rather than restating the definitions.
A reviewer checks the module doc against these sections. Missing admissible
candidates may be recorded as planned work, but a selected production winner
requires evidence.

## 6. Adding derivatives without API churn

The primal ships first (M1 through M5); derivative products arrive later (issue
#40). The contract is built so that arrival changes no existing signature.

- Derivative routines are new public names under the section 2 convention. They
  are added beside the primal, never by adding arguments to it. A caller who
  only wants the value keeps calling `foo` unchanged.
- Primal routines stay `pure` where the algorithm allows and hold no
  module-level mutable state. Purity supports autodiff and lets derivative
  wrappers call the primal as many times as a product needs.
- The reserved derivative slots already exist. `fortnum_status_t` carries the
  status side channel that derivative routines reuse without a new type, and the
  oracle CSV format already reserves a derivative column behind the
  `has_derivative` header flag (`test/oracle/data/*.csv`). When `foo_jvp` lands,
  its oracle table is the same file with `has_derivative: 1` and the column
  filled; no format change, and tables written before then set the flag to `0`
  so consumers ignore the empty column.

Adding a derivative is therefore additive: new names, the same primal, the same
status type, the same CSV layout.

## 7. Candidate generation and `fortsym`

`fortnum` uses `fortsym` for symbolic algebra and derivative code generation.
The development checkout currently resolves it at `../lazy-fortran/fortsym`.
The supported integration surface is deliberately narrow and demonstrated by
`tools/codegen/app/gen_dawson_outer.f90`: expression construction, differentiation,
engine simplification, fused kernel emission, post-CSE operation counting, and
an exact regeneration command in the generated banner.

`fortsym` remains under development. New work may extend it, with independent
tests in that repository, but must not guess an unimplemented API or duplicate a
symbolic engine inside `fortnum`. Generated kernels are build-time artifacts
committed under `src/generated/`; production `fortnum` does not dynamically
depend on `fortsym`.

### 7.1 Source ownership

The symbolic specification, generator, selected numerical kernel, provenance,
and benchmark evidence are committed. A selected generated kernel must remain
buildable without `fortsym`; normal consumers never run a CAS.

Temporary algebraic variants and routine Enzyme wrappers are generated in the
build tree. They are reproducible from committed inputs but are not committed.
Generated PNG reports are likewise excluded. Shared Enzyme timing, statistics,
RSS, and environment support is hand-written once and committed.

The following code remains hand-written when symbolic expansion would erase
useful numerical structure:

- stable recurrences and region selection
- adaptive and frozen traces
- solver, factorization, and preconditioner orchestration
- implicit tangent and adjoint solves
- kernel-specific independent validation oracles

`fortsym` generates explicit local algebra, contracted products, fused
value/product kernels, and outer chain rules around those operator boundaries.
Candidate selection still uses native complete-workload measurements.

### 7.2 Current duplication baseline

Before the shared Enzyme migration, `cmake/enzyme/hybrid/` contains 5,190
Fortran and C lines across 20 fixtures. The repeated declarations include 12
forward Enzyme interfaces, 7 reverse interfaces, 13 peak-RSS interfaces, and 9
copies each of sorting and result-report code. This is the baseline for source
reduction; numerical kernels and independent oracles are not counted as
removable boilerplate.

This documentation-only policy change affects no compiled source. The current
reference performance baseline remains the committed 32-tournament data set.
For example, the 16-product Bessel series JVP takes 917.276 ns for the selected
analytical candidate, and the recurrence JVP takes 3,788.275 ns for selected
autodiff.

## 8. Selection and dispatch

Each product is selected independently. A fast Jacobian implementation does not
automatically select the fastest JVP, VJP, gradient, or HVP.

Selection may be static at build time, generated from a benchmark registry, or
explicitly configured. Dynamic selection must be resolved before entering a
hot loop. Registry evidence records the operator, product, workload class,
hardware, validation error, runtime, and peak memory.

## 9. Hybrid interface rule

A hybrid boundary corresponds to a mathematical operator. Autodiff may run
inside the local residual, outside the operator, or both. The operator call uses
an analytical JVP, VJP, tangent solve, adjoint solve, recurrence, transform, or
frozen-trace rule.

The first hybrid implementation should use simple interoperable scalar or
explicit-array arguments. Optional active arguments, polymorphism,
allocatables, and Fortran descriptors require dedicated compiler tests before
support.
