# ADR: Derivative contract for every public module

Status: accepted (issue #37, M6.1). Normative for all later module ADRs and
implementation issues.

fortnum is primal-first. A routine computes a value, and that value is correct
on its own, before any derivative exists. But every public procedure is written
so a derivative can be added later without changing the primal signature or the
primal behavior. This document fixes the contract that makes "later" cheap: what
a derivative is called, which arguments carry it, and which policy class governs
how it is produced.

The contract is the shared vocabulary. A module ADR does not redefine these
terms; it cites the section numbers here and fills in the module-specific
choices. `CONTRIBUTING.md` summarizes the policy-class names for contributors
and names this file as the source of truth. Where the two disagree, this file
wins.

## 1. Policy classes

Each public procedure is assigned exactly one policy class. The class is
declared in the module-level doc comment and decides how the derivative is
produced, not whether one exists.

- `transparent`: Enzyme differentiates the implementation directly. The primal
  has no branch on an active variable and no hidden state, so source-level
  automatic differentiation gives the correct derivative with no extra code.
- `analytic_rule`: fortnum supplies handwritten JVP/VJP/gradient/HVP code. The
  closed form or recurrence is known and cheaper or more stable than
  differentiating the primal.
- `implicit_rule`: the derivative is defined by an implicit equation, not by
  differentiating the iteration that finds the primal. A root `f(x, p) = 0` and
  a linear solve `A(p) x = b(p)` both differentiate through the defining
  equation, not through the solver steps.
- `trace_rule`: the derivative is defined on a recorded, then frozen, primal
  trace. Adaptive methods choose step sizes, node counts, or refinement from the
  data; the derivative is taken at the schedule the primal selected, with that
  schedule held fixed.
- `primal_only`: no derivative is mathematically meaningful or supported. RNG
  draws are the canonical case: the output is not a differentiable function of
  the seed. A `primal_only` declaration carries a one-line justification.

The classes are exhaustive and disjoint. A procedure that seems to need two is
either two procedures, or its class is wrong.

## 2. Naming convention

For a primal routine `foo`, the public derivative entry points are:

- `foo_jvp(...)`: forward-mode product (Jacobian times vector).
- `foo_vjp(...)`: reverse-mode product (vector times Jacobian).
- `foo_grad(...)`: gradient, where the output is scalar and a gradient is what
  the caller wants.
- `foo_hvp(...)`: Hessian-vector product, where second order is meaningful.
- `foo_ad(...)`: optional combined convenience wrapper over the above.

Not every primal needs every name. A scalar-output routine offers `foo_grad`; a
vector-output routine offers `foo_jvp` and `foo_vjp`. Add the names the policy
class and the output shape justify, no more.

Enzyme-generated entry points stay internal. Public callers use the `foo_*`
names above and never call a raw `__enzyme_*` symbol. The internal symbols are
an implementation detail of the `analytic_rule` and `transparent` paths and may
change without notice.

## 3. Active argument rules

An argument is active if the derivative flows through it, inactive otherwise.
The classification is part of the contract and is declared per procedure.

- Active real and complex arrays use contiguous explicit-shape or assumed-size
  wrappers on the Enzyme path first. The simplest memory layout is the one the
  AD path is tested against; richer layouts are added only once they have a
  test.
- Descriptors, allocatables, polymorphism, optional active arguments, and
  derived-type components are allowed on the AD path only when a dedicated
  Enzyme test covers that shape. No test, no support; the primal may still use
  them freely.
- Integers, sizes, orders, keys, branch modes, RNG seeds, status flags, and
  workspace capacities are inactive. They select behavior or report it; they do
  not carry a derivative.

The `fortnum_status_t` object is inactive. Status reporting is a side channel,
not a differentiable output.

## 4. Per-module policy table

Default policy class for each planned module. The default is the expected case;
an individual procedure may differ when its doc comment says so and gives the
reason.

| Module      | Default policy   | Rationale |
|-------------|------------------|-----------|
| `special`   | `analytic_rule` (some `transparent`) | Derivatives follow known recurrences and identities (e.g. `Gamma'/Gamma = psi`); use the closed form where it is more stable, fall back to `transparent` for plain polynomial or rational kernels. The complex Bessel functions (`fortnum_special_complex_bessel`) and `hyperg_1f1_a1` use `analytic_rule` with the DLMF order/argument recurrences; the erf/erfc C-ABI provider (`fortnum_special_erf_cbind`) is `transparent`, active argument `x`. |
| `fft`       | `transparent`    | The transform is linear; its Jacobian is the transform itself, and the straight-line implementation differentiates cleanly. |
| `quadrature`| `transparent`    | Fixed-rule quadrature is a fixed weighted sum, so it is `transparent` with respect to integrand values; the nodes and weights are inactive parameters. `gauss_legendre` and `gauss_gen_laguerre` both build such a rule; the same `gauss_legendre_jvp`/`vjp`/`grad` products act on either weight vector, since only the integrand values are active. |
| `levin`     | `primal_only`    | Levin-u acceleration is a nonlinear rational transform of the term sequence and the selected order is data-dependent; no transparent or closed-form rule applies. Active: the terms. |
| `integrate` | `trace_rule`     | Adaptive integration picks subdivisions from the integrand; differentiate at the frozen subdivision the primal chose. |
| `ode`       | `trace_rule`     | Adaptive step control makes step sizes data-dependent; record the accepted step schedule and differentiate the frozen trace. The DOP853 integrator (`fortnum_ode_dop853`) shares the policy and the recorded-mesh carriers. |
| `roots`     | `implicit_rule`  | The root satisfies `f(x, p) = 0`; the implicit function theorem gives the derivative without differentiating the iteration. The n-dimensional solvers (`fortnum_multiroot`) follow the same rule with the Jacobian `J_x`; `deriv_central` and `argsort` in that module are `primal_only`. The complex region finder (`fortnum_roots_complex`, `complex_region_roots`) is `implicit_rule` with `differentiate_through=false`: a zero satisfies `f(z*, p) = 0`, so `dz*/dp = -f_p/f'(z*)` without differentiating the contour integral or the ZGGEV eigensolve, and no consumer differentiates through the region search, so it provides no JVP/VJP. |
| `interp`    | `transparent` (some `analytic_rule`; `grid_search` is `primal_only`) | Polynomial and spline evaluation is straight-line and `transparent`; a spline whose coefficients come from a solve uses `implicit_rule` or a hand-coded `analytic_rule` for the coefficient sensitivity. `grid_search` returns an integer cell index: it is `primal_only` because the index is inactive control flow (ad.md §3). Smooth derivatives of the interpolant are valid only inside a fixed cell; crossing a cell boundary is a non-smooth event and the caller must hold the index fixed when differentiating with respect to the evaluation point. B-spline evaluation (`fortnum_bspline`) is `transparent` inside a fixed knot span; its `bspline_span_index` is `primal_only`. |
| `rng`       | `primal_only`    | A pseudorandom draw is not a differentiable function of its seed; gradients of estimators built on draws live in the caller, not here. |

## 5. How module docs reference this contract

A module ADR or an implementation issue is normative downstream of this file. It
states, for each public procedure:

1. the policy class from section 1,
2. the active and inactive arguments per section 3,
3. the derivative entry-point names per section 2 that it provides.

It cites this document by section number rather than restating the definitions.
When a procedure departs from its module default in section 4, the doc names the
procedure and the reason. A reviewer checks the module doc against these
sections; a mismatch blocks the merge.

## 6. Adding derivatives without API churn

The primal ships first (M1 through M5); derivative products arrive later (issue
#40). The contract is built so that arrival changes no existing signature.

- Derivative routines are new public names under the section 2 convention. They
  are added beside the primal, never by adding arguments to it. A caller who
  only wants the value keeps calling `foo` unchanged.
- Primal routines stay `pure` where the algorithm allows and hold no
  module-level mutable state. Purity is the precondition for the `transparent`
  path and keeps the AD wrapper free to call the primal as many times as a
  product needs.
- The reserved derivative slots already exist. `fortnum_status_t` carries the
  status side channel that derivative routines reuse without a new type, and the
  oracle CSV format already reserves a derivative column behind the
  `has_derivative` header flag (`test/oracle/data/*.csv`). When `foo_jvp` lands,
  its oracle table is the same file with `has_derivative: 1` and the column
  filled; no format change, and tables written before then set the flag to `0`
  so consumers ignore the empty column.

Adding a derivative is therefore additive: new names, the same primal, the same
status type, the same CSV layout.
