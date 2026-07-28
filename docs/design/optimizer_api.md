# Optimizer-facing API

Status: implemented.

Optimizers consume values and derivative products without depending on the
candidate that produced them.

## Active-vector layout

`fortnum_active_layout_t` maps named blocks to one flat optimizer vector.

```fortran
type(fortnum_active_layout_t) :: layout

call layout_init(layout, 2)
call layout_add(layout, "boundary", 8, status)
call layout_add(layout, "profile", 33, status)
call pack_block(layout, x, "boundary", boundary, status)
call unpack_block(layout, x, "profile", profile, status)
```

Blocks tile the vector in declaration order. The same layout applies to primal
inputs, tangent directions, gradients, and cotangents. Shape and name errors
return status. The caller owns the layout and vectors.

## Product callbacks

`fortnum_ad_interfaces` declares:

```fortran
subroutine value_fn(n, x, y, context, status)
subroutine jvp_fn(n, x, x_dot, y, y_dot, context, status)
subroutine vjp_fn(n, x, y_bar, x_bar, context, status)
subroutine grad_fn(n, x, f, g, context, status)
subroutine hvp_fn(n, x, v, f, hv, context, status)
```

The JVP returns \(y\) and \(Jx_{\rm dot}\). The VJP returns \(J^Ty_{\rm bar}\).
Gradient and HVP callbacks specialize scalar objectives.

`context` carries immutable configuration and caller-owned workspace. It
replaces global callback state and permits concurrent objective evaluations.

## HVP consumer

The first second-order consumer is a Hessian-free Newton-CG step for the
nonlinear least-squares application in
`benchmark/bench_nonlinear_least_squares.f90`. Its inner conjugate-gradient
solve requires repeated \(H(x)v\) products to approximate the step
\(Hs=-g\). It does not require a materialized Hessian.

Candidate HVPs must return the objective and \(H(x)v\) through `hvp_fn`.
Forward-over-reverse autodiff, reverse-over-forward autodiff, and an
analytical least-squares HVP are admissible. Selection uses the complete
Newton-CG workload: validation, achieved residual reduction, wall clock, and
peak memory. A faster isolated HVP does not win if it increases outer or inner
iterations.

## Candidate opacity

The callback shape does not expose `autodiff`, `analytical`, or `hybrid`
implementation details. `fortnum_ad_status_t` reports provenance and quality
for logging and step acceptance.

A candidate registry can choose an implementation before the optimizer loop
from validated workload metadata. A mechanism name alone never selects the
callback.

## Validation

`test/ad/test_optimizer_api.f90` checks a scalar gradient with central finite
differences and a vector residual with the adjoint identity

\[
u^T(Jv)=v^T(J^Tu).
\]
