# ADR: Optimizer-facing derivative API

Status: accepted (issue #41, M6.5). Normative downstream of `docs/design/ad.md`,
which it cites by section number and does not restate.

An optimizer does not care how a derivative was produced. It wants a value, a
gradient, a Jacobian-vector product, and a flag telling it whether to trust the
result. fortnum builds those products five different ways (`ad.md §1`): an
analytic recurrence, an implicit-function rule, a frozen adaptive trace, an
Enzyme-generated pass, or a finite-difference fallback. This layer hides that
choice behind one set of interfaces so a downstream code links against the call
shape, not against the backend.

Two pieces make it work: a flat active vector that carries heterogeneous inputs
as a single `x(:)`, and a small family of abstract interfaces that every
derivative product fits.

## 1. Flat active vector

An optimizer steps a single vector `x(:)`. A physics kernel reads boundary
modes, radial profiles, and shape parameters. The flat active vector reconciles
the two. The caller declares named blocks once, packs each quantity into its
slice, and unpacks it back inside the kernel.

`fortnum_active_layout_t` records the total length `n` and the named blocks that
tile `[1, n]`. A block is a name plus a 1-based offset and a size. Blocks are
appended in declaration order; each new block starts where the previous one
ended, so the layout tiles the vector without gaps or overlap.

```fortran
type(fortnum_active_layout_t) :: layout
type(fortnum_status_t)        :: st
real(dp)                      :: x(:), modes(8), profile(33)

call layout_init(layout, 2)
call layout_add(layout, "boundary_modes", 8, st)
call layout_add(layout, "pressure_profile", 33, st)
! layout%n == 41
call pack_block(layout, x, "boundary_modes", modes, st)
call unpack_block(layout, x, "pressure_profile", profile, st)
```

The layout is data the caller owns. `layout_init`, `layout_add`, `pack_block`,
and `unpack_block` are `pure` and hold no module-level state, so two optimizer
threads packing two vectors never collide. Pack and unpack validate the vector
length and the block size at the boundary (`ad.md §3`: status is a side
channel, never a differentiable output) and report `FORTNUM_DOMAIN_ERROR` on a
mismatch.

The same layout serves the seed and the gradient. `x_dot` in a JVP, `x_bar`
from a VJP, and `g` from a gradient all share the layout of `x`, so the
optimizer reads a gradient component by the same name it used to pack the input.

## 2. Derivative-product interfaces

`fortnum_ad_interfaces` declares five abstract interfaces. A kernel implements
the ones its output shape justifies (`ad.md §2`): a scalar objective offers
`grad_fn` and `hvp_fn`; a vector residual offers `jvp_fn` and `vjp_fn`; both
offer `value_fn`.

```fortran
subroutine value_fn(n, x, y, context, status)
subroutine jvp_fn  (n, x, x_dot, y, y_dot, context, status)
subroutine vjp_fn  (n, x, y_bar, x_bar, context, status)
subroutine grad_fn (n, x, f, g, context, status)
subroutine hvp_fn  (n, x, v, f, hv, context, status)
```

`n` is the flat length. The real arrays are contiguous explicit-shape, the
layout the Enzyme path is tested against first (`ad.md §3`). The forward product
returns both the primal `y` and the directional derivative `y_dot = J(x) x_dot`
in one call, so a line search reuses the value it already paid for. The reverse
product returns `x_bar = J(x)^T y_bar`. `grad_fn` and `hvp_fn` are the
scalar-objective specializations an unconstrained optimizer calls directly.

### Context without globals

`context` is `class(*)`. The kernel carries its configuration and workspace
through it: a derived type holding the fixed grid, the frozen trace, the
profiles a residual closes over. The optimizer constructs the context once and
threads it through every call. Nothing lives at module scope, so a kernel is
reentrant and safe under parallel objective evaluations. This is the Fortran
substitute for a closure: state travels in an argument, not in a global pointer.

## 3. Backend opacity

The product is opaque; its provenance is metadata. `fortnum_ad_status_t` wraps
the primal `fortnum_status_t` and adds two tags.

`backend` records which machinery ran: `ANALYTIC`, `IMPLICIT`, `TRACE`,
`GENERATED` (Enzyme), or `FINITE_DIFF`. The first four mirror the policy classes
of `ad.md §1`; the fifth is the fallback a kernel uses when no exact product
exists.

`quality` records how good the product is, independent of the backend. `EXACT`
is correct to rounding. `APPROXIMATE` carries a controlled truncation error: a
finite difference, or a frozen trace evaluated at a point where the true
adaptive schedule would differ. `NONSMOOTH` marks a branch or event boundary
(`ad.md §3`), where the derivative is undefined and the optimizer must reject
the step.

`ad_status_ok` returns `.true.` only when the primal status is OK and the
quality is not `NONSMOOTH`. An optimizer gates a step on that single predicate.
A caller that reads only the embedded `fortnum_status_t` still sees the primal
error code, so existing status checks keep working.

The split between `backend` and `quality` is deliberate. A trust threshold keys
on quality, not backend: an `ANALYTIC` product and a `GENERATED` product that
are both `EXACT` get the same trust. The backend tag is for logging and for
diagnosing why a product is only `APPROXIMATE`.

## 4. DESC parity

DESC computes stellarator equilibria in Python and gets gradients from JAX. The
physics is rewritten in a framework that can differentiate it. That rewrite is
the cost: every residual, every objective, every constraint lives in JAX-traceable
code, and a change to the physics is a change to the differentiable model.

fortnum inverts the dependency. The physics stays in Fortran. A residual or an
objective ships as a primal `foo` plus its derivative products `foo_jvp`,
`foo_vjp`, `foo_grad` (`ad.md §2`), each produced by whichever backend its policy
class dictates. The optimizer-facing layer exposes those products through the
interfaces above. A stellarator-optimization code packs its boundary modes and
profiles into the flat active vector, calls `grad_fn` for the objective and
`jvp_fn`/`vjp_fn` for the residual block, and reads `quality` to gate each step.
It never learns whether a gradient came from a hand-coded rule or an Enzyme pass,
and it never reimplements the physics in another language.

The parity claim is concrete: a smooth residual and a smooth objective, exposed
through this layer, give an outer optimizer the same value-and-gradient contract
JAX gives DESC, with the physics kept in the language it was written in.
`test/ad/test_optimizer_api.f90` is the worked check. It drives a scalar
objective through `grad_fn` and verifies the gradient against central finite
difference, then drives a vector residual through `jvp_fn`/`vjp_fn` and verifies
the dot-product adjoint identity `u . (J v) = v . (J^T u)`. Both run without
Enzyme, so the contract holds on the default build.
