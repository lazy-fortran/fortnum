# ADR: Downstream active-kernel pattern

Status: accepted (issue #42, M6.6). Normative downstream of `docs/design/ad.md`
and `docs/design/optimizer_api.md`, which it cites by section number.

Plasma codes typically expose a large primal driver: read data, classify
orbits, run a solver loop, write output. None of that is differentiable end
to end. Derivatives are meaningful only inside a smooth subregion: a kernel
that maps floating-point parameters to a residual or an objective with no
hard branches, file I/O, or side-channel calls inside the active region.

This document fixes the pattern that every fortnum downstream code follows
to extract that smooth kernel, wire it to the optimizer-facing API (#41),
and keep the rest of the code untouched.

## 1. Identifying the active region

The active region is a pure map: real inputs to real outputs, smooth on the
domain the optimizer visits. Anything outside must stay outside.

Leave these OUTSIDE the active derivative region:

- solver loops and iteration management (bisection, Newton, mesh refinement)
- file I/O and restart logic
- MPI communication and global reductions
- orbit classification: loss/confined/stochastic labels are integer outputs;
  they are inactive (`ad.md §3`)
- event detection: a particle hitting a wall or crossing a separatrix is a
  branch boundary; the derivative is undefined at that point
- retry and fallback logic; these change the trace and break `trace_rule`
  smoothness

The active region starts where floating-point parameters enter a smooth
residual or objective computation. It ends where the result is returned to
the caller.

## 2. The pack/write contract

Each optimization variable block has a physical meaning: boundary modes,
a radial profile, shape parameters, transport coefficients. The optimizer
sees a single flat vector `x(n)`. The kernel reads the blocks it needs.

```fortran
type(fortnum_active_layout_t) :: layout
real(dp)                      :: x(n)

call layout_init(layout, nblock)
call layout_add(layout, "boundary_modes", 8, st)
call layout_add(layout, "pressure_profile", 33, st)

call pack_block(layout, x, "boundary_modes", modes, st)
call pack_block(layout, x, "pressure_profile", profile, st)
```

Inside the kernel, read the slice you need:

```fortran
call unpack_block(layout, x, "pressure_profile", profile, st)
```

The layout is data the caller owns (`optimizer_api.md §1`). It is `pure` and
holds no module-level state, so two optimizer threads packing two separate
vectors never collide.

## 3. Wiring to the abstract interfaces

A smooth kernel implements the derivative products its output shape justifies
(`ad.md §2`):

- scalar objective: `grad_fn` and optionally `hvp_fn`
- vector residual: `jvp_fn` and `vjp_fn`
- both: `value_fn`

Each procedure must match the abstract interface from `fortnum_ad_interfaces`
(`optimizer_api.md §2`). The signature is:

```fortran
subroutine grad_fn(n, x, f, g, context, status)
    integer,                   intent(in)    :: n
    real(dp),                  intent(in)    :: x(n)
    real(dp),                  intent(out)   :: f
    real(dp),                  intent(out)   :: g(n)
    class(*),                  intent(inout) :: context
    type(fortnum_ad_status_t), intent(out)   :: status
end subroutine grad_fn
```

The kernel accesses its configuration and workspace through `context`, a
`class(*)` the optimizer constructs once and threads through every call
(`optimizer_api.md §2`). Nothing lives at module scope.

## 4. Context without globals

A real physics kernel closes over grids, profiles, precomputed matrices, and
the frozen trace from a prior primal call. All of that goes into a derived type
that the caller allocates once and hands in as `context`:

```fortran
type, extends(base_context_t) :: my_kernel_ctx_t
    type(fortnum_active_layout_t) :: layout
    real(dp), allocatable         :: grid(:)
    type(integrate_result_t)      :: frozen_trace
end type my_kernel_ctx_t
```

Inside the kernel body, select-type gives access:

```fortran
select type (ctx => context)
type is (my_kernel_ctx_t)
    call unpack_block(ctx%layout, x, "alpha", alpha, us)
end select
```

No module pointer. No global procedure variable. The optimizer evaluates
multiple objectives in parallel by handing each call its own `context`.

## 5. Calling fortnum derivative products inside the kernel

The derivative entry points follow `ad.md §2`. Use the module product
directly rather than routing through the abstract interface from the call
site; the abstract interface is for the outer optimizer layer.

For a Bessel-based residual component at order `n`:

```fortran
call bessel_in_jvp(n, x_val, v_val, jv_comp)   ! analytic_rule, exact
```

For an adaptive-integration component with a frozen trace:

```fortran
call integrate_qag_jvp(dfdp, result, di_dp, st, ctx)  ! trace_rule
```

Check `ad_status_ok(st)` after each call. The optimizer gates a step on that
single predicate (`optimizer_api.md §3`): `NONSMOOTH` means the perturbation
crossed a branch and the product is unusable.

## 6. Repo-specific guidance

### libneo

Start with the math and numerical routines: Bessel functions, the lower
incomplete gamma, FFT wrappers, ODE integration, and QUADPACK. These are the
cleanest kernels and have derivative products NOW (see `docs/migration_libneo_ad.md`).

Coordinate-system routines (Boozer, Hamada, PEST transformations) come later
and only where the map from configuration parameters to transformed coordinates
is smooth and isolated. The full coordinate pipeline is not a candidate:
it reads mesh files, applies classification branches, and calls external
equilibrium solvers.

### NEO-2

NEO-2 calls an external numerical backend for special functions, quadrature,
ODE integration, B-splines, and RNG. Replace those calls with the corresponding
fortnum APIs (see `docs/migration_libneo_ad.md`). The
derivative-ready wrappers arrive as drop-in substitutes: the primal signatures
match the libneo mapping table, and the derivative entry points are new names
that leave the primal calls unchanged.

NEO-2 transport kernels call `odeint` and QUADPACK inside a smooth sensitivity
calculation. Those become `ode_integrate` + `ode_integrate_jvp` and
`integrate_qag` + `integrate_qag_jvp` under `trace_rule`. The orbit
integration that terminates on a loss boundary stays primal-only; the loss
boundary is a non-smooth event.

### KAMEL

KAMEL uses an external special-function and quadrature backend, libcerf,
SLATEC/AMOS (Bessel via Fortran 77 wrappers), and a custom ODE integrator.
Replace the following where the calling context is smooth:

- special functions: fortnum_special (Bessel, gamma, Dawson)
- quadrature: fortnum_integrate (QAG/QAGS pattern)
- libcerf Dawson: `dawson` + `dawson_jvp`/`dawson_grad`
- SLATEC/AMOS Bessel: `bessel_in`/`bessel_kn` + `bessel_in_jvp`/`bessel_kn_jvp`
- ODE integration: `ode_integrate` + `ode_integrate_jvp`/`ode_integrate_vjp`

Dielectric tensor assembly and the dispersion relation evaluation are the
target kernels: smooth maps from plasma parameters to a complex tensor, then
to a determinant. Pack plasma parameters (density, temperature, drift
velocities, magnetic field strength) into the active vector; call the tensor
assembly as the smooth kernel; the determinant is a scalar objective.

### GORILLA

GORILLA integrates guiding-centre orbits with a geometric integration scheme.
The orbit integration itself is not the derivative target; the orbit is a
primal tool, not a smooth loss function.

Use fortnum for generic RK and interpolation helpers that appear inside
GORILLA's infrastructure:

- Cash-Karp RK steps that advance the orbit: `ode_integrate` as a drop-in
  where the same step structure is used outside the orbit physics
- Lagrange interpolation of field components on the tetrahedral mesh:
  `lagrange_weights` + `lagrange_weights_jvp`/`lagrange_fval_jvp`

The orbit tracing, the tetrahedron decomposition, and the field mapping
physics stay in domain code and remain primal-only.

### SIMPLE

SIMPLE draws guiding-centre Monte Carlo orbits and classifies each as
confined or lost. The loss/confinement label is a hard output: non-differentiable
and primal-only. Monte Carlo estimators of confinement time are also
primal-only because their gradient lives in the caller, not here
(`ad.md §4`, `rng` entry).

Use fortnum for the following where they appear inside a smooth SIMPLE sub-kernel:

- The RNG stream: `fortnum_rng` (primal-only; no derivative ever)
- Generic ODE helpers used for orbit initialization or field-line following
  outside the loss-classification path: `ode_integrate`
- Any quadrature that appears in a smooth neoclassical weight: `integrate_qag`

The derivative target in SIMPLE, if one exists, is a smooth objective built
from the confined-orbit distribution BEFORE the binary classification step.
That objective is a Monte Carlo estimator; its gradient is the caller's
responsibility, not SIMPLE's.
