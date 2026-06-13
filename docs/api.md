# fortnum public API

Every routine listed here is thread-safe. State is explicit and caller-owned;
no module-level mutable state exists anywhere in fortnum. Two threads that hold
separate state objects never race.

Derivative products (`foo_jvp`, `foo_vjp`, `foo_grad`, `foo_hvp`) shipped in
issue #40 as additive routines: the primal signatures are unchanged. The
[Derivative products](#derivative-products) section lists them per module.

Derivative policy vocabulary follows `docs/design/ad.md §1`:
`transparent`, `analytic_rule`, `implicit_rule`, `trace_rule`, `primal_only`.

---

## fortnum_kinds

Kind parameters re-exported from `iso_fortran_env`.

| Symbol | Kind   | Description |
|--------|--------|-------------|
| `dp`   | `real64`  | Double precision; the primary working kind. |
| `sp`   | `real32`  | Single precision for mixed-precision interfaces. |
| `i4`   | `int32`   | 32-bit integer. |
| `i8`   | `int64`   | 64-bit integer for large index spaces. |

---

## fortnum_status

Error reporting. The `fortnum_status_t` object is inactive in the derivative
sense (ad.md §3): it is a side channel, not a differentiable output.

### Type: `fortnum_status_t`

```fortran
type :: fortnum_status_t
    integer            :: code = FORTNUM_OK
    character(120)     :: msg  = ""
end type fortnum_status_t
```

| Constant                    | Value | Meaning |
|-----------------------------|-------|---------|
| `FORTNUM_OK`                | 0     | Success. |
| `FORTNUM_DOMAIN_ERROR`      | 1     | Argument out of domain. |
| `FORTNUM_CONVERGENCE_ERROR` | 2     | Iteration did not converge. |
| `FORTNUM_NOT_IMPLEMENTED`   | 3     | Capability not yet present. |

### `status_ok(s) -> logical`

```fortran
pure logical function status_ok(s)
    type(fortnum_status_t), intent(in) :: s
```

Returns `.true.` when `s%code == FORTNUM_OK`.

### `status_set(s, code, msg)`

```fortran
pure subroutine status_set(s, code, msg)
    type(fortnum_status_t), intent(out) :: s
    integer,                intent(in)  :: code
    character(*),           intent(in)  :: msg
```

Sets both fields of an existing status object.

---

## fortnum_special

Umbrella re-export of the special function families. Import from here or from
the sub-modules directly.

Default derivative policy: `analytic_rule` (ad.md §4).
Derivative products reserved for #40: named per sub-module below.

### fortnum_special_bessel

Policy: `analytic_rule`. Active argument: `x`. Inactive: `n` (order selector).
Reserved: `bessel_in_grad`, `bessel_kn_grad`.

#### `bessel_in(n, x) -> real(dp)` (elemental)

```fortran
elemental function bessel_in(n, x) result(value)
    integer,  intent(in) :: n
    real(dp), intent(in) :: x
```

Modified Bessel function I_n(x) of the first kind, integer order. Handles
negative `n` (I_{-n} = I_n) and negative `x`.

#### `bessel_in_array(nmax, x, values)`

```fortran
pure subroutine bessel_in_array(nmax, x, values)
    integer,  intent(in)  :: nmax
    real(dp), intent(in)  :: x
    real(dp), intent(out) :: values(0:nmax)
```

Fills `values(0:nmax)` with I_0(x) through I_nmax(x) in one pass. More
efficient than `nmax+1` separate calls to `bessel_in`.

#### `bessel_kn(n, x) -> real(dp)` (pure function)

```fortran
pure function bessel_kn(n, x) result(kn)
    integer,  intent(in) :: n
    real(dp), intent(in) :: x
```

Modified Bessel function K_n(x) of the second kind, integer order `n >= 0`,
`x > 0`. K_{-n} = K_n by symmetry.

### fortnum_special_dawson

Policy: `analytic_rule`. Active argument: `x`. Reserved: `dawson_grad`.
Analytic rule: F'(x) = 1 - 2x F(x).

#### `dawson(x) -> real(dp)` (elemental)

```fortran
elemental function dawson(x) result(f)
    real(dp), intent(in) :: x
```

Dawson integral F(x) = exp(-x^2) integral_0^x exp(t^2) dt.

### fortnum_special_gamma

Policy: `analytic_rule`. Active arguments: `a` (shape), `x` (limit).
Reserved: `gamma_lower_grad`, `gamma_reg_p_grad`.

#### `gamma_lower(a, x) -> real(dp)` (pure function)

```fortran
pure function gamma_lower(a, x) result(g)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: x
```

Unnormalized lower incomplete gamma: gamma(a, x) = P(a,x) * Gamma(a). Requires
`a > 0`, `x >= 0`. Stops with `error stop` on domain violation.

#### `gamma_reg_p(a, x) -> real(dp)` (pure function)

```fortran
pure function gamma_reg_p(a, x) result(p)
    real(dp), intent(in) :: a
    real(dp), intent(in) :: x
```

Regularized lower incomplete gamma P(a, x) = gamma_lower(a, x) / Gamma(a).
Requires `a > 0`, `x >= 0`.

---

## fortnum_fft

1D discrete Fourier transforms. Mixed-radix (2,3,4,5) Stockham passes; Bluestein
chirp-z fallback for lengths with other prime factors.

Convention: the forward transform applies exp(-2 pi i j k / n), unnormalized.
The inverse applies exp(+2 pi i j k / n) without the 1/n factor.

Default derivative policy: `transparent` (ad.md §4). The DFT is linear; Enzyme
differentiates the implementation directly. Active: `z` in `fft_c2c`, `x` and
`c` in `fft_r2c`. Inactive: `sign`, plan fields, all integers.
Reserved derivative products for #40: `fft_c2c_jvp`, `fft_c2c_vjp`,
`fft_r2c_jvp`, `fft_r2c_vjp`.

### Type: `fortnum_fft_plan_t`

Caller-owned reusable plan. Holds twiddle tables for a fixed transform length.
Read-only during the transform; safe to use from multiple threads on the same
length without copying.

### `fft_plan_init(plan, n)`

```fortran
subroutine fft_plan_init(plan, n)
    type(fortnum_fft_plan_t), intent(out) :: plan
    integer,                  intent(in)  :: n
```

Initialize a plan for transforms of length `n`. Allocates twiddle tables.
Stops with `error stop` if `n < 1`.

### `fft_r2c(x, c [, plan])`

```fortran
subroutine fft_r2c(x, c, plan)
    real(dp),                            intent(in)  :: x(:)
    complex(dp),                         intent(out) :: c(:)
    type(fortnum_fft_plan_t), optional,  intent(in)  :: plan
```

Real-to-complex DFT. `c` must have size `size(x)/2 + 1`; `c(k+1)` is the k-th
bin. If `plan` is supplied it must have been built for `size(x)`.

### `fft_c2c(z, sign)`

```fortran
subroutine fft_c2c(z, sign)
    complex(dp), intent(inout) :: z(:)
    integer,     intent(in)    :: sign
```

In-place complex-to-complex DFT. `sign = -1`: forward; `sign = +1`: inverse
(unnormalized). Stops with `error stop` if `sign` is not -1 or +1.

---

## fortnum_quadrature

Fixed-rule Gauss-Legendre quadrature: node and weight generation.

Default derivative policy: `transparent` (ad.md §4). Fixed-rule quadrature is a
fixed weighted sum; Enzyme differentiates the sum directly. Active: integrand
values supplied by the caller. Inactive: `n`, `a`, `b`.

### `gauss_legendre(n, x, w)`

```fortran
pure subroutine gauss_legendre(n, x, w)
    integer,  intent(in)  :: n
    real(dp), intent(out) :: x(n), w(n)
```

Gauss-Legendre nodes and weights on [-1, 1]. Nodes in ascending order. Uses
Newton iteration from asymptotic initial estimates; exploits symmetry so only
ceil(n/2) iterations run. Requires `n >= 1`.

### `gauss_legendre_ab(n, a, b, x, w)`

```fortran
subroutine gauss_legendre_ab(n, a, b, x, w)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: a, b
    real(dp), intent(out) :: x(n), w(n)
```

Gauss-Legendre rule mapped to [a, b]. The quadrature sum
`sum_i w(i)*f(x(i))` approximates the integral of `f` over [a, b].

---

## fortnum_integrate_gk

Single-panel and simple adaptive Gauss-Kronrod quadrature. The globally adaptive
driver in `fortnum_integrate` reuses `gk_apply` per subinterval; this module
handles the simple adaptive case directly.

Default derivative policy: `trace_rule` (ad.md §4).

### Abstract interface: `gk_integrand_t`

```fortran
abstract interface
    function gk_integrand_t(x) result(fx)
        real(dp), intent(in) :: x
        real(dp)             :: fx
    end function gk_integrand_t
end interface
```

### `integrate_gk(f, a, b, epsabs, epsrel, result, abserr, ierr [, key, limit])`

```fortran
subroutine integrate_gk(f, a, b, epsabs, epsrel, result, abserr, ierr, key, limit)
    procedure(gk_integrand_t)        :: f
    real(dp), intent(in)             :: a, b, epsabs, epsrel
    real(dp), intent(out)            :: result, abserr
    integer,  intent(out)            :: ierr
    integer,  intent(in), optional   :: key, limit
```

Integrate `f` over [a, b]. `key` selects the GK pair: 15, 21 (default), 31, or
61. `limit` bounds adaptive subdivisions (default 200). `ierr`: 0 converged, 1
limit reached, 2 roundoff detected, 3 bad integrand behaviour, 6 invalid input.

### `gk_apply(f, key, a, b, result, abserr, resabs, resasc)`

```fortran
subroutine gk_apply(f, key, a, b, result, abserr, resabs, resasc)
    procedure(gk_integrand_t)  :: f
    integer,  intent(in)       :: key
    real(dp), intent(in)       :: a, b
    real(dp), intent(out)      :: result, abserr, resabs, resasc
```

Apply one GK panel over [a, b]. Used internally by `fortnum_integrate` and
available for callers who manage subdivisions themselves.

---

## fortnum_integrate

Globally adaptive integration: QAG, QAGS, QAGP, QAGIU (QUADPACK pattern).
Caller-owned workspace; no allocation inside the bisection loop.

Default derivative policy: `trace_rule` (ad.md §4). The subdivision is
data-dependent; derivatives are taken at the frozen accepted subdivision recorded
in `integrate_result_t`.

Reserved derivative products for #40: `integrate_qag_jvp`, `integrate_qag_vjp`,
`integrate_qags_jvp`, `integrate_qags_vjp`, `integrate_qag_grad`.

### Abstract interface: `integrate_integrand_t`

```fortran
abstract interface
    function integrate_integrand_t(x, ctx) result(fx)
        real(dp), intent(in)           :: x
        class(*), intent(in), optional :: ctx
        real(dp)                       :: fx
    end function integrate_integrand_t
end interface
```

`ctx` is the only channel for parameters; pass a derived type and select its
type inside the integrand. `ctx` is `intent(in)`.

### Types

| Type | Purpose |
|------|---------|
| `integrate_workspace_t` | QUADPACK alist/blist/rlist/elist work stack. Caller-owned; reuse across calls. |
| `integrate_epstab_t`    | Wynn epsilon table for QAGS/QAGP/QAGIU. Default-initialized value is a valid empty table. |
| `integrate_result_t`    | Result value with its error estimate, status, plus the frozen accepted subdivision (trace_rule schedule). |

### `integrate_qag(f, a, b, epsabs, epsrel, workspace, result, status [, key, limit, ctx])`

```fortran
subroutine integrate_qag(f, a, b, epsabs, epsrel, workspace, result, status, key, limit, ctx)
    procedure(integrate_integrand_t)          :: f
    real(dp),                    intent(in)   :: a, b, epsabs, epsrel
    type(integrate_workspace_t), intent(inout):: workspace
    type(integrate_result_t),    intent(inout):: result
    type(fortnum_status_t),      intent(out)  :: status
    integer,  intent(in), optional            :: key, limit
    class(*), intent(in), optional            :: ctx
```

Globally adaptive integration on a finite interval. No extrapolation. `key`
defaults to 21; `limit` defaults to 500.

### `integrate_qags(f, a, b, epsabs, epsrel, workspace, epstab, result, status [, limit, ctx])`

```fortran
subroutine integrate_qags(f, a, b, epsabs, epsrel, workspace, epstab, result, status, limit, ctx)
    procedure(integrate_integrand_t)          :: f
    real(dp),                    intent(in)   :: a, b, epsabs, epsrel
    type(integrate_workspace_t), intent(inout):: workspace
    type(integrate_epstab_t),    intent(inout):: epstab
    type(integrate_result_t),    intent(inout):: result
    type(fortnum_status_t),      intent(out)  :: status
    integer,  intent(in), optional            :: limit
    class(*), intent(in), optional            :: ctx
```

Adaptive bisection plus Wynn epsilon extrapolation. GK21 throughout. For
endpoint singularities and integrable algebraic or logarithmic spikes.

### `integrate_qagp(f, a, b, points, epsabs, epsrel, workspace, epstab, result, status [, limit, ctx])`

```fortran
subroutine integrate_qagp(f, a, b, points, epsabs, epsrel, workspace, epstab, result, status, limit, ctx)
    procedure(integrate_integrand_t)          :: f
    real(dp),                    intent(in)   :: a, b, epsabs, epsrel
    real(dp),                    intent(in)   :: points(:)
    type(integrate_workspace_t), intent(inout):: workspace
    type(integrate_epstab_t),    intent(inout):: epstab
    type(integrate_result_t),    intent(inout):: result
    type(fortnum_status_t),      intent(out)  :: status
    integer,  intent(in), optional            :: limit
    class(*), intent(in), optional            :: ctx
```

QAGS seeded with user break points at known interior singularities or kinks.
Break points coinciding with `a` or `b` within panel tolerance return
`FORTNUM_DOMAIN_ERROR`.

### `integrate_qagiu(f, bound, inf, epsabs, epsrel, workspace, epstab, result, status [, limit, ctx])`

```fortran
subroutine integrate_qagiu(f, bound, inf, epsabs, epsrel, workspace, epstab, result, status, limit, ctx)
    procedure(integrate_integrand_t)          :: f
    real(dp),                    intent(in)   :: bound, epsabs, epsrel
    integer,                     intent(in)   :: inf
    type(integrate_workspace_t), intent(inout):: workspace
    type(integrate_epstab_t),    intent(inout):: epstab
    type(integrate_result_t),    intent(inout):: result
    type(fortnum_status_t),      intent(out)  :: status
    integer,  intent(in), optional            :: limit
    class(*), intent(in), optional            :: ctx
```

Semi-infinite or doubly infinite interval. `inf = +1`: [bound, +inf);
`inf = -1`: (-inf, bound]; `inf = +2`: (-inf, +inf) split at `bound`. The
QUADPACK dqagi transform maps the interval to (0,1] and drives QAGS.

### `integrate(f, a, b, value, status [, epsabs, epsrel, key, ctx])`

```fortran
subroutine integrate(f, a, b, value, status, epsabs, epsrel, key, ctx)
    procedure(integrate_integrand_t)     :: f
    real(dp),               intent(in)   :: a, b
    real(dp),               intent(out)  :: value
    type(fortnum_status_t), intent(out)  :: status
    real(dp), intent(in), optional       :: epsabs, epsrel
    integer,  intent(in), optional       :: key
    class(*), intent(in), optional       :: ctx
```

Flat convenience call. Owns its workspace; discards the trace. Defaults:
`epsabs = 0`, `epsrel = 1e-8`.

---

## fortnum_ode

Adaptive ODE integration with Cash-Karp RK5(4) and PI step-size control.

Default derivative policy: `trace_rule` (ad.md §4). The accepted step schedule
is data-dependent; a sensitivity differentiates the frozen recorded mesh in
`ode_solution_t` with that schedule held fixed.

Reserved derivative products for #40: `ode_integrate_jvp`, `ode_integrate_vjp`,
`ode_at_jvp`, `ode_at_vjp`.

Event direction constants:

| Constant | Value | Meaning |
|----------|-------|---------|
| `ODE_EVENT_RISING`  | +1 | Detect g crossing zero from below. |
| `ODE_EVENT_FALLING` | -1 | Detect g crossing zero from above. |
| `ODE_EVENT_ANY`     |  0 | Detect either direction. |

### Abstract interfaces

```fortran
abstract interface
    subroutine ode_rhs_t(t, y, dydt, ctx)
        real(dp), intent(in)           :: t
        real(dp), intent(in)           :: y(:)
        real(dp), intent(out)          :: dydt(:)
        class(*), intent(in), optional :: ctx
    end subroutine ode_rhs_t
end interface

abstract interface
    function ode_event_t(t, y, ctx) result(g)
        real(dp), intent(in)           :: t
        real(dp), intent(in)           :: y(:)
        class(*), intent(in), optional :: ctx
        real(dp)                       :: g
    end function ode_event_t
end interface
```

### Types

| Type | Purpose |
|------|---------|
| `ode_problem_t`   | Problem definition: rhs pointer, t0, t1, y0, tolerances, step bounds, event pointer. |
| `ode_workspace_t` | Caller-owned stage arrays (k1..k6, temporaries). Reuse across calls to avoid reallocation. |
| `ode_solution_t`  | Recorded accepted-step mesh (t, y, h, err) plus event fields and status. |

### `ode_integrate(problem, workspace, solution, status)`

```fortran
subroutine ode_integrate(problem, workspace, solution, status)
    type(ode_problem_t),    intent(in)    :: problem
    type(ode_workspace_t),  intent(inout) :: workspace
    type(ode_solution_t),   intent(inout) :: solution
    type(fortnum_status_t), intent(out)   :: status
```

Integrate from `problem%t0` to `problem%t1`. Fills `solution` with the
accepted-step mesh. `FORTNUM_CONVERGENCE_ERROR` on `max_steps` exceeded or step
forced below `hmin`. `FORTNUM_DOMAIN_ERROR` on invalid input.

### `ode_solve(rhs, t0, t1, y0, t_out, y_out, status [, rtol, atol])`

```fortran
subroutine ode_solve(rhs, t0, t1, y0, t_out, y_out, status, rtol, atol)
    procedure(ode_rhs_t)               :: rhs
    real(dp),               intent(in) :: t0, t1
    real(dp),               intent(in) :: y0(:)
    real(dp), allocatable, intent(out) :: t_out(:)
    real(dp), allocatable, intent(out) :: y_out(:,:)
    type(fortnum_status_t), intent(out):: status
    real(dp), intent(in), optional     :: rtol, atol
```

Flat convenience call. Allocates and returns the mesh; discards the workspace.
Events are not exposed. Tolerances default to `rtol = 1e-6`, `atol = 1e-9`.

### `ode_at(problem, t_eval, workspace, y_out, status)` (from `fortnum_ode_wrapper`)

```fortran
subroutine ode_at(problem, t_eval, workspace, y_out, status)
    type(ode_problem_t),   intent(in)    :: problem
    real(dp),              intent(in)    :: t_eval(:)
    type(ode_workspace_t), intent(inout) :: workspace
    real(dp), allocatable, intent(out)   :: y_out(:,:)
    type(fortnum_status_t), intent(out)  :: status
```

Integrate and return the solution at each time in `t_eval`. `t_eval` must be
strictly monotone. `y_out` is allocated `(neq, size(t_eval))`. Policy:
`trace_rule` per segment.

---

## fortnum_roots

Scalar root-finding.

Default derivative policy: `implicit_rule` (ad.md §4). The root satisfies
f(x*, p) = 0; the implicit function theorem gives dx*/dp = -f_p/f_x without
differentiating the iteration.

Reserved derivative products for #40: `root_bisect_jvp`, `root_bisect_vjp`,
`root_newton_jvp`, `root_newton_vjp`, `root_brent_jvp`, `root_brent_vjp`.

### Abstract interfaces

```fortran
abstract interface
    pure function root_fn_t(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
    end function root_fn_t
end interface

abstract interface
    pure subroutine root_fn_df_t(x, fx, dfx)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: fx, dfx
    end subroutine root_fn_df_t
end interface
```

### `root_bisect(f, a, b, x, status [, xtol, ftol, max_iter])`

```fortran
subroutine root_bisect(f, a, b, x, status, xtol, ftol, max_iter)
    procedure(root_fn_t)               :: f
    real(dp),               intent(in) :: a, b
    real(dp),               intent(out):: x
    type(fortnum_status_t), intent(out):: status
    real(dp), intent(in), optional     :: xtol, ftol
    integer,  intent(in), optional     :: max_iter
```

Illinois-variant bisection on [a, b]. Requires f(a)*f(b) < 0. Defaults:
`xtol = 4*epsilon(1._dp)`, `ftol = 0`, `max_iter = 200`.

### `root_newton(fdf, a, b, x0, x, status [, xtol, ftol, max_iter, deriv_floor])`

```fortran
subroutine root_newton(fdf, a, b, x0, x, status, xtol, ftol, max_iter, deriv_floor)
    procedure(root_fn_df_t)            :: fdf
    real(dp),               intent(in) :: a, b, x0
    real(dp),               intent(out):: x
    type(fortnum_status_t), intent(out):: status
    real(dp), intent(in), optional     :: xtol, ftol, deriv_floor
    integer,  intent(in), optional     :: max_iter
```

Newton-Raphson with bisection guard. `fdf` returns both f(x) and f'(x).
Falls back to bisection when a Newton step leaves [a, b] or when
|f'(x)| < `deriv_floor` (default 1e-14).

### `root_brent(f, a, b, x, status [, xtol, ftol, max_iter])`

```fortran
subroutine root_brent(f, a, b, x, status, xtol, ftol, max_iter)
    procedure(root_fn_t)               :: f
    real(dp),               intent(in) :: a, b
    real(dp),               intent(out):: x
    type(fortnum_status_t), intent(out):: status
    real(dp), intent(in), optional     :: xtol, ftol
    integer,  intent(in), optional     :: max_iter
```

Brent's method: inverse quadratic interpolation, secant, and bisection combined.
Superlinear convergence on smooth functions; guaranteed linear in the worst case.
Requires f(a)*f(b) <= 0. Defaults same as `root_bisect`.

---

## fortnum_rng

Counter-based pseudorandom number generator (Threefry-2x64-20) with explicit
caller-owned state. No module-level mutable state; two threads holding separate
`rng_t` objects never race.

Derivative policy: `primal_only` (ad.md §4). A pseudorandom draw is not a
differentiable function of its seed. No derivative entry point exists and none is
reserved. Gradients of estimators built on draws live in the caller, not here.

### Type: `rng_t`

```fortran
type :: rng_t
    integer(int64) :: key(2)       = 0_int64
    integer(int64) :: counter(2)   = 0_int64
    integer(int64) :: buffer       = 0_int64
    logical        :: have_buffer  = .false.
    real(dp)       :: spare_normal = 0.0_dp
    logical        :: have_spare   = .false.
end type rng_t
```

Reproducible state: `key` and `counter`. The buffer and spare fields are derived
caches that `rng_seed` and `rng_split` clear.

### `rng_seed(g, seed, status)`

```fortran
pure subroutine rng_seed(g, seed, status)
    type(rng_t),            intent(out) :: g
    integer(int64),         intent(in)  :: seed
    type(fortnum_status_t), intent(out) :: status
```

Initialize a generator from a 64-bit seed. SplitMix64 runs twice to spread the
seed across the 128-bit key. `status` is always `FORTNUM_OK`.

### `rng_split(parent, stream, child, status)`

```fortran
pure subroutine rng_split(parent, stream, child, status)
    type(rng_t),            intent(in)  :: parent
    integer(int64),         intent(in)  :: stream
    type(rng_t),            intent(out) :: child
    type(fortnum_status_t), intent(out) :: status
```

Derive an independent substream. `child%key` differs from `parent%key` by the
stream index xored into the low word; the counter resets to zero. Does not
advance `parent`. Negative `stream` returns `FORTNUM_DOMAIN_ERROR`.

### `rng_next_u64(g, value)`

```fortran
pure subroutine rng_next_u64(g, value)
    type(rng_t),    intent(inout) :: g
    integer(int64), intent(out)   :: value
```

Next uniform 64-bit word. One Threefry-2x64 call yields two words; the second is
cached for the next call with no cipher call needed.

### `rng_uniform(g, value)`

```fortran
pure subroutine rng_uniform(g, value)
    type(rng_t), intent(inout) :: g
    real(dp),    intent(out)   :: value
```

Uniform real in [0, 1). Keeps the top 53 bits of a u64 and scales by 2^-53.

### `rng_normal(g, value)`

```fortran
pure subroutine rng_normal(g, value)
    type(rng_t), intent(inout) :: g
    real(dp),    intent(out)   :: value
```

Standard normal deviate by Box-Muller. Two uniforms per pair; the second is
cached. `u1` is floored before the log so the [0,1) lower endpoint never
overflows.

### `rng_threefry2x64(key, counter, out)`

```fortran
pure subroutine rng_threefry2x64(key, counter, out)
    integer(int64), intent(in)  :: key(2)
    integer(int64), intent(in)  :: counter(2)
    integer(int64), intent(out) :: out(2)
```

Bare Threefry-2x64-20 block cipher. Exposed for the determinism gate KAT against
published Random123 vectors.

---

## fortnum_interp

Binary search in a sorted grid.

Derivative policy: `primal_only` for `grid_search` (ad.md §4). The returned
cell index is an integer produced by control flow branching on the grid values;
it carries no derivative. Smooth derivatives of the interpolant are valid only
inside a fixed cell. Callers must hold the cell index fixed when differentiating
with respect to the evaluation point.

### `grid_search(p, nmin, nmax, xi, i)`

```fortran
pure subroutine grid_search(p, nmin, nmax, xi, i)
    integer,  intent(in)  :: nmin, nmax
    real(dp), intent(in)  :: p(nmin:nmax)
    real(dp), intent(in)  :: xi
    integer,  intent(out) :: i
```

Binary search in strictly increasing `p(nmin:nmax)`. Returns `i` in
[nmin+1, nmax] such that `p(i-1) < xi <= p(i)`. Clamped at boundaries:
`xi <= p(nmin)` gives `i = nmin+1`; `xi >= p(nmax)` gives `i = nmax`.
O(log2(nmax-nmin)) comparisons.

---

## fortnum_polynomial

Lagrange interpolation weights on an arbitrary node set.

Default derivative policy: `analytic_rule` (ad.md §4). Value weights are
linear in the nodal values; derivative weights are a closed-form analytic rule.

Reserved derivative products for #40: `lagrange_weights_jvp`.

### `lagrange_weights(n, x, xp, coef)`

```fortran
pure subroutine lagrange_weights(n, x, xp, coef)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: x
    real(dp), intent(in)  :: xp(n)
    real(dp), intent(out) :: coef(n)
```

Compute Lagrange value weights at `x` over nodes `xp(1:n)`. The interpolant is
`p(x) = sum_i f(xp(i)) * coef(i)`. Nodes must be distinct.

### `lagrange_deriv_weights(n, x, xp, dcoef)`

```fortran
pure subroutine lagrange_deriv_weights(n, x, xp, dcoef)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: x
    real(dp), intent(in)  :: xp(n)
    real(dp), intent(out) :: dcoef(n)
```

Compute Lagrange derivative weights at `x`. The derivative of the interpolant is
`p'(x) = sum_i f(xp(i)) * dcoef(i)`. Weights depend only on `x` and `xp`, not
on `f`; precompute once and reuse over many function-value vectors.

---

## fortnum_oracle

Testing infrastructure. Reads CSV reference tables produced by
`test/oracle/gen_oracle.py` and checks a Fortran primal against them.

### Type: `oracle_table_t`

```fortran
type :: oracle_table_t
    character(:),  allocatable :: name
    real(dp),      allocatable :: x(:)
    real(dp),      allocatable :: primal(:)
    real(dp),      allocatable :: derivative(:)
    logical                    :: has_derivative = .false.
end type oracle_table_t
```

`has_derivative` mirrors the CSV header flag. When `.false.` the derivative
column is present but carries placeholder values.

### Abstract interface: `oracle_primal_fn`

```fortran
abstract interface
    pure function oracle_primal_fn(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
    end function oracle_primal_fn
end interface
```

### `oracle_read(path, table, status)`

```fortran
subroutine oracle_read(path, table, status)
    character(*),           intent(in)  :: path
    type(oracle_table_t),   intent(out) :: table
    type(fortnum_status_t), intent(out) :: status
```

Read a CSV oracle table from `path`. Header lines start with `#`; the
`has_derivative:` header sets the flag. `FORTNUM_DOMAIN_ERROR` on open failure
or malformed row.

### `oracle_check(table, f, atol, rtol, status)`

```fortran
subroutine oracle_check(table, f, atol, rtol, status)
    type(oracle_table_t),    intent(in)  :: table
    procedure(oracle_primal_fn)           :: f
    real(dp),                intent(in)  :: atol, rtol
    type(fortnum_status_t),  intent(out) :: status
```

Check `f` against every row: `|got - expected| <= atol + rtol*|expected|`. On
the first failure writes a report to `error_unit`; continues to scan all rows.
Sets `FORTNUM_CONVERGENCE_ERROR` if any row fails.

---

## Derivative products

Additive routines from issue #40. Primal signatures are unchanged. Each is
verified against central finite difference and, where both directions exist,
the adjoint identity `uᵀ(Jv) = vᵀ(Jᵀu)`. Active and inactive arguments follow
`docs/design/ad.md §3`: integer orders, node counts, tolerances, step counts,
and seeds are inactive.

### fortnum_special (`analytic_rule`)

```fortran
subroutine dawson_jvp(x, v, jv)        ! F'(x)v, F'(x) = 1 - 2 x F(x)
subroutine dawson_grad(x, u, jtu)      ! scalar adjoint u F'(x)
subroutine gamma_lower_jvp(x, v, jv)   ! d/dx gamma_lower(a,x) = x^(a-1) e^(-x); a inactive
subroutine bessel_in_jvp(n, x, v, jv)  ! dI_n/dx = (I_{n-1}+I_{n+1})/2  (DLMF 10.29.2)
subroutine bessel_kn_jvp(n, x, v, jv)  ! dK_n/dx = -(K_{n-1}+K_{n+1})/2 (DLMF 10.29.4)
```

Derivative of `gamma_lower` with respect to the shape `a` is deferred (digamma
series). HVP is deferred for all three: no scalar loss context in the module.

### fortnum_fft (`transparent`, linear)

```fortran
subroutine fft_c2c_jvp(dz, sign)       ! Jv = forward DFT of the tangent
subroutine fft_c2c_vjp(u, sign)        ! Jᵀu, real-adjoint under <a,b> = Re(sum conjg(a) b)
subroutine fft_r2c_jvp(dx, dc, plan)
subroutine fft_r2c_vjp(u, xbar, plan)
```

The DFT is C-linear: the Hessian is zero, so HVP carries no information and is
omitted. No `grad`: the transforms are vector-valued.

### fortnum_quadrature (`transparent` in integrand values)

```fortran
subroutine gauss_legendre_jvp(w, v, jv)   ! dI = sum_i w_i v_i for tangent integrand values
subroutine gauss_legendre_vjp(w, u, jtu)  ! jtu_i = u w_i
subroutine gauss_legendre_grad(w, grad)   ! dI/df_i = w_i (the weight vector)
```

The map `f -> I` is linear; HVP is zero and omitted. Nodes and order are
inactive.

### fortnum_integrate (`trace_rule`)

```fortran
subroutine integrate_qag_jvp(dfdp, result, di_dp, status, ctx)
```

Differentiates at the frozen accepted subdivision the primal chose: `dI/dp` is
the sum over accepted panels of the Gauss-Kronrod quadrature of `df/dp` at the
frozen nodes. `status` reports non-smoothness when a perturbation would change
the subdivision. The scalar output makes the reverse product equal the forward
sensitivity, so no separate `vjp`/`grad` is shipped.

### fortnum_ode (`trace_rule`)

```fortran
subroutine ode_var_rhs_t(t, y, s, dsdt, ctx)   ! abstract interface: variational/tangent RHS
subroutine ode_integrate_jvp(problem, var_rhs, s0, solution, s1, status)      ! forward sensitivity
subroutine ode_integrate_vjp(problem, var_rhs_adj, u, solution, jtu, status)  ! discrete adjoint
```

Both products re-run the Cash-Karp stepper over the recorded `solution%t` /
`solution%h` schedule held fixed. The forward pass propagates primal and tangent
in lockstep; the adjoint walks the same trace backward. HVP needs a
caller-defined scalar loss on `y(t1)` and is deferred.

### fortnum_roots (`implicit_rule`)

```fortran
subroutine root_jvp(f_x, f_p, tp, dx, status, deriv_floor)   ! dx* = -(f_p . tp)/f_x
subroutine root_vjp(f_x, f_p, u, jtu, status, deriv_floor)   ! jtu_i = -(f_p_i/f_x) u
subroutine root_grad(f_x, f_p, dxdp, status, deriv_floor)    ! dx*/dp = -f_p/f_x
```

The implicit function theorem gives the derivative from the converged root
without differentiating the iteration. `status` reports unreliability when
`|f_x|` is near zero (near-multiple root).

### fortnum_polynomial / fortnum_interp (`transparent`, fixed cell)

```fortran
subroutine lagrange_weights_jvp(n, x, xp, f, vx, jv)   ! d p(x)/dx using derivative weights
subroutine lagrange_weights_vjp(n, x, xp, f, u, jtu)
subroutine lagrange_fval_jvp(n, x, xp, vf, jv)         ! d p(x)/df_i using value weights
subroutine lagrange_fval_vjp(n, x, xp, u, jtu)
```

Derivatives are valid inside a fixed cell. Crossing a cell boundary is a
non-smooth event; the caller holds the index fixed when differentiating with
respect to the evaluation point.

### fortnum_rng (`primal_only`)

A pseudorandom draw is not a differentiable function of its seed. No derivative
products.
