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

### fortnum_special_complex_bessel

Complex-argument Bessel functions of integer order: J_n(z), I_n(z), K_n(z).
Clean-room from DLMF chapter 10. The `scaled` flag on I and K factors out
e^{Re z} (KODE scaling), needed where Re z is large. Negative order is handled
internally: J_{-n} = (-1)^n J_n, I_{-n} = I_n, K_{-n} = K_n.

Policy: `analytic_rule`. Active argument: `z` (complex). Inactive: `order`,
`scaled`. Derivatives use the recurrences J_n'(z) = (J_{n-1}-J_{n+1})/2 (DLMF
10.6.1), I_n'(z) = (I_{n-1}+I_{n+1})/2 (DLMF 10.29.2), K_n'(z) =
-(K_{n-1}+K_{n+1})/2 (DLMF 10.29.4); forward complex products below.

#### `bessel_j_complex(order, z, result, status)`

```fortran
subroutine bessel_j_complex(order, z, result, status)
    integer,                intent(in)  :: order
    complex(dp),            intent(in)  :: z
    complex(dp),            intent(out) :: result
    type(fortnum_status_t), intent(out) :: status
```

J_order(z) by power series near the origin, Hankel asymptotic for large |z|,
and downward (Miller) recurrence in the near-real strip.

#### `bessel_i_complex(order, z, scaled, result, status)`

```fortran
subroutine bessel_i_complex(order, z, scaled, result, status)
    integer,                intent(in)  :: order
    complex(dp),            intent(in)  :: z
    logical,                intent(in)  :: scaled
    complex(dp),            intent(out) :: result
    type(fortnum_status_t), intent(out) :: status
```

I_order(z). With `scaled = .true.` returns e^{-Re z} I_order(z).

#### `bessel_i_complex_array(order0, nseq, z, scaled, result, status)`

```fortran
subroutine bessel_i_complex_array(order0, nseq, z, scaled, result, status)
    integer,                intent(in)  :: order0
    integer,                intent(in)  :: nseq
    complex(dp),            intent(in)  :: z
    logical,                intent(in)  :: scaled
    complex(dp),            intent(out) :: result(nseq)
    type(fortnum_status_t), intent(out) :: status
```

I_{order0}(z) through I_{order0+nseq-1}(z) in one downward-recurrence pass.

#### `bessel_k_complex(order, z, scaled, result, status)`

```fortran
subroutine bessel_k_complex(order, z, scaled, result, status)
    integer,                intent(in)  :: order
    complex(dp),            intent(in)  :: z
    logical,                intent(in)  :: scaled
    complex(dp),            intent(out) :: result
    type(fortnum_status_t), intent(out) :: status
```

K_order(z) by the integral representation 10.32.18 on the right half-plane
(`Re z > 0`). With `scaled = .true.` returns e^{z} K_order(z).

#### `bessel_k_complex_array(order0, nseq, z, scaled, result, status)`

```fortran
subroutine bessel_k_complex_array(order0, nseq, z, scaled, result, status)
    integer,                intent(in)  :: order0
    integer,                intent(in)  :: nseq
    complex(dp),            intent(in)  :: z
    logical,                intent(in)  :: scaled
    complex(dp),            intent(out) :: result(nseq)
    type(fortnum_status_t), intent(out) :: status
```

K_{order0}(z) through K_{order0+nseq-1}(z) by upward recurrence in order.

### fortnum_special_hypergeometric_1f1

Confluent hypergeometric function 1F1(a;b;z) (Kummer M), complex a, b, z.
Clean-room from DLMF chapter 13. Holomorphic in a and z; poles in b at the
non-positive integers. One selector switches between Kummer transformation
(Re z < 0), Taylor series (small |z|), and large-z asymptotics; only the
converged value is exposed.

Policy: `analytic_rule`. Active argument: `z`. Inactive: `a`, `b`.
d/dz M = (a/b) M(a+1,b+1,z) (DLMF 13.3.15). `hyperg_1f1_a1` fixes a = 1, the
form the FLR / plasma-dispersion consumers need; its forward and reverse
products treat z as a real 2-vector (Re z, Im z).

#### `hyperg_1f1(a, b, z, result, status)`

```fortran
subroutine hyperg_1f1(a, b, z, result, status)
    complex(dp),            intent(in)  :: a, b, z
    complex(dp),            intent(out) :: result
    type(fortnum_status_t), intent(out) :: status
```

M(a, b, z). `FORTNUM_DOMAIN_ERROR` when b is at or near a non-positive integer.

#### `hyperg_1f1_a1(b, z, result, status)`

```fortran
subroutine hyperg_1f1_a1(b, z, result, status)
    complex(dp),            intent(in)  :: b, z
    complex(dp),            intent(out) :: result
    type(fortnum_status_t), intent(out) :: status
```

M(1, b, z); thin specialization forwarding to `hyperg_1f1` with a = 1.

### fortnum_special_erf_cbind

erf/erfc provider for the C ABI. Fortran callers use the F2008 intrinsics
directly; this module forwards to them with no reimplementation so the ABI
layer has C-callable symbols.

Policy: `transparent`. Active argument: `x` (real scalar). d/dx erf(x) =
2/sqrt(pi) e^{-x^2}, d/dx erfc(x) = -2/sqrt(pi) e^{-x^2} (DLMF 7.2.1).

#### `fortnum_erf(x) -> real(dp)` (elemental)

```fortran
elemental function fortnum_erf(x) result(y)
    real(dp), intent(in) :: x
```

Error function erf(x). Forwards to the intrinsic `erf`.

#### `fortnum_erfc(x) -> real(dp)` (elemental)

```fortran
elemental function fortnum_erfc(x) result(y)
    real(dp), intent(in) :: x
```

Complementary error function erfc(x) = 1 - erf(x). Forwards to the intrinsic
`erfc`.

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

### `gauss_gen_laguerre(n, alpha, x, w)`

```fortran
pure subroutine gauss_gen_laguerre(n, alpha, x, w)
    integer,  intent(in)  :: n
    real(dp), intent(in)  :: alpha
    real(dp), intent(out) :: x(n), w(n)
```

Generalized Gauss-Laguerre nodes and weights for the weight
`w(x) = x^alpha exp(-x)` on `[0, inf)`, `alpha > -1`. The quadrature sum
`sum_i w(i)*f(x(i))` approximates the integral of `x^alpha exp(-x) f(x)` over
`[0, inf)`. Nodes ascending. Golub-Welsch on the symmetric tridiagonal Jacobi
matrix, solved by clean-room implicit-shift QL (no LAPACK), with zeroth moment
`Gamma(alpha+1)`. Requires `n >= 1` and `alpha > -1`. Same `transparent`
derivative policy as `gauss_legendre`: the linear map `f -> I` reuses
`gauss_legendre_jvp` / `gauss_legendre_vjp` / `gauss_legendre_grad` on the
returned weight vector.

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

## fortnum_levin

Levin-u sequence acceleration. Given a finite recorded sequence of series terms
a_0..a_{n-1}, it forms the partial sums and returns the Levin-u accelerated
limit with an error estimate, summing slowly convergent or divergent series
such as the Kummer-ratio sums in the confluent hypergeometric consumers.

Algorithm: Levin (1973), Weniger (1989) eqs 7.2-8 / 7.3-9, Fessler-Ford-Smith
(1983). Remainder estimate omega_n = (1+n) a_n (the u variant). The table grows
one term at a time; the diagonal value that moves least, scaled by a rounding
floor, is the best estimate.

Policy: `primal_only`. The accelerated value is a nonlinear rational transform
of the terms and the selected order is data-dependent, so no derivative product
applies. Active: `terms`. Inactive: `n`.

### `levin_u_accel(terms, n, sum_accel, abserr, status)`

```fortran
subroutine levin_u_accel(terms, n, sum_accel, abserr, status)
    integer,                intent(in)  :: n
    real(dp),               intent(in)  :: terms(n)
    real(dp),               intent(out) :: sum_accel
    real(dp),               intent(out) :: abserr
    type(fortnum_status_t), intent(out) :: status
```

Accelerate the series with terms `terms(1:n)` (Levin-u sequence acceleration).

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

### fortnum_ode_dop853

Adaptive Prince-Dormand RK8(7)-13M integrator. Twelve explicit stages give an
eighth-order solution; embedded order-5 and order-3 estimators drive the PI
step control. Same carriers (`ode_problem_t`, `ode_workspace_t`,
`ode_solution_t`) and the same `trace_rule` policy as the Cash-Karp path; use
it for the very tight tolerances where the higher order pays off. Active: `y0`,
`ctx` parameters. Inactive: `rtol`, `atol`, step bounds.

#### `ode_integrate_dop(problem, workspace, solution, status)`

```fortran
subroutine ode_integrate_dop(problem, workspace, solution, status)
    type(ode_problem_t),    intent(in)    :: problem
    type(ode_workspace_t),  intent(inout) :: workspace
    type(ode_solution_t),   intent(inout) :: solution
    type(fortnum_status_t), intent(out)   :: status
```

Integrate `problem` over its span, recording the accepted-step mesh in
`solution`. Drop-in replacement for `ode_integrate` using the DOP853 stepper.

#### `ode_solve_dop(rhs, t0, t1, y0, t_out, y_out, status [, rtol, atol])`

```fortran
subroutine ode_solve_dop(rhs, t0, t1, y0, t_out, y_out, status, rtol, atol)
    procedure(ode_rhs_t)               :: rhs
    real(dp),               intent(in) :: t0, t1
    real(dp),               intent(in) :: y0(:)
    real(dp), allocatable, intent(out) :: t_out(:)
    real(dp), allocatable, intent(out) :: y_out(:,:)
    type(fortnum_status_t), intent(out) :: status
    real(dp), intent(in), optional     :: rtol, atol
```

Flat DOP853 call. Allocates and returns the accepted mesh `t_out` and states
`y_out(neq, size(t_out))`. Mirrors `ode_solve`.

#### `dop853_step(rhs, t, y, h, ...)`

```fortran
subroutine dop853_step(rhs, t, y, h, have_k1, k1, k2, k3, k4, k5, k6, &
    k7, k8, k9, k10, k11, k12, ytmp, y8, err5, err3, nfev)
```

One RK8(7)-13M stage block: from `(t, y)` it fills the order-8 update `y8` and
the order-5 and order-3 error vectors `err5`, `err3` for the step controller.
The stage slots and `nfev` are caller-owned workspace; `have_k1` reuses the
first stage (FSAL) across accepted steps. Used by `ode_integrate_dop`; exposed
for callers driving a custom mesh.

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

### fortnum_multiroot

Multidimensional root finding F(x) = 0 in R^n, plus two scalar utilities the
same consumer needs. Both solvers run a Newton step with a backtracking line
search on 1/2 |F|^2; the Newton system is solved by Gaussian elimination with
partial pivoting. A singular Jacobian reports `FORTNUM_DOMAIN_ERROR`.

Solver policy: `implicit_rule`; the root satisfies F(x*, p) = 0 and the implicit
function theorem gives dx*/dp = -J_x^{-1} J_p. `deriv_central` and `argsort` are
`primal_only`. Reserved n-dim products: `multiroot_jvp`, `multiroot_vjp`,
`multiroot_grad`.

Callback abstract interfaces: `multiroot_fdf_t` returns F and the analytic
Jacobian; `multiroot_fn_t` returns F only (Jacobian by finite difference);
`deriv_fn_t` is the scalar function differenced by `deriv_central`. Each takes
the optional unlimited-polymorphic `ctx`.

```fortran
abstract interface
    function deriv_fn_t(x, ctx) result(fx)
        real(dp),           intent(in) :: x
        class(*), optional, intent(in) :: ctx
        real(dp)                       :: fx
    end function deriv_fn_t
end interface
```

#### `multiroot_hybrid(fdf, n, x0, x, status [, xtol, ftol, max_iter, ctx])`

```fortran
subroutine multiroot_hybrid(fdf, n, x0, x, status, xtol, ftol, max_iter, ctx)
    procedure(multiroot_fdf_t)         :: fdf
    integer,                intent(in) :: n
    real(dp),               intent(in) :: x0(n)
    real(dp),               intent(out):: x(n)
    type(fortnum_status_t), intent(out):: status
    real(dp), intent(in), optional     :: xtol, ftol
    integer,  intent(in), optional     :: max_iter
    class(*), intent(in), optional     :: ctx
```

Hybrid solve with an analytic Jacobian supplied by the `fdf` callback
(`multiroot_fdf_t`); Powell hybrid dogleg with analytic Jacobian.

#### `multiroot_hybrids(fn, n, x0, x, status [, xtol, ftol, max_iter, ctx])`

```fortran
subroutine multiroot_hybrids(fn, n, x0, x, status, xtol, ftol, max_iter, ctx)
    procedure(multiroot_fn_t)          :: fn
    integer,                intent(in) :: n
    real(dp),               intent(in) :: x0(n)
    real(dp),               intent(out):: x(n)
    type(fortnum_status_t), intent(out):: status
    real(dp), intent(in), optional     :: xtol, ftol
    integer,  intent(in), optional     :: max_iter
    class(*), intent(in), optional     :: ctx
```

Same iteration with the Jacobian built column by column by central differences
of the residual `fn` (`multiroot_fn_t`); finite-difference-Jacobian variant of
the Powell hybrid dogleg.

#### `deriv_central(f, x, h, result, abserr, status [, ctx])`

```fortran
subroutine deriv_central(f, x, h, result, abserr, status, ctx)
    procedure(deriv_fn_t)               :: f
    real(dp),               intent(in)  :: x, h
    real(dp),               intent(out) :: result, abserr
    type(fortnum_status_t), intent(out) :: status
    class(*), intent(in), optional      :: ctx
```

Central finite-difference first derivative of `f` at `x` with step `h`, with a
Richardson truncation/round-off error estimate in `abserr`. The operation is
finite differencing, so it is `primal_only`: use the value, do not AD through it.

#### `argsort(x, perm)`

```fortran
pure subroutine argsort(x, perm)
    real(dp), intent(in)  :: x(:)
    integer,  intent(out) :: perm(size(x))
```

Ascending stable index permutation: `x(perm)` is sorted. `primal_only`.

---

## fortnum_roots_complex

Distinct zeros of an analytic function inside an axis-aligned complex rectangle,
with multiplicities. Replaces the ZEAL ICON=3 region search that KIM's
`wkb_dispersion` used to locate dispersion-relation roots.

The count of zeros in the box (with multiplicity) is the winding number
N = (1/2 pi i) oint_C f'/f dz. The power-sum moments
s_p = (1/2 pi i) oint_C z^p f'/f dz are the Newton sums of the zeros; the
distinct zeros are the generalized eigenvalues of the Hankel pencil built from
the moments (Kravanja and Van Barel 2000), recovered with LAPACK ZGGEV and
polished by complex Newton. Multiplicities solve the confluent Vandermonde
system. When the box holds more than `m_max` zeros it is bisected and each half
recursed, as ICON=3 does. The contour integrals run on the four edges with the
package Gauss-Kronrod driver; f'/f uses a complex central difference, so no
branch cut of f is crossed.

Derivative policy: `implicit_rule` with `differentiate_through=false` (ad.md §4).
A zero satisfies f(z*, p) = 0, so dz*/dp = -f_p/f'(z*) by the implicit function
theorem, without differentiating the contour integral or the eigensolve. No
consumer differentiates through the region search, so this module provides no
JVP/VJP; the scalar implicit rule in `fortnum_roots` covers per-root
sensitivity. Inactive: `m_max`, the tolerances, `max_refine`, status,
multiplicities.

### Abstract interface: `complex_root_fn_t`

```fortran
abstract interface
    subroutine complex_root_fn_t(kr, fk, ctx)
        complex(dp), intent(in)  :: kr
        complex(dp), intent(out) :: fk
        class(*),    intent(in), optional :: ctx
    end subroutine complex_root_fn_t
end interface
```

Analytic f whose zeros are sought. `ctx` forwards parameters; the finder never
inspects it.

### `complex_region_roots(f, ll, ur, roots, fvals, mult, nfound, status [, m_max, newtonz, newtonf, max_refine, ctx])`

```fortran
subroutine complex_region_roots(f, ll, ur, roots, fvals, mult, nfound, &
        status, m_max, newtonz, newtonf, max_refine, ctx)
    procedure(complex_root_fn_t)            :: f
    complex(dp),             intent(in)     :: ll, ur
    complex(dp), allocatable, intent(out)   :: roots(:), fvals(:)
    integer,     allocatable, intent(out)   :: mult(:)
    integer,                 intent(out)    :: nfound
    type(fortnum_status_t),  intent(out)    :: status
    integer,     intent(in), optional       :: m_max, max_refine
    real(dp),    intent(in), optional       :: newtonz, newtonf
    class(*),    intent(in), optional       :: ctx
```

Find the distinct zeros of `f` in the rectangle with corners `ll` (lower-left)
and `ur` (upper-right). On return `roots(1:nfound)` are the distinct zeros,
`fvals` is `f` at each, and `mult` their multiplicities; the multiplicities sum
to the winding number. `m_max` (default 5) caps the zeros handled per subregion
before bisection; `newtonz` (default 5e-8) and `newtonf` (default 1e-14) are the
Newton stopping tolerances on `|dz|` and `|f|`; `max_refine` (default 60) bounds
the Newton iterations. A degenerate box or `m_max < 1` reports
`FORTNUM_DOMAIN_ERROR`; a zero or pole on the contour, a failed edge integral, a
failed eigensolve, or unconverged Newton reports `FORTNUM_CONVERGENCE_ERROR`.

This solver links LAPACK (ZGGEV) and its BLAS dependency. It is the only
fortnum module with an external numerical dependency.

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

## fortnum_bspline

B-spline basis functions and derivatives over a breakpoint set, in the clamped-knot
convention: order `k = degree + 1`, `nbreak` breakpoints, end breakpoints
repeated with multiplicity `k` (clamped basis). The coefficient count is
`ncoef = nbreak + k - 2`. Cox-de Boor recursion (de Boor 2001; Piegl & Tiller).

Policy: `transparent` inside a fixed knot span; within one span the spline value
sum_i c_i B_{i,k}(x) is a fixed polynomial in x. The span index from
`bspline_span_index` is `primal_only`: crossing a knot is non-smooth, so hold
the span fixed when differentiating in x. Active: `x`. Inactive: order `k`, the
breakpoint array, `nderiv`.

### Type: `bspline_workspace_t`

Caller-owned workspace holding the order, breakpoints, augmented knot vector,
and `ncoef`. Built by `bspline_init` and `bspline_set_knots`.

### `bspline_init(ws, order, nbreak, status)`

```fortran
subroutine bspline_init(ws, order, nbreak, status)
    type(bspline_workspace_t), intent(out) :: ws
    integer,                   intent(in)  :: order
    integer,                   intent(in)  :: nbreak
    type(fortnum_status_t),    intent(out) :: status
```

Allocate the workspace for spline order `order` over `nbreak` breakpoints.

### `bspline_set_knots(ws, breakpts, status)`

```fortran
subroutine bspline_set_knots(ws, breakpts, status)
    type(bspline_workspace_t), intent(inout) :: ws
    real(dp),                  intent(in)     :: breakpts(:)
    type(fortnum_status_t),    intent(out)    :: status
```

Build the clamped augmented knot vector from strictly increasing `breakpts`.

### `bspline_eval_basis(ws, x, values, status)`

```fortran
subroutine bspline_eval_basis(ws, x, values, status)
    type(bspline_workspace_t), intent(in)  :: ws
    real(dp),                  intent(in)  :: x
    real(dp),                  intent(out) :: values(:)
    type(fortnum_status_t),    intent(out) :: status
```

Scatter the `ncoef` basis-function values at `x` into `values(1:ncoef)`; the
nonzero entries are the `k` functions on x's span.

### `bspline_eval_deriv(ws, x, nderiv, dvalues, status)`

```fortran
subroutine bspline_eval_deriv(ws, x, nderiv, dvalues, status)
    type(bspline_workspace_t), intent(in)  :: ws
    real(dp),                  intent(in)  :: x
    integer,                   intent(in)  :: nderiv
    real(dp),                  intent(out) :: dvalues(0:, :)
    type(fortnum_status_t),    intent(out) :: status
```

Fill `dvalues(0:nderiv, 1:ncoef)` with the basis functions and their
derivatives up to order `nderiv` at `x`.

### `bspline_span_index(ws, x) -> integer` (pure function)

```fortran
pure function bspline_span_index(ws, x) result(span)
    type(bspline_workspace_t), intent(in) :: ws
    real(dp),                  intent(in) :: x
```

Knot span containing `x`, by binary search. `primal_only`.

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
subroutine bessel_j_complex_jvp(order, z, v, jv, status)  ! J_n'(z)v (DLMF 10.6.1)
subroutine bessel_i_complex_jvp(order, z, scaled, v, jv, status)  ! I_n'(z)v (DLMF 10.29.2)
subroutine bessel_k_complex_jvp(order, z, scaled, v, jv, status)  ! K_n'(z)v (DLMF 10.29.4)
subroutine hyperg_1f1_a1_jvp(z, b, v, jv)   ! d/dz M(1,b,z) v; z as (Re,Im) 2-vector
subroutine hyperg_1f1_a1_vjp(z, b, u, jtu)  ! reverse product of the same map
subroutine fortnum_erf_jvp(x, v, jv)        ! 2/sqrt(pi) e^{-x^2} v (DLMF 7.2.1)
subroutine fortnum_erfc_jvp(x, v, jv)       ! -2/sqrt(pi) e^{-x^2} v
subroutine fortnum_erf_grad(x, u, jtu)      ! scalar adjoint (VJP)
subroutine fortnum_erfc_grad(x, u, jtu)
```

Derivative of `gamma_lower` with respect to the shape `a` is deferred (digamma
series). HVP is deferred for all: no scalar loss context in the module. The
complex Bessel products carry `status` because the primal evaluates a status
path. `hyperg_1f1_a1` differentiates the analytic map z -> M, exercised as the
Cauchy-Riemann Jacobian on (Re z, Im z); `erf`/`erfc` are `transparent`.

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
inactive. `gauss_gen_laguerre` shares this policy: its weight vector drops
into the same three products, since the rule (Legendre or generalized Laguerre)
only sets the inactive weights.

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

### fortnum_bspline (`transparent`, fixed span)

```fortran
subroutine bspline_eval_jvp(ws, x, coef, vx, jv, status)   ! d/dx [sum c_i B_i(x)] vx
subroutine bspline_eval_vjp(ws, x, coef, u, jtu, status)   ! (d/dx [sum c_i B_i(x)])^T u
```

The spline value is differentiated in `x` at the fixed span; the coefficients
`coef` are inactive. Crossing a knot is non-smooth, so the caller holds the span
fixed (the AD test checks this). Both directions match the analytic recurrence
weights from `bspline_eval_deriv`.

### fortnum_rng (`primal_only`)

A pseudorandom draw is not a differentiable function of its seed. No derivative
products.
