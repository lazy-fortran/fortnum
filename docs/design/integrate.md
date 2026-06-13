# ADR: fortnum_integrate API

Status: accepted (issue #23, M4.1). Normative for the implementation issues
#24 (QAG/QAGS primal) and #25 (QAGP/QAGIU primal), and for #40 (derivative
products).

This ADR is subordinate to `docs/design/ad.md`. The derivative vocabulary,
policy classes, naming convention, and active-argument rules come from that
document; this one cites it by section and fixes the module-specific Fortran.
The default policy class for `integrate` is `trace_rule` (ad.md §4): the
adaptive integrator picks its subdivisions from the integrand, so a derivative
differentiates the accepted subdivision with that subdivision held fixed.

The module builds on the single-panel Gauss-Kronrod kernel that already ships
in `src/quadrature/fortnum_integrate_gk.f90`. That kernel owns the GK node and
weight tables and the per-panel error estimate; this module owns the global
adaptive driver, the bisection bookkeeping, the epsilon extrapolation, the
interval transforms, and the recorded trace. An implementer reuses `gk_apply`
for every panel evaluation and does not re-derive the rule. Where this file
underspecifies a constant, the QUADPACK routines it names (`dqag`, `dqags`,
`dqagp`, `dqagi`) fix it. Piessens, de Doncker-Kapenga, Ueberhuber, Kahaner,
"QUADPACK" (Springer 1983).

## 1. Scope and method

Four globally adaptive integrators over a finite or transformed interval:

- `QAG`: globally adaptive, selectable GK rule, for integrands with no known
  singularity.
- `QAGS`: adaptive bisection plus Wynn epsilon extrapolation, for endpoint
  singularities and integrable algebraic or logarithmic spikes.
- `QAGP`: QAGS with user-supplied break points that seed the initial
  subdivision at known interior singularities or kinks.
- `QAGIU`: semi-infinite and doubly infinite intervals, mapped to `(0,1]` by a
  documented variable transform, then driven by the QAGS machinery.

All four share one driver. The driver bisects the subinterval carrying the
largest local error estimate, re-evaluates the GK pair on the two halves, and
repeats until the summed error estimate meets the tolerance or a stop condition
fires. QAGS, QAGP, and QAGIU add the epsilon table on top of the same bisection
loop; QAG runs the loop alone. This is the QUADPACK structure: one `qag`-style
driver, with `qags` differing only by the extrapolation step.

INACTIVE controls for M4.1: only QAG with the default rule and QAGS are wired
in #24. QAGP break points and QAGIU transforms ship in #25. Singular-mode
selection (the QAGS epsilon path versus the plain QAG path) is selected by the
entry point the caller invokes, not by a runtime flag. The rule `key` is the
only runtime selector and defaults to 21, matching the kernel.

## 2. Integrand interface

The integrand is an abstract interface. The caller passes the abscissa in and
gets the value back, and may thread its own data through an optional unlimited
polymorphic context. No module-level state and no global procedure pointer
carry the integrand; this is what lets a derivative product call the integrand
as many times as it needs without a hidden side channel (ad.md §6).

```fortran
abstract interface
    function integrate_integrand_t(x, ctx) result(fx)
        import :: dp
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
    end function integrate_integrand_t
end interface
```

The optional `ctx` is the only channel for parameters: a caller with
parameters defines a derived type, passes it as `ctx`, and selects its type
inside the integrand. `ctx` is `intent(in)`; the integrand reads parameters, it
does not mutate caller state. A parameter-free integrand ignores the argument.

The kernel's own `gk_integrand_t` takes no `ctx`. The driver bridges the two
with a thin internal closure: it holds the caller's integrand and `ctx` in a
local derived type for the duration of one `integrate_qag` call and presents a
`ctx`-free wrapper to `gk_apply`. The bridge lives on the call stack, so the
no-global-state rule holds. The kernel signature does not change.

## 3. Derived types

Three caller-owned types carry the work stack, the extrapolation table, and the
recorded result. None of them is module-level state; the caller allocates them,
passes them in, and may keep them between calls to avoid reallocation.

### 3.1 Work stack

The subinterval list is the QUADPACK `alist/blist/rlist/elist` workspace held as
one type. Each entry is one subinterval and its GK result.

```fortran
type :: integrate_workspace_t
    integer               :: limit = 0        ! capacity, set on first use
    integer               :: last  = 0         ! subintervals in use
    real(dp), allocatable :: a(:)              ! left endpoints,  a(limit)
    real(dp), allocatable :: b(:)              ! right endpoints, b(limit)
    real(dp), allocatable :: r(:)              ! GK result per subinterval
    real(dp), allocatable :: e(:)              ! GK abserr per subinterval
    integer,  allocatable :: iord(:)           ! error-ordered indices
    integer,  allocatable :: level(:)          ! bisection depth per entry
    integer,  allocatable :: ndin(:)           ! QAGP: 1 if break-seeded
end type integrate_workspace_t
```

`a(i)..b(i)` is the i-th accepted subinterval, `r(i)` its GK result, `e(i)` its
GK error estimate; `iord` keeps the subintervals ordered by descending `e` so
the driver pops the worst one in O(1) and reinserts the two halves. `level`
caps the bisection depth for the QAGS roundoff test. `ndin` is used only by
QAGP to mark break-seeded intervals; QAG and QAGS leave it zero. The driver
allocates all arrays to `limit` once, when `workspace%limit` does not match the
requested limit, and reuses them otherwise. No allocation happens inside the
bisection loop.

### 3.2 Epsilon table

Wynn's epsilon algorithm extrapolates the sequence of partial integral sums to
accelerate convergence at an endpoint singularity. The table is the QUADPACK
`qelg` workspace: a fixed-width column of the last few sums plus the running
extrapolated value and its error estimate.

```fortran
type :: integrate_epstab_t
    integer  :: n = 0                  ! entries currently in the table
    real(dp) :: tab(52) = 0.0_dp       ! epsilon column (QUADPACK qelg width)
    real(dp) :: result = 0.0_dp        ! best extrapolated value so far
    real(dp) :: abserr = huge(0.0_dp)  ! its error estimate
    real(dp) :: res3la(3) = 0.0_dp     ! last three extrapolations, for the
                                       ! roundoff guard
    integer  :: nres = 0               ! count of calls to the extrapolator
end type integrate_epstab_t
```

The width 52 matches QUADPACK `dqelg`. QAG never touches this type; QAGS, QAGP,
and QAGIU pass it to the extrapolator after every bisection once at least three
subintervals exist. The caller owns it; a fresh default-initialized
`integrate_epstab_t` is a valid empty table.

### 3.3 Result and trace

The result type holds the value, the error estimate, the status, and the frozen
accepted subdivision. The subdivision is the ordered list of subintervals the
primal kept, each with its GK result. This is the `trace_rule` schedule
(ad.md §1, §4): once the primal has chosen and frozen it, the derivative of the
integral with respect to the integrand values is a fixed weighted sum over the
panels on that subdivision.

```fortran
type :: integrate_result_t
    real(dp)               :: value  = 0.0_dp   ! the integral estimate
    real(dp)               :: abserr = 0.0_dp   ! final error estimate
    integer                :: neval  = 0        ! integrand evaluations
    integer                :: nsub   = 0        ! accepted subintervals
    real(dp), allocatable  :: sub_a(:)          ! accepted left endpoints
    real(dp), allocatable  :: sub_b(:)          ! accepted right endpoints
    real(dp), allocatable  :: sub_r(:)          ! GK result on each subinterval
    real(dp), allocatable  :: sub_e(:)          ! GK error on each subinterval
    integer                :: key   = 21        ! GK rule used on every panel
    logical                :: extrapolated = .false.  ! epsilon path was used
    type(fortnum_status_t) :: status
end type integrate_result_t
```

`sub_a`, `sub_b`, `sub_r`, `sub_e` are the trimmed, final-length copies of the
work-stack columns in interval order (left to right), not error order. `value`
is their `sub_r` sum for QAG, or the extrapolated value for QAGS when
extrapolation improved on the direct sum. `key` records which GK pair ran so a
derivative product reuses the same panel weights. `extrapolated` records whether
the reported value came from the epsilon table, because a derivative over an
extrapolated value differs from one over the plain panel sum (section 6).

`integrate_result_t` is `intent(inout)` on the primal calls so the caller may
reuse its allocations; the driver reallocates the trace arrays only when `nsub`
exceeds their current length.

## 4. Primal routines (#24, #25)

All four take an integrand, the interval, the convergence controls, a
workspace, an epsilon table, and the result. The workspace and table are always
present in the signature so the four routines stay uniform; QAG ignores the
table. `status` mirrors `result%status` so a caller that keeps only the result
still sees the outcome.

```fortran
subroutine integrate_qag(f, a, b, epsabs, epsrel, workspace, result, &
                         status, key, limit, ctx)
    procedure(integrate_integrand_t)        :: f
    real(dp),                  intent(in)    :: a, b, epsabs, epsrel
    type(integrate_workspace_t), intent(inout) :: workspace
    type(integrate_result_t),  intent(inout) :: result
    type(fortnum_status_t),    intent(out)   :: status
    integer,  intent(in), optional :: key, limit
    class(*), intent(in), optional :: ctx
end subroutine integrate_qag

subroutine integrate_qags(f, a, b, epsabs, epsrel, workspace, epstab, &
                          result, status, limit, ctx)
    procedure(integrate_integrand_t)        :: f
    real(dp),                  intent(in)    :: a, b, epsabs, epsrel
    type(integrate_workspace_t), intent(inout) :: workspace
    type(integrate_epstab_t),  intent(inout) :: epstab
    type(integrate_result_t),  intent(inout) :: result
    type(fortnum_status_t),    intent(out)   :: status
    integer,  intent(in), optional :: limit
    class(*), intent(in), optional :: ctx
end subroutine integrate_qags

subroutine integrate_qagp(f, a, b, points, epsabs, epsrel, workspace, &
                          epstab, result, status, limit, ctx)
    procedure(integrate_integrand_t)        :: f
    real(dp),                  intent(in)    :: a, b, epsabs, epsrel
    real(dp),                  intent(in)    :: points(:)
    type(integrate_workspace_t), intent(inout) :: workspace
    type(integrate_epstab_t),  intent(inout) :: epstab
    type(integrate_result_t),  intent(inout) :: result
    type(fortnum_status_t),    intent(out)   :: status
    integer,  intent(in), optional :: limit
    class(*), intent(in), optional :: ctx
end subroutine integrate_qagp

subroutine integrate_qagiu(f, bound, inf, epsabs, epsrel, workspace, &
                           epstab, result, status, limit, ctx)
    procedure(integrate_integrand_t)        :: f
    real(dp),                  intent(in)    :: bound, epsabs, epsrel
    integer,                   intent(in)    :: inf
    type(integrate_workspace_t), intent(inout) :: workspace
    type(integrate_epstab_t),  intent(inout) :: epstab
    type(integrate_result_t),  intent(inout) :: result
    type(fortnum_status_t),    intent(out)   :: status
    integer,  intent(in), optional :: limit
    class(*), intent(in), optional :: ctx
end subroutine integrate_qagiu
```

`key` defaults to 21 and is forwarded to `gk_apply` unchanged; only 15, 21, 31,
61 are valid. `limit` defaults to 500 and bounds the work-stack capacity.
QAGS, QAGP, and QAGIU run the GK21 pair on every panel as QUADPACK does, so they
take no `key`.

`integrate_qagp` reads the interior break points from `points`; entries outside
`(a, b)` are dropped, the rest seed the initial subdivision so the driver starts
with a panel boundary at each known singularity. An empty `points` reduces QAGP
to QAGS.

`integrate_qagiu` maps the requested interval to `t in (0, 1]` and integrates
the transformed integrand with the QAGS machinery. `inf` selects the interval
and the transform:

- `inf = +1`: `[bound, +inf)`, `x = bound + (1 - t)/t`,
  `dx = dt / t**2`.
- `inf = -1`: `(-inf, bound]`, `x = bound - (1 - t)/t`,
  `dx = dt / t**2`.
- `inf = +2`: `(-inf, +inf)`, split at `bound` and apply the `+1` and `-1`
  transforms to the two halves, then sum.

The Jacobian `1/t**2` multiplies the integrand inside the transformed call; the
endpoint singularity it introduces at `t = 0` is exactly what the QAGS epsilon
path is built to handle. This is the QUADPACK `dqagi` transform.

### 4.1 Convergence control

`epsabs` and `epsrel` are the absolute and relative tolerances; the driver stops
when the summed error estimate is at most `max(epsabs, epsrel*|value|)`.
`limit` caps the number of subintervals. All three are inactive controls
(ad.md §3): they select when to stop, they do not carry a derivative.

Invalid input is rejected before any panel runs: `b <= a` for the finite
routines, non-positive `epsabs` with `epsrel` below the kernel's floor,
`limit < 1`, an unsupported `key`, or an `inf` outside `{-1, +1, +2}`. Each
reports `FORTNUM_DOMAIN_ERROR` with a message naming the offending argument.

### 4.2 Status semantics

The status path distinguishes ordinary convergence failure from a non-smooth
derivative case. This separation is mandatory: a derivative product must be able
to tell a value that converged on a stable subdivision from one whose
subdivision is fragile.

Ordinary convergence failures map to `FORTNUM_CONVERGENCE_ERROR`, with the
message naming the cause:

- max subdivisions: `limit` reached before the tolerance,
- roundoff: the error estimate stopped improving while still above tolerance
  (the QUADPACK `iroff1/iroff2/iroff3` detectors and the extrapolation-table
  roundoff guard),
- slow convergence: the epsilon table diverged or stalled, so QAGS fell back to
  the direct panel sum.

In every such case `result%value` and the trace are still the best primal
estimate and remain usable; only the tolerance was not met.

Non-smooth derivative cases are distinct and map to `FORTNUM_DOMAIN_ERROR`.
These are perturbations of the integrand or the interval that change the
accepted subdivision or the break structure, so the frozen-subdivision
derivative of section 6 does not exist:

- a user break point in `points` that coincides with `a` or `b` to within the
  panel tolerance, so the break structure is ambiguous,
- a QAGIU transform whose finite image collapses (`bound` non-finite, or the
  `inf = +2` split point degenerate),
- an integrand the driver detects as non-integrable on the transformed interval
  (the QUADPACK "bad integrand behaviour" code 3 condition).

The primal still returns its best value when one exists, but the status tells a
later derivative product that differentiating at this subdivision is not
well-posed. A derivative routine under #40 checks for `FORTNUM_DOMAIN_ERROR`
before propagating a sensitivity and refuses the non-smooth case rather than
returning a wrong number.

## 5. Procedural wrapper

```fortran
subroutine integrate(f, a, b, value, status, epsabs, epsrel, key, ctx)
    procedure(integrate_integrand_t)     :: f
    real(dp),               intent(in)    :: a, b
    real(dp),               intent(out)   :: value
    type(fortnum_status_t), intent(out)   :: status
    real(dp), intent(in), optional :: epsabs, epsrel
    integer,  intent(in), optional :: key
    class(*), intent(in), optional :: ctx
end subroutine integrate
```

`integrate` is the flat call for a caller who wants a value on a finite interval
and does not want to manage the three types. It allocates a local workspace and
result, calls `integrate_qag`, and copies out `value`. Tolerances default to
`epsabs = 0` and `epsrel = 1.0e-8`. The trace stays inside the local result and
is discarded; a caller who needs the trace for a derivative uses the object
form. The name `integrate` is reserved for this primal call and does not collide
with the ad.md §2 derivative names.

## 6. Trace-based derivatives (contract only, #40)

No derivative code ships now. This section fixes names and hook points so #40 is
additive (ad.md §6).

The derivative is taken at the frozen subdivision recorded in
`integrate_result_t` (section 3.3). Once `sub_a`, `sub_b`, and `key` are fixed,
the integral is a fixed weighted sum of GK panel values, so its derivative with
respect to the integrand values is the same weighted sum applied to the
integrand's derivative values. The panel weights are the kernel's `wgk` for the
recorded `key`, evaluated at the panel nodes mapped onto each frozen
subinterval. The derivative product re-walks `sub_a(i)..sub_b(i)` and reuses
`gk_apply`'s node layout; it does not re-adapt and does not move a boundary.

Derivative classification (ad.md §3):

- Active: the integrand values, the parameters carried in `ctx`, and the
  interval endpoints `a`, `b` (a finite endpoint contributes the boundary term
  `f(b) db - f(a) da`).
- Inactive: `epsabs`, `epsrel`, `limit`, `key`, `inf`, the break-point list
  `points` (control knobs and structure selectors, ad.md §3), and `neval`,
  `nsub`, `status` (counts and the status side channel report behavior, they do
  not carry a derivative; ad.md §3 final paragraph).

Entry-point names follow ad.md §2:

- `integrate_qag_jvp`, `integrate_qags_jvp`: forward-mode product. The seed
  tangent flows through the integrand and `ctx`; the routine sums the GK panel
  derivatives over the frozen subdivision and adds the endpoint boundary terms.
- `integrate_qag_vjp`, `integrate_qags_vjp`: reverse-mode product, the same
  frozen-subdivision weighted sum contracted in reverse.
- `integrate_qag_grad`: gradient with respect to `ctx` parameters when the
  integral is the scalar objective.

The extrapolated QAGS value needs care: when `result%extrapolated` is true, the
reported value is `qelg`'s extrapolation of the panel-sum sequence, and its
derivative is the extrapolation applied to the sequence of panel-sum
derivatives, not the derivative of the panel sum alone. The derivative product
reads `result%extrapolated` and reuses the same epsilon recurrence on the
derivative sequence. A non-smooth case (section 4.2, `FORTNUM_DOMAIN_ERROR`) has
no frozen-subdivision derivative; the derivative routine returns that status
rather than a number.

The callbacks carry no global state. A forward product evaluates the integrand's
tangent through the same optional `ctx` channel, and a reverse product
accumulates into a caller-provided adjoint of `ctx`; neither needs a module
pointer (ad.md §6). The primal signatures in sections 4 and 5 are final;
derivatives arrive as new public names beside them.
