# Adaptive integration contract

Status: implemented.

## Scope

`fortnum_integrate` implements globally adaptive Gauss-Kronrod drivers with
caller-owned state:

| Routine | Behavior |
| --- | --- |
| `integrate_qag` | adaptive bisection with selectable GK rule |
| `integrate_qags` | GK21 plus Wynn epsilon extrapolation |
| `integrate_qagp` | QAGS initialized at caller breakpoints |
| `integrate_qagiu` | transformed semi-infinite or doubly infinite interval |
| `integrate` | convenience allocation around QAG |

`fortnum_integrate_gk` owns the fixed panel rules. The adaptive module owns
interval selection, error ordering, extrapolation, and the accepted trace.

## Integrand

```fortran
function integrate_integrand_t(x, ctx) result(value)
    real(dp), intent(in) :: x
    class(*), intent(in), optional :: ctx
    real(dp) :: value
end function integrate_integrand_t
```

The optional context carries immutable caller data. The module stores no global
callback pointer.

## State

`integrate_workspace_t` stores interval bounds, estimates, errors, ordering,
levels, and breakpoint flags. `limit` is the maximum active interval count.

`integrate_epstab_t` stores the Wynn epsilon table and its short result history.

`integrate_result_t` stores:

- value and absolute error estimate
- function-evaluation and subinterval counts
- accepted panel bounds, results, and errors
- Gauss-Kronrod rule key
- extrapolation flag
- primal status

The panel arrays are the frozen derivative trace.

## Input and status rules

Bounds and tolerances must be finite. The driver accepts a positive absolute
tolerance or a relative tolerance above its floating-point floor. Rule keys and
workspace limits must be supported.

`integrate_qagp` sorts interior breakpoints. Duplicate points and points at or
outside the endpoints are domain errors.

`integrate_qagiu` uses:

| `inf` | Interval |
| ---: | --- |
| `1` | `[bound, +infinity)` |
| `-1` | `(-infinity, bound]` |
| `2` | `(-infinity, +infinity)` |

The transformed coordinate is \(t\in(0,1]\). Non-finite integrand values,
invalid inputs, and unusable traces return a domain error. Exhausted
subdivision or tolerance failure returns a convergence error.

## Analytical integral products

For fixed bounds,

\[
\dot I = \int_a^b \dot f(x)\,dx.
\]

`integrate_fixed_jvp` evaluates this integral. Active lower and upper bounds
use the Leibniz terms:

\[
\dot I =
\int_a^b \dot f(x)\,dx
- f(a)\dot a
+ f(b)\dot b.
\]

`integrate_moving_lower_jvp` and `integrate_moving_upper_jvp` expose the two
boundary terms separately.

## Frozen-trace products

`integrate_qag_jvp`, `integrate_qags_jvp`, and `integrate_qagp_jvp` replay the
accepted panels and the recorded Gauss-Kronrod rule on the supplied contracted
integrand tangent. They do not re-run adaptive decisions.

`integrate_qagiu_jvp` replays the trace in transformed coordinates and includes
the transformation Jacobian. It supports the single-trace `inf=-1` and `inf=1`
cases. A doubly infinite primal stores only one of its two internal traces, so
`inf=2` is a domain error for this derivative interface.

A non-successful primal trace has no valid frozen linearization. Derivative
routines propagate its status without forming a product.

The frozen-trace contract differentiates the accepted numerical map. A caller
that needs the derivative of the continuous integral can use differentiation
under the integral sign and its own derivative integrand.

## Validation

The test suite covers exact polynomial integrals, singular and infinite
reference tables, Leibniz boundary terms, finite-difference directions, frozen
trace replay, and trace-change detection.
