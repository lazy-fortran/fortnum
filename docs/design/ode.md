# ODE contract

Status: implemented.

## Main solver

`fortnum_ode` integrates

\[
\dot y=f(t,y,p)
\]

with an adaptive Cash-Karp RK5(4) method and PI step-size control. The caller
owns `ode_problem_t`, `ode_workspace_t`, and `ode_solution_t`.

```fortran
problem%rhs => rhs
problem%t0 = t0
problem%t1 = t1
problem%y0 = y0
call ode_integrate(problem, workspace, solution, status)
```

`ode_solve` builds this state for a simple allocating call.

The RHS interface is:

```fortran
subroutine ode_rhs_t(t, y, dydt, ctx)
    real(dp), intent(in) :: t
    real(dp), intent(in) :: y(:)
    real(dp), intent(out) :: dydt(:)
    class(*), intent(in), optional :: ctx
end subroutine ode_rhs_t
```

## Problem and solution state

`ode_problem_t` contains:

- RHS callback
- start and end times
- initial state
- relative and absolute tolerances
- initial, minimum, and maximum step sizes
- maximum accepted/rejected step budget
- optional event callback and direction
- terminal-event flag and event tolerance

`ode_workspace_t` owns stage and error arrays.

`ode_solution_t` records accepted times, states, signed step sizes, local error
norms, counters, status, and an optional terminal event. The trace contains the
initial point followed by accepted endpoints.

The solver supports forward and backward time integration. It clips the final
step to the requested endpoint. Domain errors cover missing callbacks, empty
states, invalid tolerances, and invalid controls. Step exhaustion or a forced
step below `hmin` returns a convergence error.

## Other methods

The package also exposes:

- DOP853 step and driver interfaces
- stateful RK8PD evolution
- stateful Adams DDEABM integration
- stateful VODE-style integration
- `ode_at` for requested output times

Each method owns its documented state. They do not share hidden solver globals.

## Matrix-free implicit BDF integration

`fortnum_ode_bdf` integrates `F(t,y,ydot) = M(t,y) ydot + R(t,y) = 0`
with a nonsingular mass matrix `M`, using first-order backward differentiation
(BDF1, backward Euler) at startup and variable-step BDF2 thereafter. `R` is
the spatial or physical residual and `ydot` is the time derivative. The caller
extends `bdf_problem_t`, which reuses `real_linear_operator_t` from
`fortnum_krylov`, and implements these methods:

| Method | Contract |
| --- | --- |
| `evaluate(t,y,ydot,residual,status)` | Evaluate `F` without changing the current linearization. Zero `status` means success. |
| `linearize(t,y,ydot,alpha,status)` | Retain the point and prepare the shifted Jacobian and preconditioner. Zero `status` means success. |
| `apply(input,output)` | Apply `F_y input + alpha F_ydot input` at the retained point. |
| `precondition(r,z,status)` | Optional right inverse action for the shifted Jacobian; the default is the identity. Nonzero `status` rejects the attempt. |

The shifted Jacobian includes `(dM/dy input) ydot` for a state-dependent mass
matrix, as well as `alpha M input` and `R_y input`. Omitting the mass derivative
changes Newton's equation. The integrator never assembles or inverts `M` and
never assumes a symmetric Jacobian. Callback state belongs to the problem
object; each trajectory owns a separate `bdf_state_t`.

```fortran
type(bdf_options_t) :: options
type(bdf_state_t) :: state
type(bdf_report_t) :: report
type(fortnum_status_t) :: status
! problem is a caller-defined extension of bdf_problem_t.
options%rtol = 1.0e-5_dp
options%atol = 1.0e-8_dp
options%initial_step = 1.0e-3_dp
call bdf_init(state, 0.0_dp, initial_y, options, status)
if (status%code == FORTNUM_OK) then
    call bdf_integrate_to(problem, state, final_time, options, report, status)
end if
```

`state%y` and `state%t` are the last accepted endpoint. Continuation to another
target retains the history; changing direction restarts with BDF1. Endpoint
clipping lands on the requested time without interpolation. `next_step` stores
a positive proposed magnitude, while `previous_step` is signed. Do not modify
the accepted state or history between calls; reinitialize after a discontinuity.

For current step `h`, previous step `k` and ratio `r=h/k`, the BDF2 derivative is

\[
\dot y_{n+1}=\frac{1+2r}{(1+r)h}(y_{n+1}-y_n)
             -\frac{r^2}{(1+r)h}(y_n-y_{n-1}).
\]

Step growth is capped at two, below the classical BDF2 zero-stability ratio
`1+sqrt(2)`. Each attempt solves both the BDF2 equation and an embedded BDF1
equation at the same endpoint. The error indicator is the root-mean-square
norm of their difference divided componentwise by
`atol + rtol*max(abs(y_n),abs(y_{n+1}))`. At startup it is twice the difference
between one BDF1 step and two half BDF1 steps. The accepted startup value is
the single full step. The indicator estimates the lower-order local error,
which is proportional to `h**2` for smooth solutions; it is conservative for
BDF2, rather than an estimate of its `h**3` local error. It controls local
differences, not a rigorous bound on the global trajectory error.

An indicator above one rejects the attempt. The elementary controller uses
`0.9/sqrt(error)`, limited to `[0.2,2]`; zero error selects two. Nonlinear
failure reduces the step by four. `adaptive=.false.` disables the embedded
solve and error test, retains the proposed step magnitude after acceptance,
and still rejects failed nonlinear solves. Set `initial_step=max_step` for
constant-step convergence studies. `min_step` is enforced even on a clipped
final step; a step that cannot advance floating-point time also fails.
`max_steps` bounds all attempts per call, including rejections. Invalid data
return `FORTNUM_DOMAIN_ERROR`; underflow, minimum-step or budget exhaustion
return `FORTNUM_CONVERGENCE_ERROR`, retaining the last accepted endpoint.

The inner inexact Newton solve uses fortnum's existing right-preconditioned
FGMRES, residual backtracking and
[Eisenstat--Walker choice 2](https://doi.org/10.1137/0917003), with multiplier
0.9, exponent 1.5 and its safeguard against a sudden forcing decrease. The
initial forcing is 0.1, the upper limit 0.9, and `forcing_min` (default
`1e-4`, allowed in `(0,0.1]`) prevents oversolving near the linear rounding
floor. The normalized linear right-hand side preserves relative forcing
despite FGMRES's absolute stopping floor. Newton converges on a zero residual
or a full correction whose weighted root-mean-square norm is at most
`newton_tolerance` (default 0.03), using the same state weights as the error
test. A converged correction need not decrease a residual already at roundoff.

The report counts accepted steps, error and nonlinear rejections, Newton and
Krylov iterations, and evaluations for the current call. `last_error` is the
last error indicator. `last_failure` describes the latest implicit solve as
`BDF_SOLVED`,
`BDF_NEWTON_LIMIT`, `BDF_LINEAR_FAILURE`, `BDF_EVALUATION`, `BDF_LINE_SEARCH`
or `BDF_LINEARIZE`; it is reset by a later successful implicit solve.
`callback_status` and `linear_status` retain the last callback and FGMRES
status. Nonfinite callback output is a failure. Rejections change the proposed
step, never the accepted solution or multistep history. This interface is
host-side and primal-only, without event location or a sensitivity trace.

`fortnum_ode_bdf` tests use exact stiff triangular-system and manufactured
nonlinear trajectories (including a state-dependent mass), polynomial
exactness on nonuniform grids, a refined dense trapezoidal Robertson reference,
and periodic diffusion Fourier modes. Diffusion includes a slow constant mode
and excited fast modes, with stiffness ratios starting at `1e8` and increasing
under refinement; the test preconditioner uses a dense factorization. It tests
bounded time-step counts, not the scalability of that preconditioner.

## Four-state GPU reverse communication

`fortnum_ode_rk54_device` provides allocation-free Cash-Karp 5(4) and
Dormand-Prince 5(4) stepping for four-state accelerator workloads. The caller
asks for a stage state, evaluates its application-specific RHS on the device,
and supplies the derivative. This avoids procedure pointers, polymorphism, and
host callbacks in accelerator regions while leaving the tableau and adaptive
controller in fortnum.

The compatibility Dormand-Prince mode includes the FIRM3D-compatible
four-component maximum norm, `error**(-1/3)` update, bounded growth, and
minimum-step acceptance. The tuned Dormand-Prince mode uses the same tableau,
norm, FSAL reuse, and minimum-step policy with the fifth-order first/retry
controller and Fortnum's PI history. Absolute tolerances are per component so
velocity-like states can use the same scaling as the external reference.
Cash-Karp retains Fortnum's RMS norm and PI controller. Generated stage and
embedded-error leaves are pure, elemental, OpenACC device routines and OpenMP
declare-target procedures.

## Events

An event callback returns a scalar residual \(g(t,y)\). Direction is
`ODE_EVENT_RISING`, `ODE_EVENT_FALLING`, or `ODE_EVENT_ANY`. Accepted steps are
scanned for a qualifying sign change, then the event is localized to
`event_tol`.

`ode_event_time_jvp` differentiates the defining event equation. For a
transversal event,

\[
\dot \tau =
-\frac{g_y\dot y_{\mathrm{fixed}}+g_p\dot p}
       {g_t+g_y f}.
\]

The denominator must be nonzero at the requested reliability threshold.
`ode_event_state_jvp` adds the event-time movement
\(f(\tau)\dot\tau\) to the fixed-time state tangent.

## Forward sensitivity

`ode_integrate_jvp` advances

\[
\dot s=f_y s+f_p\dot p
\]

over the accepted primal step schedule. The caller supplies `ode_var_rhs_t`,
which returns the contracted right-hand side for the requested direction.
`s0` carries the initial-state tangent.

`ode_integrate_jvp_multi` carries several directions through one replay. The
columns of `s0` are the seeds, and its `ode_var_rhs_multi_t` callback returns
all contracted right-hand sides of a stage in one call, so a caller can build
\(f_y\) and \(f_p\) once per stage instead of once per direction. Column \(k\)
of the result equals `ode_integrate_jvp` for direction \(k\).

The schedule, tolerances, status, and event choices are inactive. The primal
stage states and signed step sizes are replayed.

## Reverse sensitivity

`ode_integrate_vjp` walks the discrete Cash-Karp maps backward. Its callback
applies the transpose state Jacobian. The returned cotangent is with respect to
the initial state.

`ode_integrate_parameter_vjp` also calls `ode_param_vjp_t` at each stage. That
callback adds \(f_p^T\bar f\) into the caller's parameter layout.

Two memory candidates are available:

- `ode_integrate_vjp_checkpointed` stores selected states and recomputes each
  segment
- `ode_integrate_vjp_recomputed` stores the time/step schedule and initial
  state, then recomputes prefixes

The full trace, checkpoint, and recomputation candidates implement the same
discrete transpose. Complete-workload wall clock and peak memory select among
them.

## Continuous and discrete meaning

The JVP recurrence is a numerical approximation to the continuous variational
equation. At a fixed accepted schedule, it is also the exact tangent of the
frozen Cash-Karp step composition.

The VJP is the exact transpose of that frozen discrete composition. A
continuous adjoint and a discrete adjoint need not agree at finite step size.
Callers must choose the contract required by their objective.

## Implicit stages

For a stage defined by

\[
r(z,p)=z-b-\alpha f(t,z,p)=0,
\]

`ode_implicit_stage_jvp` solves

\[
(I-\alpha f_z)\dot z=\dot b+\alpha f_p\dot p.
\]

`ode_implicit_stage_vjp` solves the transposed system and returns cotangents for
the base state and active parameters. One factorization is reused across
multiple directions or cotangents.

## Validation

Tests compare against closed-form trajectories, matrix exponentials,
manufactured stiff systems, finite differences of the frozen map, adjoint
identities, event-time formulas, implicit residual equations, and refinement
of continuous sensitivities.
