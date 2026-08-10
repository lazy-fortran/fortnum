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

## Four-state GPU reverse communication

`fortnum_ode_rk54_device` provides allocation-free Cash-Karp 5(4) and
Dormand-Prince 5(4) stepping for four-state accelerator workloads. The caller
asks for a stage state, evaluates its application-specific RHS on the device,
and supplies the derivative. This avoids procedure pointers, polymorphism, and
host callbacks in accelerator regions while leaving the tableau and adaptive
controller in fortnum.

The Dormand-Prince controller includes the FIRM3D-compatible four-component
maximum norm, `error**(-1/3)` update, bounded growth, and minimum-step
acceptance. Absolute tolerances are per component so velocity-like states can
use the same scaling as the external reference. Cash-Karp retains fortnum's RMS
norm and PI controller. Generated stage and embedded-error leaves are pure,
elemental, OpenACC device routines and OpenMP declare-target procedures.

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
