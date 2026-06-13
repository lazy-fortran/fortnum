# ADR: fortnum_ode API

Status: accepted (issue #14, M2.1). Normative for the implementation issues
#15 (primal integrator), #16 (event detection), #17 (procedural wrapper), and
#40 (derivative products).

This ADR is subordinate to `docs/design/ad.md`. The derivative vocabulary,
policy classes, naming convention, and active-argument rules come from that
document; this one cites it by section and fixes the module-specific Fortran.
The default policy class for `ode` is `trace_rule` (ad.md §4): adaptive step
control makes the step sizes data-dependent, so a derivative differentiates the
recorded step schedule with that schedule held fixed.

An implementer should be able to write the solver from this file without
further design decisions. Where this file underspecifies a numeric constant,
the Cash-Karp tableau and the PI controller cited in section 4 fix it.

## 1. Scope and method

One adaptive explicit method: Cash-Karp RK5(4) with embedded error estimate and
PI step-size control. Cash and Karp, ACM TOMS 16 (1990) 201-222. The fifth-order
solution advances; the fourth-order solution gives the local error estimate. The
six stages share node and weight tables; stage `k1` is reused as the next step's
first stage only when the previous step was accepted (FSAL does not apply to
Cash-Karp, so `k1` is recomputed each step).

The state vector `y` has fixed length `neq` for the whole integration. The
solver rejects a system whose size changes mid-integration. The independent
variable advances
monotonically from `t0` toward `t1`; `t1 < t0` integrates backward and is
allowed (the controller works on `sign(t1 - t0)`).

## 2. RHS and event interfaces

The right-hand side is an abstract interface. The caller passes state in, writes
the derivative out, and may thread its own data through an optional unlimited
polymorphic context. No module-level state carries problem data.

```fortran
abstract interface
    subroutine ode_rhs_t(t, y, dydt, ctx)
        import :: dp
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
    end subroutine ode_rhs_t
end interface
```

`y` and `dydt` are assumed-shape and have the same length `neq`. The optional
`ctx` is the only channel for parameters; a caller with parameters defines a
derived type, passes it as `ctx`, and selects its type inside `rhs`. `ctx` is
`intent(in)`: the RHS reads parameters, it does not mutate caller state. A
parameter-free RHS ignores the argument.

The event function returns a scalar whose sign change marks the event:

```fortran
abstract interface
    function ode_event_t(t, y, ctx) result(g)
        import :: dp
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: g
    end function ode_event_t
end interface
```

The same `ctx` convention applies. The event is detected when `g` crosses zero
between two accepted steps; see section 6.

## 3. Derived types

Three types carry the problem, the workspace, and the result. The problem is
read-only input. The workspace is caller-owned scratch, reused across steps and
across calls. The solution holds the recorded trace plus metadata.

### 3.1 Problem

```fortran
type :: ode_problem_t
    procedure(ode_rhs_t), pointer, nopass :: rhs => null()
    real(dp)              :: t0   = 0.0_dp
    real(dp)              :: t1   = 0.0_dp
    real(dp), allocatable :: y0(:)
    real(dp)              :: rtol = 1.0e-6_dp
    real(dp)              :: atol = 1.0e-9_dp
    real(dp)              :: h0   = 0.0_dp     ! 0 => auto initial step
    real(dp)              :: hmin = 0.0_dp     ! 0 => no floor
    real(dp)              :: hmax = 0.0_dp     ! 0 => |t1 - t0|
    integer               :: max_steps = 100000
    procedure(ode_event_t), pointer, nopass :: event => null()  ! null => none
    integer               :: event_direction = ODE_EVENT_ANY
    logical               :: terminal_event  = .true.
    real(dp)              :: event_tol = 1.0e-12_dp
end type ode_problem_t
```

`neq = size(y0)`. A null `event` pointer means no event detection; the event
fields are then ignored. `event_direction` is one of the module parameters

```fortran
integer, parameter :: ODE_EVENT_RISING  =  1   ! g - to +
integer, parameter :: ODE_EVENT_FALLING = -1   ! g + to -
integer, parameter :: ODE_EVENT_ANY     =  0
```

A procedure pointer holds the RHS, not an allocatable, because a procedure is
not data; `nopass` keeps the interface free of a passed-object dummy. The
defaults give a usable problem once `rhs`, `t0`, `t1`, and `y0` are set.

### 3.2 Workspace

```fortran
type :: ode_workspace_t
    integer               :: neq = 0
    real(dp), allocatable :: k1(:), k2(:), k3(:), k4(:), k5(:), k6(:)
    real(dp), allocatable :: ytmp(:), y5(:), y4(:), yerr(:)
end type ode_workspace_t
```

The six stage arrays and the temporaries are sized `neq`. The integrator
allocates them once, on the first call or whenever `workspace%neq` does not
match the problem, and reuses them on every later step and call. No allocation
happens inside the step loop. The caller owns the workspace and may keep it
between integrations to avoid reallocation; passing a default-initialized
`ode_workspace_t` is correct and triggers a one-time allocation.

### 3.3 Solution and trace

```fortran
type :: ode_solution_t
    integer               :: nsteps     = 0   ! accepted steps
    integer               :: nrejected  = 0
    integer               :: nfev       = 0   ! RHS evaluations
    real(dp), allocatable :: t(:)             ! t(0:nsteps), t(1)=t0
    real(dp), allocatable :: y(:,:)           ! y(neq, nsteps+1)
    real(dp), allocatable :: h(:)             ! accepted step sizes, h(nsteps)
    real(dp), allocatable :: err(:)           ! scaled error estimate per step
    logical               :: event_found = .false.
    real(dp)              :: t_event = 0.0_dp
    real(dp), allocatable :: y_event(:)       ! state at the located event
    type(fortnum_status_t) :: status
end type ode_solution_t
```

`t` and the columns of `y` record the accepted mesh: `t(1) = t0`, `y(:,1) = y0`,
and `t(nsteps+1)` is the last accepted point (either `t1` or the terminal event
time). `h(i)` is the step that advanced from `t(i)` to `t(i+1)`; `err(i)` is the
scaled error estimate that accepted that step. This recorded mesh is the frozen
schedule a `trace_rule` derivative differentiates (ad.md §1, §4): the step sizes
are fixed data, and the sensitivity solves the variational system on that mesh.

The arrays grow geometrically as steps are accepted; the integrator trims them
to the final length before returning. `status` mirrors the `status` argument so
a caller that keeps only the solution still sees the outcome.

## 4. Primal integrate routine

```fortran
subroutine ode_integrate(problem, workspace, solution, status)
    type(ode_problem_t),    intent(in)    :: problem
    type(ode_workspace_t),  intent(inout) :: workspace
    type(ode_solution_t),   intent(inout) :: solution
    type(fortnum_status_t), intent(out)   :: status
end subroutine ode_integrate
```

The routine integrates from `t0` to `t1` with adaptive Cash-Karp RK5(4). Each
step computes the six stages, forms the fifth-order update `y5` and the
fourth-order update `y4`, and takes the error estimate `yerr = y5 - y4`. The
error norm is the RMS over components scaled by `atol + rtol*max(|y|,|y5|)`; a
step is accepted when the norm is at most one.

Step-size control is a PI controller on the accepted-step error history:
`h_new = h * fac * err_norm**(-alpha) * err_prev**(beta)`, with `alpha = 0.7/5`,
`beta = 0.4/5`, safety `fac = 0.9`, and growth clamped to `[0.2, 5.0]`. On the
first step and after any rejection the controller falls back to the standard
`I` rule (`beta = 0`). A rejected step shrinks `h` and retries the same `t`
without recording a trace entry; only accepted steps append to the trace and
increment `nsteps`. When `h0 = 0` the initial step is chosen from the scaled
RHS magnitude at `t0` (Hairer, Norsett, Wanner, "Solving ODEs I", II.4,
starting-step estimate).

The step is clipped so the integration lands exactly on `t1`. `hmax` caps the
step; `hmin` floors it, and a step forced below `hmin` reports
`FORTNUM_CONVERGENCE_ERROR`. Exceeding `max_steps` reports the same code. On
success `status%code = FORTNUM_OK`. Invalid input (`neq < 1`, non-positive
tolerances, null `rhs`, unallocated `y0`) reports `FORTNUM_DOMAIN_ERROR` before
any stepping.

Derivative classification (ad.md §3):

- Active: `y0`, the RHS parameters carried in `ctx`, and the recorded trace
  (`t`, `y`, `h`). A `trace_rule` sensitivity differentiates the frozen mesh.
- Inactive: `rtol`, `atol`, `h0`, `hmin`, `hmax`, `max_steps`, `event_direction`
  (control knobs select behavior, ad.md §3), and `nsteps`, `nrejected`, `nfev`,
  `err`, `status` (counts and the status side channel report behavior, they do
  not carry a derivative; ad.md §3 final paragraph).

## 5. Procedural wrapper

```fortran
subroutine ode_solve(rhs, t0, t1, y0, t_out, y_out, status, rtol, atol)
    procedure(ode_rhs_t)               :: rhs
    real(dp),               intent(in) :: t0, t1
    real(dp),               intent(in) :: y0(:)
    real(dp), allocatable, intent(out) :: t_out(:)
    real(dp), allocatable, intent(out) :: y_out(:,:)
    type(fortnum_status_t), intent(out) :: status
    real(dp), intent(in), optional     :: rtol, atol
end subroutine ode_solve
```

`ode_solve` is the flat call for callers who do not want to manage the three
types. It builds an `ode_problem_t`, allocates a local workspace and solution,
calls `ode_integrate`, and moves the trace mesh into `t_out`, `y_out`. Tolerances
default as in `ode_problem_t`. Events are not exposed here; a caller who needs an
event uses the object form.

The name `ode_solve` is reserved for the primal convenience call and does not
collide with the ad.md §2 derivative names. Sensitivity wrappers added under #40
keep the primal signature untouched and take the `ode_*_jvp`, `ode_*_vjp`,
`ode_*_grad` forms of section 7. `ode_solve` never grows a derivative argument.

## 6. Event detection (#16)

When `problem%event` is associated, the integrator evaluates `g(t, y, ctx)` at
each accepted point. A sign change between two consecutive accepted points,
consistent with `event_direction`, brackets a root. The root is located in `t`
on the step's dense interpolant (the Cash-Karp fourth-order interpolant, Hairer
II.6) by Brent's method to `event_tol`. The located time goes to
`solution%t_event`, the interpolated state to `solution%y_event`, and
`solution%event_found` is set. With `terminal_event = .true.` the integration
stops at the event and the last trace point is the event point; otherwise
integration continues and only the first event is reported.

Event derivative rule (ad.md §1, `trace_rule`): an event time and event state
are differentiable only at a smooth, transversal crossing. Transversal means
`dg/dt = grad_y g . f` is nonzero at the root, so the implicit function theorem
gives `d t_event` from `d g = 0`. When the crossing is tangential
(`dg/dt` within `event_tol` of zero) or `g` is non-smooth at the root, the event
time is not a differentiable function of `y0` or the parameters; the integrator
reports `FORTNUM_DOMAIN_ERROR` with a message naming the non-transversal event,
and the trace up to the event stays valid as a primal. A derivative product
under #40 checks transversality before propagating an event sensitivity.

## 7. Forward, reverse, and checkpointing (contract only, #40)

No derivative code ships now. This section fixes names and hook points so #40 is
additive (ad.md §6).

- Forward sensitivity: `ode_integrate_jvp` solves the variational system
  `dot(S) = (df/dy) S + df/dp` on the frozen trace mesh from section 3.3,
  seeded by the input tangent. It reuses `ode_problem_t` and a sensitivity
  workspace; the primal trace is the schedule it integrates against.
- Reverse/adjoint: `ode_integrate_vjp` integrates the adjoint
  `dot(lambda) = -(df/dy)^T lambda` backward over the recorded mesh, then
  contracts against `df/dp`. It consumes `solution%t`, `solution%y`,
  `solution%h` as the fixed backward schedule.
- Checkpointing: the adjoint needs the forward state at each mesh point. The
  default keeps the full `solution%y`; a checkpointed variant stores a strided
  subset and re-integrates between checkpoints (Griewank-Walther revolve). The
  checkpoint stride is an inactive control on the sensitivity workspace, not on
  `ode_problem_t`, so the primal types in section 3 do not change when #40 lands.

These names follow ad.md §2 (`foo_jvp`, `foo_vjp`). The primal `ode_integrate`
and `ode_solve` signatures in sections 4 and 5 are final; derivatives arrive as
new public names beside them.
