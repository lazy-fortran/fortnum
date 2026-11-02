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

## Matrix-free phi-function actions

`fortnum_phi` evaluates a linear combination

\[
w=\sum_{k=0}^{p}\varphi_k(\tau A)v_k,
\qquad
\varphi_k(z)=\sum_{j=0}^{\infty}\frac{z^j}{(j+k)!}.
\]

Here `A` is a real operator, `tau` is a signed real scalar, and `v_k` are
real vectors of equal length. The integer `p` is the highest function order;
`z` is a dimensionless scalar argument. This includes the exponential (`k=0`) and the
functions needed by exponential Rosenbrock methods (`k=1,2,3`). No inverse
of `A` is required, including when `A` is singular or `tau` is zero.

The caller extends `real_linear_operator_t` and supplies its `apply` method.
The implementation builds an Arnoldi basis for the augmented operator

\[
 B=\begin{pmatrix}
 \tau A & v_p/\sigma & \cdots & v_1/\sigma\\
 0 &&K&
 \end{pmatrix},
 \qquad
 x_0=\begin{pmatrix}v_0/\sigma\\e_p\end{pmatrix},
\]

where `K` is the `p` by `p` shift matrix with ones on its first
superdiagonal, `e_p` is its last coordinate vector, and `sigma` is the
largest absolute input entry. The physical part of `sigma*exp(B)*x_0`
is the requested combination. Zero inputs return zero without applying `A`.
The `p=0` case omits the auxiliary block. Input normalization avoids
squaring extremely large or small vector entries in the Euclidean norm.
The full augmented matrix is never assembled.

Arnoldi uses the existing parallel vector kernels and two passes of
classical Gram-Schmidt. The basis grows from `min_dimension` in increments
of four up to `max_dimension`, or stops at an invariant subspace. A
rejected evaluation halves the internal substep. Accepted substeps may
grow by at most a factor of two. Both basis storage and the substep and
rejection budgets are bounded by caller controls. The small projected
exponential uses degree-13 diagonal Padé scaling and squaring with a
LAPACK linear solve, following
[Higham, SIAM J. Matrix Anal. Appl. 26, 1179–1193 (2005)](https://doi.org/10.1137/04061101X).

Two additional projected coordinates integrate the Arnoldi defect once and
twice. The larger magnitude estimates the omitted contribution and reduces
sensitivity to cancellation of the first term. A substep is accepted when
this estimate is at most one tenth of its allocated tolerance. For a
substep occupying fraction `ds` of the interval, that allocation is
`ds*(atol + rtol*scale)`, where `scale` is the larger of the maximum input
vector norm and the candidate result norm. Reported error is the sum of
accepted local estimates in input units. It excludes Padé and floating-point
roundoff, and is an estimate, not a rigorous forward-error bound. Strongly
nonnormal growth can amplify earlier local errors; cancellation in the
requested combination can also make its relative error exceed `rtol`.

This is an augmented-exponential Arnoldi implementation with adaptive
substeps and basis dimensions. It does not implement the particular
cost-minimizing controller of
[Niesen and Wright, Algorithm 919, ACM TOMS 38(3), article 22 (2012)](https://doi.org/10.1145/2168773.2168781).
That paper supplies the related polynomial-inhomogeneity formulation and
adaptive Krylov context.

The action is intended for exponential integration of `y' = f(t,y)`.
A constant mass matrix is not part of this interface. A consumer must
supply an operator in the explicit evolution equation's coordinates.
Exponential Rosenbrock stepping and its embedded time-error controller
are separate from the internal phi-action error control.

Validation compares orders zero through three against the closed-form
Dirichlet Laplacian spectrum and against a quadruple-precision Taylor
scaling-and-squaring exponential of a full nonnormal augmented matrix.
Both cover matrix 1-norms through `1e4` with relative error at most `1e-10`.
Additional cases cover mixed combinations, forced substepping, oscillatory
operators, signed and zero time, singular operators, strided inputs,
large and small vector scales, and explicit failure statuses. These tests
do not establish a stiffness-independent time integrator.
