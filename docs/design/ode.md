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

Each attempt makes one implicit time-step solve. Its Newton initial guess is
the quadratic Hermite predictor through `y_{n-1}`, `y_n` and the retained
endpoint derivative `d_n`:

\[
y^P_{n+1}=y_n+h d_n+r\{h d_n-r(y_n-y_{n-1})\}.
\]

For exact starting data and fixed positive `r`, Taylor expansion gives the
predictor and BDF2 local solution errors

\[
y^P_{n+1}-y(t_n+h)=-p(r)h^3 y'''(t_n)+O(h^4),\quad
y_{n+1}-y(t_n+h)=b(r)h^3 y'''(t_n)+O(h^4),
\]
\[
p(r)=\frac{1+r}{6r},\qquad
b(r)=\frac{(1+r)^2}{6r(1+2r)}.
\]

Here the derivative defect of BDF2 is
`-h*(h+k)*y'''/6`; negation and division by its leading shifted coefficient
`alpha=(1+2r)/((1+r)*h)` gives `b(r)`. Eliminating the common third derivative
is Milne's device:

\[
e_{n+1}=C(r)(y_{n+1}-y^P_{n+1}),\qquad
C(r)=\frac{b(r)}{b(r)+p(r)}=\frac{1+r}{2+3r}.
\]

This estimates the BDF2 local solution error of order `h**3`, rather than the
derivative residual or a BDF1 error. The derivation above is for this Hermite
predictor; an Adams--Bashforth predictor has different constants. For the
standard BDF and multistep error-estimation framework, see E. Hairer and
G. Wanner, *Solving Ordinary Differential Equations II: Stiff and
Differential-Algebraic Problems*, second revised edition, Springer (1996),
[Chapter V](https://doi.org/10.1007/978-3-642-05221-7).
The formulas are retained in the FortSym generator
`tools/codegen/app/gen_bdf2.f90` and emitted into `src/generated/`.

At startup, the predictor is explicit Euler `y^P_1=y_0+h*d_0`. The estimate
`(y_1-y^P_1)/2` is asymptotically the BDF1 local solution error of order
`h**2`; the accepted value is the single backward Euler step. There are no
embedded or half-step implicit solves. An optional `derivative` argument to
`bdf_init` supplies a consistent initial `d_0`. Otherwise the first attempt
obtains it once from the linear equation `M(t_0,y_0)*d_0=-R(t_0,y_0)` using
unpreconditioned FGMRES with normalized right-hand side and relative tolerance
`1e-8`. Its mass action uses `F(t_0,y_0,v)-F(t_0,y_0,0)`, since `F` is affine
in `ydot`. This setup does not solve another implicit time step. Supplying
`d_0` avoids that linear setup, useful for a difficult mass matrix. Subsequent
derivatives come directly from the accepted implicit solve, including across
continuation and direction changes.

The error indicator is the root-mean-square norm of `e_{n+1}` divided
componentwise by `atol + rtol*max(abs(y_n),abs(y_{n+1}))`. It is an asymptotic
local estimate for smooth resolved solutions and sufficiently accurate Newton
solves, not a rigorous bound on the global trajectory error or a uniform stiff
error bound.

An indicator above one rejects the attempt. The elementary controller uses
`h_new = abs(h)*0.6*error**(-1/3)` for BDF2 and exponent `-1/2` for startup.
The multiplier is limited to `[0.2,2]`; zero error selects two. Every attempted
step is also capped at twice the previous accepted step magnitude, below the classical
BDF2 zero-stability ratio `1+sqrt(2)`, even after endpoint clipping or a retry.
Nonlinear failure reduces the step by four. `adaptive=.false.` disables the
error test, retains the proposed step magnitude after acceptance, and still
rejects failed nonlinear solves. Set `initial_step=max_step` for
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
Krylov iterations, and evaluations for the current call. `implicit_solves`
counts all attempted implicit time-step solves, including rejected ones;
`mass_iterations` counts the separate initial mass solve's Krylov iterations.
`last_error` is the last error indicator. `last_failure` describes the latest
solve as
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
and periodic diffusion Fourier modes. Exact-start single steps check the local
error estimate at ratios `0.5`, `1` and `2` in both time directions. A smooth
nonlinear trajectory checks cube-root tolerance scaling. Robertson on `[0,1]`
retains the error/rejection regression, while `[0,4]` checks achieved error
within a factor ten of tolerance and approximately one implicit solve per
accepted step. Diffusion includes a slow constant mode
and excited fast modes, with stiffness ratios starting at `1e8` and increasing
under refinement; the test preconditioner uses a dense factorization. It tests
bounded time-step counts, not the scalability of that preconditioner.

## Autonomous exponential Rosenbrock integration

`fortnum_ode_exprb` integrates `y' = f(y)` with exponential Rosenbrock–Euler
(order two, `EXPRB_EULER`) or the embedded third/second-order pair `EXPRB32`
(the default). These are the autonomous methods of Hochbruck, Ostermann and
Schweitzer, [SIAM J. Numer. Anal. 47, 786–803 (2009), Sections 2 and 5–6](https://doi.org/10.1137/080717717).
The caller extends `exprb_problem_t`, which inherits `real_linear_operator_t`:

| Method | Contract |
| --- | --- |
| `evaluate(y,f,status)` | Evaluate `f(y)` without changing the retained Jacobian. Zero status means success. |
| `linearize(y,status)` | Retain `y` and prepare its Jacobian. Zero status means success. |
| `apply(input,output)` | Apply the retained Jacobian to `input`; no assembly is required. |

For non-autonomous equations, append time as a state with derivative one.
The caller's Jacobian action must then include the derivative with respect to
that appended time. There is no separate time-derivative callback. This is an
ordinary differential equation interface, without a mass matrix; conversion
from a mass-matrix residual belongs to the caller.

For accepted state `y`, step `h`, and Jacobian `J = f'(y)`, the stages are

\[
U=y+h\phi_1(hJ)f(y),\qquad D=f(U)-f(y)-J(U-y),
\]
\[
y_{\rm new}=U+2h\phi_3(hJ)D.
\]

Here `phi_k` is the entire function defined in the
[phi-action contract](../api.md#linear-algebra). Euler accepts `U`; exprb32
accepts `y_new`. Each stage makes one `fortnum_phi` action with its coefficient
vector already scaled by `h` or `2h`. Fixed-step Euler needs only the first
stage. Adaptive Euler also computes the correction for its error estimate.
The Jacobian is frozen for both stages and refreshed on every attempt.

```fortran
type(exprb_options_t) :: options
type(exprb_state_t) :: state
type(exprb_report_t) :: report
type(fortnum_status_t) :: status
! problem extends exprb_problem_t.
call exprb_init(state, initial_time, initial_y, options, status)
if (status%code == FORTNUM_OK) then
    call exprb_integrate_to(problem, state, final_time, options, report, status)
end if
```

The caller owns options, state and report. `state%y` and `state%t` remain the
last accepted endpoint on failure. `next_step` is a positive magnitude;
forward, backward, clipped-endpoint and continuation calls are supported.
Reinitialize after changing the state discontinuously.

The estimate is the correction `2h phi_3(hJ) D`, computed directly to avoid
subtracting nearby endpoints. Its weighted root-mean-square norm uses
`atol + rtol*max(abs(y),abs(accepted_candidate))`. An indicator above one
rejects the attempt. Both methods use the cube-root controller
`0.8*error**(-1/3)`, limited to `[0.2,2]`, with two for zero error. The scalar
controller is retained in `tools/codegen/app/gen_exprb.f90` and emitted by
FortSym. `adaptive=.false.` disables the error test and keeps the proposed
step magnitude; choose `initial_step=max_step` for constant-step studies.

A failed callback, nonfinite output or failed phi action rejects the attempt
and reduces its magnitude by four. `max_steps` bounds all attempts per call.
Invalid controls or state return `FORTNUM_DOMAIN_ERROR`; exhausted attempts,
a step below `min_step` (including a clipped final step), or an inability to
advance floating-point time return `FORTNUM_CONVERGENCE_ERROR`.

`options%phi` holds independent phi-action controls, defaulting to
`phi_options_t()` with relative tolerance `1e-12`. They are not automatically
relaxed with the time tolerance. The embedded estimate measures time error,
not phi-action error, so callers must keep the phi error small enough for
their accuracy target. Neither estimate is a rigorous global bound, and no
uniform accuracy or work bound for arbitrary stiff problems is claimed.

Reports reset per call. They count accepted steps, error/failure rejections,
function evaluations, linearizations, Jacobian actions (including the explicit
remainder action), phi calls, phi substeps/rejections, and maximum Krylov
dimension. `last_failure` is `EXPRB_SOLVED`, `EXPRB_EVALUATION`,
`EXPRB_LINEARIZE` or `EXPRB_PHI_FAILURE`. `callback_status` retains the latest
nonzero callback status; `phi_status` and `last_error` describe the latest
phi action and completed error test. No events, device execution or
sensitivity trace are provided.

The registered `fortnum_ode_exprb` test checks nonlinear manufactured orders,
a Dirichlet reaction–diffusion system made autonomous by a clock state,
stiff van der Pol tolerance control against refined classical RK4, and
periodic diffusion with excited fast Fourier modes under three refinements.
It also compares both exponential methods and BDF2 with independent oracles.
The diffusion gate bounds accepted steps, not Krylov work or memory; small
meshes permit full-space Krylov projections. BDF comparisons use dense
factorization as a test preconditioner, not a scalable PDE solver.

### Measured comparison

The `fortnum_ode_exprb` test calibrates tolerances against the independent
endpoint oracle to reach a maximum absolute error between `1e-7` and `1e-6`
on each nonlinear problem. Calibration is excluded from timing. Periodic
linear diffusion uses the same `1e-6` error ceiling; the exponential methods
are exact in time there, up to phi-action and floating-point error, so their
errors cannot be matched to BDF2 by adjusting the time tolerance.

The following CPU samples use GNU Fortran 16.2, Release, on an Apple M3 Ultra,
with one OpenMP and one OpenBLAS thread. Time is the median of three complete
solves, including setup and allocation; brackets give the minimum and maximum.
Function counts include the supplied initial BDF derivative. Jacobian actions
include the BDF preconditioner's unit-vector columns and shifted operator
actions. All BDF preconditioners are dense factors; these small tests do not
measure production preconditioner scalability. The state and time variables
in these synthetic problems are dimensionless.

| Problem | Method | Maximum error | Steps | Function evaluations | Jacobian actions | Wall time, ms [min, max] |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| Scalar decay | Rosenbrock–Euler | 6.04e-07 | 143 | 286 | 1001 | 1.106 [1.095, 1.108] |
| Scalar decay | exprb32 | 5.29e-07 | 21 | 42 | 147 | 0.160 [0.159, 0.163] |
| Scalar decay | BDF2 | 6.03e-07 | 354 | 1057 | 2808 | 0.910 [0.901, 0.915] |
| Dirichlet reaction–diffusion | Rosenbrock–Euler | 5.99e-07 | 209 | 418 | 4389 | 7.713 [7.669, 7.978] |
| Dirichlet reaction–diffusion | exprb32 | 6.19e-07 | 20 | 40 | 480 | 1.114 [1.103, 1.127] |
| Dirichlet reaction–diffusion | BDF2 | 8.56e-07 | 160 | 472 | 5909 | 1.220 [1.210, 1.239] |
| Van der Pol | Rosenbrock–Euler | 6.01e-07 | 7946 | 15898 | 71541 | 80.390 [80.219, 81.362] |
| Van der Pol | exprb32 | 4.13e-07 | 911 | 1828 | 8226 | 9.638 [9.197, 9.845] |
| Van der Pol | BDF2 | 5.92e-07 | 15055 | 45169 | 150555 | 44.238 [43.867, 45.051] |
| Periodic diffusion | Rosenbrock–Euler | 5.35e-14 | 10 | 20 | 56 | 0.122 [0.121, 0.124] |
| Periodic diffusion | exprb32 | 2.33e-14 | 10 | 20 | 56 | 0.122 [0.121, 0.123] |
| Periodic diffusion | BDF2 | 6.01e-07 | 1459 | 4396 | 55613 | 8.737 [8.736, 9.023] |

Scalar decay solves `y'=-y^2`, `y(0)=1`, to time `0.5`. Reaction–diffusion
uses 15 interior points on `[0,1]`, zero Dirichlet boundaries, unit diffusivity,
quadratic reaction and manufactured forcing for
`y_i(t)=sin(pi*x_i)/(1+t)`, with an appended clock, to time `0.5`.
Van der Pol uses stiffness parameter `mu=100`, initial state `(2,0)` and
endpoint `100`; its independent classical fourth-order Runge–Kutta reference
uses 1.6 million uniform steps and is checked against 0.8 million steps
(endpoint difference `5.64e-10`). Diffusion uses 16 periodic points, decay
rate one, diffusion coefficient `1e8/(4*16^2)`, an initial constant plus cosine
and alternating modes, and endpoint one. Its exact endpoint differs from
the decayed constant by exponentially negligible fast modes.

The order tests measure Euler slopes `2.004` and `2.011`, and exprb32 slopes
`3.006` and `2.996`, on scalar decay and reaction–diffusion respectively.
For van der Pol, tolerances `1e-3` through `1e-7` give endpoint errors
`2.01e-5`, `4.73e-6`, `6.00e-7`, `1.27e-7`, `1.45e-8`.
Diffusion grids of 16, 32 and 64 points have stiffness ratios starting at
`1e8` and increasing quadratically. All take seven accepted steps, with
maximum Krylov dimensions 12, 20, 36 and total phi substeps 8, 8, 9.
These observations concern these oracles, not arbitrary stiff systems.

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
