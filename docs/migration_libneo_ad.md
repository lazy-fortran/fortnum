# libneo derivative migration

Status: current mapping guide.

Use the product required by the downstream algorithm. Keep the same active
arguments and numerical convention as the primal call.

## Special functions

| Primal | JVP or gradient | VJP |
| --- | --- | --- |
| `bessel_in` | `bessel_in_jvp` | scalar output uses the same derivative times the cotangent |
| `bessel_in_array` | `bessel_in_array_jvp` | `bessel_in_array_vjp` |
| `bessel_kn` | `bessel_kn_jvp` | scalar contraction |
| `dawson` | `dawson_jvp`, `dawson_grad` | scalar contraction |
| `gamma_lower` | `gamma_lower_jvp`, `gamma_lower_jvp_da` | form from requested active parameters |
| `gamma_reg_p` | `gamma_reg_p_jvp`, `gamma_reg_p_grad` | scalar contraction |
| complex Bessel functions | corresponding complex JVP | contract manually until a public VJP is available |
| `hyperg_1f1_a1` | `hyperg_1f1_a1_jvp` | `hyperg_1f1_a1_vjp` |
| `fortnum_erf`, `fortnum_erfc` | corresponding JVP or gradient | scalar contraction |

Special-function products retain stable primal regime selection and neighboring
function identities. Match scaled and unscaled conventions.

## FFT and quadrature

`fft_c2c_jvp` and `fft_r2c_jvp` apply the primal linear operator to a tangent.
`fft_c2c_vjp` and `fft_r2c_vjp` apply the documented adjoint transform.

`gauss_legendre_jvp` and `gauss_legendre_vjp` contract a fixed weighted sum.

For a parameterized integral:

- use `integrate_fixed_jvp` for inactive bounds
- add `integrate_moving_lower_jvp` or `integrate_moving_upper_jvp` for active
  bounds
- use `integrate_qag*_jvp` to differentiate a recorded adaptive trace

The mathematical-integral and frozen-trace products are different contracts.

## Roots and solves

Use analytical implicit products:

| State | JVP | VJP |
| --- | --- | --- |
| scalar root | `root_implicit_jvp` | `root_implicit_vjp` |
| vector root | `multiroot_implicit_jvp` | `multiroot_implicit_vjp` |
| fixed point | `fixed_point_jvp` | `fixed_point_vjp` |
| linear solve | `linear_solve_jvp` | `linear_solve_vjp` |
| fitted spline coefficients | `bspline_fit_jvp_factored` | `bspline_fit_vjp_factored` |

Return and reuse converged Jacobians or factorizations. Do not differentiate
the nonlinear or linear iterations by default.

## ODEs

Run the primal solver first to record its accepted trace.

| Need | Interface |
| --- | --- |
| forward state sensitivity | `ode_integrate_jvp` |
| initial-state cotangent | `ode_integrate_vjp` |
| initial-state and parameter cotangents | `ode_integrate_parameter_vjp` |
| bounded reverse storage | checkpointed or recomputed VJP interface |
| implicit-stage products | `ode_implicit_stage_jvp`, `ode_implicit_stage_vjp` |
| event time and state | `ode_event_time_jvp`, `ode_event_state_jvp` |

The caller supplies RHS JVP or transpose actions. The discrete adjoint is the
transpose of the frozen Cash-Karp map. Verify whether the downstream objective
requires this discrete contract or a continuous sensitivity.

## Interpolation

Lagrange and B-spline products expose separate activity for evaluation points,
sampled values or coefficients, support nodes or breakpoints, and combined
arguments.

Call `grid_search_derivative_status` or
`bspline_span_derivative_status` when a direction can cross a cell or knot
boundary. A fixed-cell product is not valid across a changed index.

## Random estimators

`fortnum_rng` has no derivative with respect to seeds, counters, or draws.
Differentiate estimator parameters at the caller layer through a
reparameterized sample, score estimator, or analytical expectation.

## Downstream acceptance

For every migrated product:

1. state active and inactive arguments
2. match the primal convention
3. test a directional finite difference or exact sensitivity
4. test the adjoint identity for JVP/VJP pairs
5. test residual linearization for implicit products
6. benchmark complete value-plus-product wall clock and peak memory
7. measure forward/reverse scaling for the downstream input and output counts
