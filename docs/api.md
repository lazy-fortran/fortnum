# API guide

This guide groups the supported public interfaces. Production module sources
under `src/` define exact argument order, optional arguments, shapes, and
status behavior.

All examples use:

```fortran
use, intrinsic :: iso_fortran_env, only: dp => real64
```

## Kinds and status

`fortnum_kinds` exports `dp`, `sp`, `i4`, and `i8`.

`fortnum_status` exports:

```fortran
type(fortnum_status_t) :: status

logical = status_ok(status)
call status_set(status, code, message)
```

Status codes are:

| Code | Meaning |
| --- | --- |
| `FORTNUM_OK` | Successful result |
| `FORTNUM_DOMAIN_ERROR` | Invalid input or undefined operation |
| `FORTNUM_CONVERGENCE_ERROR` | Iteration stopped without convergence |
| `FORTNUM_NOT_IMPLEMENTED` | Requested supported surface is not implemented |

Numerical routines return status explicitly when failure is part of the
contract. They do not use mutable global error state.

## Special functions

`fortnum_special` re-exports the common real-valued surface:

| Family | Primal | Products |
| --- | --- | --- |
| modified Bessel | `bessel_in`, `bessel_in_array`, `bessel_kn` | `bessel_in_jvp`, `bessel_kn_jvp` |
| Dawson | `dawson` | `dawson_jvp`, `dawson_grad` |
| incomplete gamma | `gamma_lower`, `gamma_reg_p` | `gamma_lower_jvp` |
| confluent hypergeometric | `hyperg_1f1`, `hyperg_1f1_a1` | `hyperg_1f1_a1_jvp`, `hyperg_1f1_a1_vjp` |
| Jacobi/simplex polynomials | `jacobi_p`, `scaled_jacobi_p`, `triangle_dubiner`, `tetrahedron_koornwinder` | `jacobi_p_derivative` |
| Ferrers associated Legendre | `legendre_p` | `legendre_p_derivative` |
| toroidal associated Legendre | `toroidal_p`, `toroidal_q` | `toroidal_p_derivative`, `toroidal_q_derivative` |

Domain modules expose additional products:

- `fortnum_special_bessel`: array JVP and VJP products
- `fortnum_special_complex_bessel`: complex `J`, `I`, and `K`, including
  scaled variants and JVPs
- `fortnum_special_dawson`: primal, generated outer product, JVP, and gradient
- `fortnum_special_gamma`: argument and parameter JVPs plus gradients
- `fortnum_special_hypergeometric_1f1`: `1F1`, specialized `a=1`, and
  `1F1M`
- `fortnum_special_jacobi`: Jacobi \(P_n^{(\alpha,\beta)}(x)\), its
  derivative, a homogeneous scaled form with removable collapsed-coordinate
  limits, and orthogonal Dubiner/Koornwinder modes on reference simplices
- `fortnum_special_legendre`: Ferrers \(P_\ell^m(x)\) for integer degree and
  order on \([-1,1]\), with the Condon-Shortley phase
- `fortnum_special_toroidal`: Hobson \(P_{n-1/2}^m(x)\) and
  \(Q_{n-1/2}^m(x)\), plus \(x\)-derivatives, for nonnegative integer
  \(n,m\) and \(x>1\)
- `fortnum_special_erf_cbind`: C-interoperable `erf`, `erfc`, JVPs, and
  gradients

Example:

```fortran
use fortnum_special, only: bessel_in, bessel_in_jvp

real(dp) :: value, tangent
value = bessel_in(2, 0.75_dp)
call bessel_in_jvp(2, 0.75_dp, 1.0_dp, tangent)
```

### Legendre and toroidal normalization

`legendre_p(l,m,x)` is the Ferrers function on the cut. Negative orders use
\[
P_l^{-m}(x)=(-1)^m\frac{(l-m)!}{(l+m)!}P_l^m(x).
\]
Values outside \([-1,1]\) are NaN. The derivative entry point is defined on
the open interval because endpoint derivatives may be singular.

`toroidal_p(n,m,x)` and `toroidal_q(n,m,x)` use degree \(n-\tfrac12\), not
\(n+\tfrac12\), and return Hobson-normalized functions. They are therefore
directly compatible with the conventional toroidal harmonics used after
separating Laplace's equation in toroidal coordinates. Invalid degree, order,
or \(x\) is reported as NaN.

The implementation follows
[DLMF 14.3](https://dlmf.nist.gov/14.3) for the hypergeometric definitions,
[DLMF 14.10](https://dlmf.nist.gov/14.10) for recurrence and derivative
relations, and [DLMF 14.19](https://dlmf.nist.gov/14.19) for the toroidal
specialization. In particular, DLMF's Olver-normalized
\(\boldsymbol{Q}_\nu^\mu\) is converted to Hobson \(Q_\nu^\mu\); the two must
not be interchanged.

The present Gauss-series implementation targets moderate nonnegative degree
and order at ordinary torus aspect ratios. It does not claim the near-\(x=1\),
large-order envelope of the continued-fraction and uniform-asymptotic
algorithm in
[Gil, Segura, and Temme (2000)](https://ir.cwi.nl/pub/1181/1181D.pdf).
The recurrence coefficients and series-term update are emitted by fortsym;
the exact generator revisions and regeneration commands are recorded in the
generated source banners.

## Quadrature and integration

`fortnum_quadrature` provides fixed rules:

- `gauss_legendre`
- `gauss_legendre_ab`
- `gauss_gen_laguerre`
- `gauss_legendre_jvp`
- `gauss_legendre_vjp`
- `gauss_legendre_grad`

`fortnum_integrate_gk` provides one Gauss-Kronrod panel through `gk_apply` and
a finite-interval driver through `integrate_gk`.

`fortnum_integrate` provides caller-owned adaptive state:

```fortran
type(integrate_workspace_t) :: workspace
type(integrate_epstab_t) :: epsilon_table
type(integrate_result_t) :: result
```

Primal drivers:

| Routine | Domain |
| --- | --- |
| `integrate_qag` | finite interval, selectable Gauss-Kronrod rule |
| `integrate_qags` | finite interval with epsilon extrapolation |
| `integrate_qagp` | finite interval with caller breakpoints |
| `integrate_qagiu` | semi-infinite or doubly infinite interval |
| `integrate` | allocating convenience wrapper around QAG |
| `integrate_cquad` | CQUAD-style adaptive integration in `fortnum_cquad` |
| `levin_u_accel` | Levin u-transform in `fortnum_levin` |

The integrand interface accepts an optional caller context:

```fortran
function f(x, ctx) result(value)
    real(dp), intent(in) :: x
    class(*), intent(in), optional :: ctx
    real(dp) :: value
end function f
```

Analytical products:

- `integrate_fixed_jvp`
- `integrate_moving_lower_jvp`
- `integrate_moving_upper_jvp`
- `integrate_qag_jvp`
- `integrate_qags_jvp`
- `integrate_qagp_jvp`
- `integrate_qagiu_jvp`

The first three differentiate the mathematical integral. The adaptive products
replay the accepted subdivision stored in `integrate_result_t`.

## FFT

`fortnum_fft` exports:

- `fortnum_fft_plan_t`
- `fft_plan_init`
- `fft_c2c`
- `fft_r2c`
- `fft_c2c_jvp`, `fft_c2c_vjp`
- `fft_r2c_jvp`, `fft_r2c_vjp`

`sign=-1` applies the unnormalized forward transform and `sign=+1` applies the
unnormalized inverse transform. The JVP uses the same sign. The VJP uses the
opposite sign with no extra `1/n` factor. Applying forward then inverse and
dividing by `n` recovers the input. Independent direct-DFT tests enforce these
sign and normalization conventions for the primal, JVP, and VJP.

Length-eight complex transforms dispatch to a fixed production execution leaf
shared by the analytical product and its autodiff benchmark candidate.
That internal leaf lives in `fortnum_fft8_kernel`; callers should continue
through `fortnum_fft`.

## ODEs

The main stateful API is in `fortnum_ode`:

```fortran
type(ode_problem_t) :: problem
type(ode_workspace_t) :: workspace
type(ode_solution_t) :: solution

problem%rhs => rhs
problem%t0 = 0.0_dp
problem%t1 = 1.0_dp
problem%y0 = [1.0_dp]
call ode_integrate(problem, workspace, solution, status)
```

`ode_problem_t` owns tolerances, step limits, event configuration, and callback
pointers. `ode_solution_t` owns the accepted trace and optional terminal event.
`ode_solve` is the allocating convenience wrapper. `fortnum_ode_wrapper`
provides `ode_at` for requested output times.

Method modules expose:

| Module | Surface |
| --- | --- |
| `fortnum_ode_cash_karp` | one Cash-Karp step |
| `fortnum_ode_dop853` | DOP853 step and drivers |
| `fortnum_ode_rk8pd` | stateful RK8PD evolution |
| `fortnum_ode_ddeabm` | stateful Adams method |
| `fortnum_ode_vode` | stateful VODE-style integration |
| `fortnum_ode_events` | event scan and event derivatives |

First-order products in `fortnum_ode`:

- `ode_integrate_jvp`
- `ode_integrate_vjp`
- `ode_integrate_parameter_vjp`
- `ode_integrate_vjp_checkpointed`
- `ode_integrate_vjp_recomputed`
- `ode_implicit_stage_jvp`
- `ode_implicit_stage_vjp`

The adaptive trace is inactive. Callers provide RHS JVP or transpose actions.
See [design/ode.md](design/ode.md) for continuous and discrete semantics.

## Roots

`fortnum_roots` provides:

- `root_bisect`
- `root_newton`
- `root_brent`
- `root_jvp`, `root_vjp`, `root_grad`
- `root_implicit_jvp`, `root_implicit_vjp`

The generic implicit products accept caller callbacks for residual products.
They differentiate the residual equation at the converged root.

`fortnum_multiroot` provides:

- `multiroot_hybrid`, `multiroot_hybrids`
- `multiroot_jvp`, `multiroot_vjp`, `multiroot_grad`
- factored JVP and VJP variants
- generic implicit JVP and VJP boundaries
- derivative and sorting helpers used by compatibility paths

`fortnum_fixed_point` differentiates a fixed-point residual through
`fixed_point_jvp` and `fixed_point_vjp`.

`fortnum_multiroot_rc` is a fixed-capacity reverse-communication solver.
`fortnum_roots_complex` finds zeros inside a complex rectangle.

## Linear algebra

`fortnum_linalg` provides fixed-size primitives and a fixed-capacity LU object:

- `det2`, `det3` and JVP/VJP products
- `inv2`, `inv3` and fused value/JVP/VJP products
- `jacobian_ok3`
- `lu_factor`, `lu_solve_factored`, `lu_solve`
- `lu_factorization_t`
- `linear_solve_jvp`, `linear_solve_vjp`
- factored and multiple-right-hand-side JVP/VJP variants

Solve products use the implicit equations

\[
A\,dx = db-dA\,x,
\qquad
A^T\lambda=u.
\]

They never form an inverse for the derivative of a solve.

`fortnum_krylov` provides restarted matrix-free complex GMRES through
`complex_gmres_operator`. The caller supplies a stateless matrix-vector
procedure; the solver uses reorthogonalized modified Gram--Schmidt, complex
Givens rotations, and reports iterations plus the true final residual.

## Interpolation and splines

`fortnum_interp` exports grid search and a directional status check for cell
crossings.

`fortnum_polynomial` exports Lagrange weights and products for active:

- evaluation point
- sampled values
- support nodes
- evaluation point and sampled values together

`fortnum_bspline` exports:

- `bspline_workspace_t`
- initialization and knot setup
- basis and derivative evaluation
- knot-span lookup and crossing status
- JVP/VJP products for the evaluation point, coefficients, and breakpoints
- combined evaluation-point/coefficient products
- factorization-reusing implicit products for fitted coefficients

Products involving location arguments are defined on a fixed cell or knot
span. Use the status guards when a perturbation can cross a boundary.

## Random numbers

`fortnum_rng` uses caller-owned `rng_t` state:

```fortran
type(rng_t) :: generator
call rng_seed(generator, seed, status)
call rng_uniform(generator, value)
call rng_normal(generator, normal_value)
```

`rng_split` derives a child stream without advancing its parent.
`rng_next_u64` returns a raw word. `rng_threefry2x64` exposes the block
function for known-answer tests.

Seeds, keys, counters, and draws have no derivative products. Differentiate
distribution-level estimators in the caller when needed.

## Optimizer-facing interfaces

`fortnum_active_vector` maps named array blocks to one flat active vector:

- `fortnum_active_layout_t`
- `layout_init`, `layout_add`, `layout_index`, `layout_block`
- `pack_block`, `unpack_block`

`fortnum_ad_interfaces` defines backend-independent value, JVP, VJP, gradient,
and HVP callbacks plus derivative provenance and quality.

`fortnum_derivative_registry` selects a validated candidate before hot loops
using committed workload metadata.

## C interface

`include/fortnum.h` declares the installed C ABI for selected special
functions, quadrature, integration, roots, ODEs, and B-splines. Opaque handles
own state for B-spline and RK8PD interfaces. C callers must check integer status
codes and release every created handle.

## Testing helpers

`fortnum_oracle` reads CSV reference tables through `oracle_read` and validates
a callback through `oracle_check`. It is a test-support module, not a
production reference-data dependency.
