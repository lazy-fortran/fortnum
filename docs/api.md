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

`fortnum_special` re-exports the common special-function surface:

| Family | Primal | Products |
| --- | --- | --- |
| modified Bessel | `bessel_in`, `bessel_in_array`, `bessel_kn` | `bessel_in_jvp`, `bessel_kn_jvp` |
| Dawson | `dawson` | `dawson_jvp`, `dawson_grad` |
| incomplete gamma | `gamma_lower`, `gamma_reg_p` | `gamma_lower_jvp` |
| confluent hypergeometric | `hyperg_1f1`, `hyperg_1f1_a1` | `hyperg_1f1_a1_jvp`, `hyperg_1f1_a1_vjp` |
| Jacobi/simplex polynomials | `jacobi_p`, `scaled_jacobi_p`, `triangle_dubiner`, `tetrahedron_koornwinder` | `jacobi_p_derivative` |
| Ferrers associated Legendre | `legendre_p` | `legendre_p_derivative` |
| ordinary Legendre second kind | `legendre_q` | `legendre_q_derivative` |
| complex spherical harmonics | `spherical_harmonic` | angular derivatives, `spherical_harmonic_product_coefficient` |
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
  order on \([-1,1]\), with the Condon-Shortley phase, and real ordinary
  \(Q_\ell(x)\) on the \(x>1\) branch
- `fortnum_special_spherical`: standard orthonormal complex \(Y_\ell^m\) on
  \(0\le\theta\le\pi\), with analytical theta and phi derivatives
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

`legendre_q(l,x)` is the real ordinary Legendre function of the second kind
on \(x>1\), with
\[
Q_0(x)=\tfrac12\log\frac{x+1}{x-1},\qquad
(l+1)Q_{l+1}(x)=(2l+1)xQ_l(x)-lQ_{l-1}(x).
\]
`legendre_q_derivative` uses the corresponding DLMF derivative recurrence.
Invalid degrees or \(x\le1\) return NaN. The present upward recurrence is
intended for moderate degrees; a uniform large-degree continuation remains a
separate roadmap item.

`spherical_harmonic(l,m,theta,phi)` uses the Condon-Shortley phase and the
orthonormal convention of DLMF 14.30. The azimuth is periodic; the polar
angle must satisfy \(0\le\theta\le\pi\). The angular derivative entry points
are analytical and are intended away from the coordinate poles.

`spherical_harmonic_product_coefficient(l1,m1,l2,m2,L,M)` returns the real
Gaunt coefficient (G_{l_1m_1,l_2m_2}^{LM}) in

\[
Y_{l_1}^{m_1}Y_{l_2}^{m_2}
 = \sum_{L,M}G_{l_1m_1,l_2m_2}^{LM}Y_L^M.
\]

The implementation uses the two Wigner-3j symbols and the finite Racah sum
from [DLMF 34.3](https://dlmf.nist.gov/34.3). Triangle, parity, and azimuthal
selection rules return an exact zero; invalid degree/order indices return NaN.
Logarithmic factorial scaling is intended for moderate degrees and avoids the
pointwise Legendre evaluator, so products can be checked against an
independent angular-quadrature oracle.

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

The zero-order \(Q\) branch uses the DLMF zero-balanced continuation in
\(1-x^{-2}\) near \(x=1\), and associated orders are generated from its
analytical derivative and the generated order recurrence. This avoids the raw
hypergeometric series' near-cut nonconvergence. For half-integer degrees at
ordinary torus aspect ratios, (P) is continued upward in degree and the
recessive (Q) branch is continued with a scaled Miller backward recurrence;
the degree-80/order-4 values are checked against an independent 50-digit
reference and the three-term recurrence. Very close to the cut, the
zero-balanced direct branch remains the preferred path; a uniform asymptotic
envelope for arbitrarily large degree is still outside this API. The
continued-fraction and uniform-asymptotic literature is
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
- `fft_c2c_plan_init`
- `fft_c2c` (optionally with a caller-owned complex-transform plan)
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

## Structured tensor products

`fortnum_tensor_product` provides `tensor_product_operator_t` for a caller-owned
array of square `tensor_factor_t` matrices. Factor 1 is the fastest-varying
dimension, so the represented matrix is (A_d\otimes\cdots\otimes A_1).
`matvec`, `matmat`, and `diagonal` apply the factors without assembling the
full Kronecker matrix. Invalid factor shapes and vector or RHS dimensions
return `FORTNUM_DOMAIN_ERROR` through `fortnum_status_t`.

For OpenACC, `enter_data(status, n_rhs)` copies the factors and allocates
persistent vector and optional multi-RHS workspaces; `exit_data(status)` ends
that lifetime. `matvec_device` and `matmat_device` then apply the contractions
inside a caller-owned OpenACC data region whose input and output arrays are
present. This keeps factor and work arrays resident across repeated products;
the independent device test checks the results against the same dense oracle as
the host path. OpenMP-target and CUDA Fortran backends remain separate roadmap
items.

`fortnum_toeplitz` provides `toeplitz_operator_t` for one-dimensional grid
operators. Its first column and optional first row define the Toeplitz matrix;
omitting the row gives the symmetric covariance case. Initialization caches a
circulant embedding spectrum, and `matvec`/`matmat` use FFT products without
forming the dense matrix. The FFT convention is unnormalized in both
directions, so the operator applies the required inverse-transform scaling
internally. The current implementation is host-resident; accelerator-resident
FFT products remain a separate integration item.

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

`fortnum_ode_geometric` provides structure-preserving geometric time
integrators and their analytical products for caller-supplied Hamiltonian
systems.

Method modules expose:

| Module | Surface |
| --- | --- |
| `fortnum_ode_cash_karp` | one Cash-Karp step |
| `fortnum_ode_dop853` | DOP853 step and drivers |
| `fortnum_ode_extrapolation` | adaptive Gragg-Bulirsch-Stoer extrapolation integration |
| `fortnum_ode_gauss_radau` | adaptive 15th-order Gauss-Radau integration and drivers |
| `fortnum_ode_tdrk` | two-derivative Runge-Kutta and general Runge-Kutta-Nystrom integration |
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

`fortnum_krylov` provides matrix-free Krylov products. The caller supplies a
stateless matrix-vector procedure; `real_conjugate_gradient_operator` solves
real symmetric positive-definite systems with an optional preconditioner and
reports iterations plus the true final residual. The
`real_conjugate_gradient_matmat_operator` variant applies independent CG
recurrences to multiple right-hand sides while batching active directions in
one matrix-matrix callback. `complex_gmres_operator`
provides restarted complex GMRES using reorthogonalized modified
Gram--Schmidt and complex Givens rotations. Structured kernel operators and
device-resident callbacks are planned consumers of these contracts.

`fortnum_cholesky` provides `cholesky_factorization_t` and the corresponding
`cholesky_factorize`, `cholesky_solve_vector`, `cholesky_solve_matrix`,
`cholesky_solve_lower_matrix`, and `cholesky_log_determinant` procedures.
Factorization and solve failures are returned through `fortnum_status_t`; the
lower-triangular solve is exposed separately for covariance and posterior
calculations that need `L^-1 b` without a second triangular solve.

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

`fortnum_sobol` provides caller-owned low-discrepancy state through `sobol_t`.
Use `sobol_initialize`, `sobol_next`, `sobol_skip`, and `sobol_fill` to emit
reproducible points in the unit cube. The implementation supports dimensions
through `SOBOL_MAX_DIMENSION`; `SOBOL_TABULATED_DIMENSION` identifies the
dimensions using the published Joe--Kuo initial direction values.

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

## Choosing the differentiation backend

`fortnum_ad_backend` selects which engine's kernels the library calls.

fortnum carries two sets of generated derivatives for the same operators.
`FORTNUM_AD_FORTSYM` is the symbolic set: for a small closed-form operator it
can beat anything a differentiator produces, because it simplifies the
expression rather than the program. `FORTNUM_AD_FORTAD` is the
source-transformation set, which scales to operators with loops and branches
where a symbolic form would blow up, and covers every operator the Enzyme
oracle covers without needing an external toolchain.

Enzyme is not one of the choices. It is the independent second answer the
fixtures compare against, and it never appears in a library build - see
`design/enzyme_toolchain.md`. An operator the oracle cannot differentiate is a
limit on what can be cross-checked, not on what fortnum computes.

`FORTNUM_AD_ENGINE` names the one in use and defaults to `FORTNUM_AD_FORTAD`.
It is a named constant, so the compiler folds the branch and the unused call
costs nothing at runtime. Both sets stay compiled and both stay tested against
each other. An operator with no fortad kernel falls back to fortsym on its own,
so the choice is per operator rather than per library.
