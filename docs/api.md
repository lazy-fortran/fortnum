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
| Faddeeva and plasma dispersion | `faddeeva_w`, `plasma_dispersion_z` | `plasma_dispersion_z_derivative` |
| incomplete gamma | `gamma_lower`, `gamma_reg_p` | `gamma_lower_jvp` |
| complex digamma | `digamma_complex` | `trigamma_complex` |
| confluent hypergeometric | `hyperg_1f1`, `hyperg_1f1_a1` | `hyperg_1f1_a1_jvp`, `hyperg_1f1_a1_vjp` |
| Jacobi/simplex polynomials | `jacobi_p`, `scaled_jacobi_p`, `triangle_dubiner`, `tetrahedron_koornwinder` | `jacobi_p_derivative` |
| Ferrers associated Legendre | `legendre_p` | `legendre_p_derivative` |
| normalized Legendre table | `legendre_p_normalized_table` | none |
| ordinary Legendre second kind | `legendre_q` | `legendre_q_derivative` |
| Riccati-Bessel tables | `riccati_bessel` | none |
| complex spherical harmonics | `spherical_harmonic` | angular derivatives, `spherical_harmonic_product_coefficient` |
| toroidal associated Legendre | `toroidal_p`, `toroidal_q` | `toroidal_p_derivative`, `toroidal_q_derivative` |

Domain modules expose additional products:

- `fortnum_special_bessel`: array JVP and VJP products
- `fortnum_special_complex_bessel`: complex `J`, `I`, and `K`, including
  scaled variants and JVPs
- `fortnum_special_dawson`: primal, generated outer product, JVP, and gradient
- `fortnum_special_faddeeva`: \(w(z)=e^{-z^2}\operatorname{erfc}(-iz)\) and
  the Fried-Conte \(Z(\zeta)=i\sqrt\pi\,w(\zeta)\) on the whole complex
  plane, with \(Z'=-2(1+\zeta Z)\)
- `fortnum_special_gamma`: argument and parameter JVPs plus gradients, and
  complex \(\psi(z)\), \(\psi'(z)\)
- `fortnum_special_hypergeometric_1f1`: `1F1`, specialized `a=1`, and
  `1F1M`
- `fortnum_special_jacobi`: Jacobi \(P_n^{(\alpha,\beta)}(x)\), its
  derivative, a homogeneous scaled form with removable collapsed-coordinate
  limits, and orthogonal Dubiner/Koornwinder modes on reference simplices
- `fortnum_special_legendre`: Ferrers \(P_\ell^m(x)\) for integer degree and
  order on \([-1,1]\), with the Condon-Shortley phase, and real ordinary
  \(Q_\ell(x)\) on the \(x>1\) branch
- `fortnum_special_riccati_bessel`: \(\hat j_l(x)=xj_l(x)\) and
  \(\hat y_l(x)=xy_l(x)\) for \(0\le l\le l_{\max}\), real \(x>0\)
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
Invalid degrees or \(x\le1\) return NaN. Because \(Q_l\) is the minimal
solution of the recurrence on \(x>1\), degrees with
\((l+1)\operatorname{arccosh}x>1\) use the backward ratio continued fraction
anchored at \(Q_0\); smaller degrees use the upward recurrence, whose
amplification is modest there. Close to the cut both directions are nearly
neutral and rounding grows like \(l^{3/2}\varepsilon\) (about \(3\times
10^{-12}\) relative at \(l=2000\)).

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

`legendre_p_normalized_table(lmax,mmax,x,p,condon_shortley)` fills the
caller array `p(0:lmax,0:mmax)` with fully normalized functions
\[
\bar P_l^m(x)=\sqrt{\tfrac{2l+1}{2}\tfrac{(l-m)!}{(l+m)!}}\,P_l^m(x),
\qquad \int_{-1}^{1}\bar P_l^m\bar P_k^m\,dx=\delta_{lk},
\]
with zeros for \(l<m\). Without the optional flag the Condon-Shortley phase
is omitted; with `condon_shortley=.true.`,
\(Y_l^m=\bar P_l^m(\cos\theta)e^{im\phi}/\sqrt{2\pi}\) equals
`spherical_harmonic`. The seed and degree recurrence follow Holmes and
Featherstone (2002) on values carried with a separate binary exponent, so
\(l_{\max}\) of several thousand neither overflows nor loses the
underflowing diagonal seed. The tested sum-rule residual is below
\(25\,l_{\max}\,\epsilon\) for \(l_{\max}=2000\). The routine allocates
nothing; \(|x|>1\) yields NaN.

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

### Plasma dispersion, digamma and Riccati-Bessel conventions

`faddeeva_w(z)` and `plasma_dispersion_z(zeta)` are elemental and entire.
\(Z(\zeta)=\pi^{-1/2}\int e^{-t^2}/(t-\zeta)\,dt\) for \(\operatorname{Im}\zeta>0\),
so the same call is the Landau continuation for \(\operatorname{Im}\zeta\le0\);
\(Z(0)=i\sqrt\pi\). The algorithm is Gautschi's, with the Poppe-Wijers
(TOMS 680) regions, at about \(10^{-14}\) relative accuracy. In the lower half
plane \(2e^{-z^2}\) overflows to IEEE infinity once
\(\operatorname{Im}(z)^2-\operatorname{Re}(z)^2\gtrsim709\).

`digamma_complex(z)` and `trigamma_complex(z)` are elemental on the whole
plane: reflection for \(\operatorname{Re}z<1/2\), upward shift to \(|z|\ge16\)
and the Bernoulli series through \(B_{14}\). Poles \(z=0,-1,\dots\) return
\(+\infty\). \(\operatorname{Re}\psi(1+i\eta)\) is
`real(digamma_complex(cmplx(1, eta, dp)))`.

`call riccati_bessel(lmax, x, jhat, yhat)` fills `jhat(0:lmax)` and
`yhat(0:lmax)` without allocation, with \(\hat j_0=\sin x\) and
\(\hat y_0=-\cos x\). `yhat` uses the stable upward recurrence; for
\(l_{\max}>x\), `jhat` uses Miller's backward recurrence normalized by the
Wronskian \(\hat j_1\hat y_0-\hat j_0\hat y_1=1\). Invalid arguments fill both
arrays with NaN.

## Quadrature and integration

`fortnum_quadrature` provides fixed rules:

- `gauss_legendre`
- `gauss_legendre_ab`
- `gauss_lobatto_legendre`
- `gauss_gen_laguerre`
- `gauss_legendre_jvp`
- `gauss_legendre_vjp`
- `gauss_legendre_grad`

`gauss_lobatto_legendre(n,x,w)` returns the \(n\ge2\) point
Gauss-Lobatto-Legendre rule on \([-1,1]\): ascending nodes including
\(\pm1\), the zeros of \(P_{n-1}'\) inside, and weights
\(2/(n(n-1)P_{n-1}(x_i)^2)\). It is exact through degree \(2n-3\); nodes
are exactly antisymmetric and agree with a quad-precision reference to
\(4\epsilon\) up to \(n=201\), and weights to \(8n\epsilon\) relative.
The linear products of `gauss_legendre` apply to these weights.

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
| `fortnum_ode_rk54_device` | allocation-free four-state Cash-Karp and Dormand-Prince reverse communication for GPU callers |
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

`fortnum_cholesky` provides a reusable positive-definite factorization with
vector and matrix solves and a log-determinant operation.

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

`fortnum_symmetric_eigen` provides real symmetric eigensolvers.
`symmetric_eigen` returns ascending eigenvalues and orthonormal eigenvectors of
a dense matrix by Householder tridiagonalization and implicit QL.
`block_lanczos_lowest` returns the lowest eigenpairs of a matrix-free symmetric
operator (the `fortnum_krylov` callback contract) by block Lanczos with full
reorthogonalization and Rayleigh--Ritz; the block size bounds the eigenvalue
multiplicity it resolves. Reported residual norms are recomputed from operator
images, so each Ritz value lies within its residual of an eigenvalue; that the
values are the lowest is a candidate claim without a separate gap argument.

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

It also exports `barycentric_weights(n,xp,bw)`, for the second-kind
barycentric interpolation formula (weights scaled to \(\max|bw|=1\)), and
`lagrange_differentiation_matrix(n,xp,d)`, the nodal matrix
\(d_{ij}=L_j'(x_i)\) with negative-row-sum diagonal. With
`gauss_lobatto_legendre` nodes these give the spectral-element/DVR basis.

`fortnum_bspline` exports:

- `bspline_workspace_t`
- initialization and knot setup
- basis and derivative evaluation
- knot-span lookup and crossing status
- JVP/VJP products for the evaluation point, coefficients, and breakpoints
- combined evaluation-point/coefficient products
- factorization-reusing implicit products for fitted coefficients

`fortnum_bspline_lsq` exports `bspline_1d_lsq_cgls`. It fits sampled values to
the existing `bspline_workspace_t` with matrix-free CGLS, returning a status,
iteration count, and residual norm. The basis is cached once per fit; normal
equations are not formed. Tensor-product fitting and direct interpolation are
not part of this API.

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

## Verified computing

The verified modules return enclosures: every result contains the exact real
or complex result for all inputs in the argument enclosures. They use only
correctly rounded IEEE operations with outward rounding by the successor
formula, so no rounding-mode change and no libm accuracy claim enters. The
[verified design](design/verified.md) states the semantics, trust
assumptions, and error models.

| Module | Surface |
| --- | --- |
| `fortnum_rounding` | `round_up`, `round_down`, directed `add_up`/`mul_up`/`div_up`/`sqrt_up` and their downward forms, `sum_up`, `sum_down`, `gamma_up` |
| `fortnum_interval` | `interval_t`, `cinterval_t`, arithmetic operators, `sqrt`, `exp`, `log`, `sin`, `cos`, `sincos`, `sinh`, `cosh`, `abs`, `sqr`, `linear_interpolation_remainder`, `interval_pi`, `hull`, `intersect`, `mid`, `rad`, `width`, `mag`, `mig`, `contains`, `subset`, `disjoint`, `cabs_up`, `cabs_down` |
| `fortnum_interval_qp` | binary128 `qinterval_t` with `qrat`, `qsqrt`, `qinterval_pi`, `to_interval` |
| `fortnum_ball` | complex balls `ball_t`: operators, `binv`, `bsqrt`, `bpowi`, modulus and component bounds, box conversions |
| `fortnum_idual` | `idual_t` interval forward AD with up to `idual_max_vars` seeded variables |
| `fortnum_cdual` | `cdual_t` complex-box forward AD with up to `cdual_max_vars` seeded complex variables, real-interval scaling, `cdual_sincos` |
| `fortnum_taylor_series` | order-by-order Taylor arithmetic over `interval_t`, `idual_t`, and `cdual_t` coefficient arrays: `ts_mul`, `ts_div`, `ts_sincos`, `ts_exp`, `ts_sqrt`, `ts_lin`, `ts_zero` |
| `fortnum_bernstein` | rigorous Bernstein-form evaluation on real (`bernstein_eval_real_box`) and complex (`bernstein_eval_complex_box`) boxes by local Taylor re-expansion, `bernstein_hull` convex-hull value/derivative bounds, `bernstein_grid_t` cached-centre evaluation, `bernstein_lsq_fit` |
| `fortnum_validated_ode` | validated ODE integration with an abstract right-hand side `ode_rhs_t`: `picard_apriori`, `gronwall_bound`, `eps_ball_matrix`, `lohner_step_jacobian` (first-order Q = I + h A (1 + eps)), `taylor_lohner_predictor` (order-q centre predictor; optional `ay`/`q_jac` arguments give the order-q variational step Jacobian Q = sum_k D y_k(Y) h^k + D y_{q+1}(Y) M h^(q+1)), `lohner_qr_frame`, `lohner_inverse_enclosure`, `lohner_state_t`/`lohner_step`/`lohner_integrate` (optional `variational` flag, default true, selects the order-q Jacobian over the first-order fallback; optional `apriori_y`/`apriori_f`/`ybx_out` pass-through of the step's a priori box, right-hand side there, and box Taylor coefficients), `event_crossing_newton`, `section_fn_if`/`section_crossing` (validated first-return map to a section g(y) = 0 with a prescribed crossing direction, by adaptive-step Lohner propagation, a node-based sign-change state machine, and interval-Newton crossing localization on the detected step's box-valid local Taylor polynomial), `stop_unresolved`/`stop_crossed`/`stop_avoided`/`stopping_enclosure` (enclosure of m_T = max_{t<=T} g(y(t)) for a box of initial conditions: a lower bound from each step's refined grid-time box, an upper bound from each step's a priori box, classified crossed/avoided/unresolved) |
| `fortnum_cfft_rigorous` | `rigorous_fft_plan_t` with certified twiddles, `rigorous_fft_apply`, `rigorous_fft2_apply`, `conv_error_bound`, `conv_error_bound_2d`, `conv_nonneg`, `conv_nonneg_2d` |
| `fortnum_fseries` | ball Fourier series with l1 tails: `fseries_t` (1D) and `fseries2d_t` (2D, sigma-weighted) with sums, direct and FFT products, derivatives, weighted inner products, and Wiener reciprocals `fs_inverse`, `fs2_inverse` |
| `fortnum_verified_linalg` | matrix balls `mball_t`, `fro_up`, `norm1_up`, `norminf_up`, `norm2_up`, `norm2_tight_up`, `matmul_err`, `approx_inverse`, `verified_inverse`, Cholesky-certified `sym_lower_bound`, `eig_lower_bound`, `eig_upper_bound` (real symmetric and Hermitian), `interval_matmul`, `interval_matvec` |
| `fortnum_ode_residual_certificate` | a posteriori trajectory certificate for stored ODE candidates: `hermite_cubic_enclosure`, `linear_residual_bound` (interval action callback `residual_action_if`, optional cubic forcing), `nonlinear_residual_bound` (box callback `residual_rhs_box_if`, first order), `residual_radius_step` (log-norm stability model), `certify_linear_trajectory`, `certify_nonlinear_trajectory` |
| `fortnum_verified_matvec` | immutable binary64 matrix acting on interval vectors with a proved aggregate rounding bound: `matvec_bound_t`, `prepare_matvec_bound`, `enclose_matvec`, `matvec_row_bound` |
| `fortnum_sampled_norm` | sampled vector-norm upper bound valid with a declared probability law: `sampled_norm_upper` |
| `fortnum_verified_quadrature` | `composite_midpoint(callback, edges, integral, ok)` and `composite_midpoint_batch` with cellwise rigorous second-derivative bounds |

Verified midpoint quadrature accepts strictly increasing finite binary64 edges
as exact real endpoints. The callback receives the whole cell and an outward
enclosure of its mathematical midpoint, then returns an interval integrand value,
a nonnegative interval bounding the magnitude of its second derivative throughout
that cell, and success status. Each cell must have an absolutely continuous first
derivative and an essentially bounded second derivative; place nonsmooth knots
among the edges. The integral packet includes midpoint remainders and evaluation,
geometry, weighting and accumulation rounding. Nonuniform cells are supported.
The batch companion accepts an output interval array and a callback returning
matching value/curvature arrays, sharing integrand work across components.
Its scratch arrays exist once per call; neither interface allocates within
the cell loop or chooses a hidden tolerance/refinement schedule. Invalid inputs,
callback failure or nonfinite results return `ok=.false.` and an empty packet.
The caller owns the truth of its value/derivative bounds, refinement, improper
tails and root-bracket schedules. See [the verified design](design/verified.md).

The ODE residual certificate checks an untrusted stored trace: strictly
increasing binary64 nodes, values and slopes define an exact C1 cubic Hermite
reconstruction. Each cell's residual `p'/h - f(p)` is bounded over the whole
cell by interval Bernstein controls, and the Euclidean error radius is
accumulated with the caller's bound `mu` on the logarithmic 2-norm
(`mu <= 0`: amplification 1). Complex states use the real/imaginary split.
Callers supply scratch; the cell loop does not allocate. Failures return a
nonzero `fortnum_status_t` and infinite radii. See
[the verified design](design/verified.md).

The interval and ball modules also export the runtime interface that
`fortsym`-emitted rigorous kernels call: `ipoint`, `ienclose`, `iadd`,
`isub`, `imul`, `idiv`, `ineg`, `iinv`, `isqrt`, `ipowi`, `iscale` for
`interval_t`, and `bpoint`, `bcpoint`, `benclose`, `badd`, `bsub`, `bmul`,
`bdiv`, `bneg`, `binv`, `bsqrt`, `bpowi`, `bscale` for `ball_t`.

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
