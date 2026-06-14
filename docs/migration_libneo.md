# Migration: libneo GSL/FGSL and math-kit to fortnum

This document maps the math routines used in libneo's GSL/FGSL calls and
math-kit to their fortnum equivalents. For each mapping it states whether
derivative products are available now, planned, or excluded.

The fortnum derivative policy vocabulary is defined in `docs/design/ad.md`.
No derivative products ship yet; all are reserved for issue #40.

---

## Special functions

### Modified Bessel functions

| libneo / GSL | fortnum |
|--------------|---------|
| `gsl_sf_bessel_In(n, x)` | `bessel_in(n, x)` |
| `gsl_sf_bessel_In_array(nmin, nmax, x, result_array)` | `bessel_in_array(nmax, x, values(0:nmax))` |
| `gsl_sf_bessel_Kn(n, x)` | `bessel_kn(n, x)` |

Interface differences:
- `bessel_in` is `elemental`; pass a scalar or an array directly.
- `bessel_in_array` fills `values(0:nmax)`; it does not take `nmin`. For a
  subrange starting above zero, call `bessel_in_array` and index the result.
- `bessel_kn` is a pure function, not a subroutine with an error-handler
  argument. Domain errors (x <= 0) stop with `error stop` in the primal; a
  status path can be added later without changing the signature.

Derivative status (policy `analytic_rule`):
- Now: primal only.
- Issue #40: `bessel_in_grad`, `bessel_kn_grad` using the recurrences
  d/dx I_n = (I_{n-1} + I_{n+1}) / 2 and d/dx K_n = -(K_{n-1} + K_{n+1}) / 2.

### Dawson integral

| libneo / math-kit | fortnum |
|-------------------|---------|
| `dawson(x)` (math-kit `src/math/dawson.f90`) | `dawson(x)` |

The fortnum implementation is a clean-room port of the math-kit routine with
style conformance (`0.0d0` -> `_dp`). The algorithm is identical.

Derivative status (policy `analytic_rule`):
- Now: primal only.
- Issue #40: `dawson_grad` using F'(x) = 1 - 2x F(x).

### Incomplete gamma

| libneo / GSL | fortnum |
|--------------|---------|
| `gsl_sf_gamma_inc(a, x)` (lower incomplete, unnormalized) | `gamma_lower(a, x)` |
| `gsl_sf_gamma_inc_P(a, x)` (regularized P) | `gamma_reg_p(a, x)` |

Both routines are pure functions. Domain violations (`a <= 0` or `x < 0`) stop
with `error stop`; the primal signature carries no status argument by design.
A status path may be added later without changing the signature.

Derivative status (policy `analytic_rule`):
- Now: primal only.
- Issue #40: `gamma_lower_grad`, `gamma_reg_p_grad` using
  d/dx gamma(a,x) = x^{a-1} exp(-x).

### Complex-argument Bessel functions

| libneo / AMOS (via KiLCA/KAMEL) | fortnum |
|---------------------------------|---------|
| `zbesj` (J_n, complex z) | `bessel_j_complex(order, z, result, status)` |
| `zbesi` (I_n, complex z, KODE scaling) | `bessel_i_complex(order, z, scaled, result, status)` |
| `zbesi` sequence | `bessel_i_complex_array(order0, nseq, z, scaled, result, status)` |
| `zbesk` (K_n, complex z, KODE scaling) | `bessel_k_complex(order, z, scaled, result, status)` |
| `zbesk` sequence | `bessel_k_complex_array(order0, nseq, z, scaled, result, status)` |

Clean-room from DLMF chapter 10; no AMOS source is used. The `scaled` logical
replaces AMOS `KODE`: `.false.` is the unscaled value, `.true.` factors out
e^{Re z} for I and e^{z} for K. Negative integer order is handled internally.
Errors are a `fortnum_status_t` argument, not the AMOS `NZ`/`IERR` integers.

Derivative status (policy `analytic_rule`): forward complex products
`bessel_j_complex_jvp`, `bessel_i_complex_jvp`, `bessel_k_complex_jvp` ship now,
using the order recurrences (DLMF 10.6.1, 10.29.2, 10.29.4).

### Confluent hypergeometric 1F1 (Kummer M)

| libneo / consumer | fortnum |
|-------------------|---------|
| MEPHIT `hyper1F1` (Kummer M, complex args) | `hyperg_1f1(a, b, z, result, status)` |
| KAMEL KiLCA `hyper1F1` with a = 1 | `hyperg_1f1_a1(b, z, result, status)` |

Clean-room from DLMF chapter 13. Replaces the multi-variant selectors of MEPHIT
`src/hyper1F1.c` and KiLCA `math/hyper/hyper1F1.cpp` with one internal selector
(Kummer transformation, Taylor series, large-z asymptotics) exposing the single
converged value. `b` at or near a non-positive integer reports
`FORTNUM_DOMAIN_ERROR`.

Derivative status (policy `analytic_rule`): `hyperg_1f1_a1_jvp` /
`hyperg_1f1_a1_vjp` ship now, using d/dz M = (a/b) M(a+1,b+1,z) (DLMF 13.3.15).

### Error function (C ABI)

| KAMEL / GSL | fortnum |
|-------------|---------|
| `gsl_sf_erf(x)` | `fortnum_erf(x)` |
| `gsl_sf_erfc(x)` | `fortnum_erfc(x)` |

Fortran callers should use the F2008 intrinsics `erf`/`erfc` directly. The
`fortnum_erf`/`fortnum_erfc` wrappers exist to give the C ABI a callable symbol;
they forward to the intrinsics with no reimplementation.

Derivative status (policy `transparent`): `fortnum_erf_jvp`, `fortnum_erfc_jvp`,
`fortnum_erf_grad`, `fortnum_erfc_grad` ship now.

---

## Series acceleration

| libneo / GSL | fortnum |
|--------------|---------|
| `gsl_sum_levin_u_alloc` + `gsl_sum_levin_u_accel` | `levin_u_accel(terms, n, sum_accel, abserr, status)` |

Clean-room Levin-u transform (Levin 1973; Weniger 1989; Fessler-Ford-Smith
1983). The caller passes the recorded series terms `terms(1:n)`; the workspace
the GSL `gsl_sum_levin_u_workspace` held is internal, so no alloc/free pair is
needed. `sum_accel` and `abserr` replace the GSL out-parameters.

Derivative status (policy `primal_only`): none. The accelerated value is a
data-dependent nonlinear rational transform of the term sequence.

---

## Multidimensional root finding

| libneo / GSL | fortnum |
|--------------|---------|
| `gsl_multiroot_fdfsolver_hybridj` (analytic Jacobian) | `multiroot_hybrid(fdf, n, x0, x, status, ...)` |
| `gsl_multiroot_fdfsolver_hybridsj` (finite-diff Jacobian) | `multiroot_hybrids(fn, n, x0, x, status, ...)` |
| `gsl_deriv_central` | `deriv_central(f, x, h, result, abserr, status, ctx)` |
| index sort of a real array | `argsort(x, perm)` |

Interface differences:
- The residual and Jacobian come through abstract interfaces with an optional
  `ctx`, not a `gsl_multiroot_function_fdf` struct: `multiroot_fdf_t` returns F
  and the analytic Jacobian; `multiroot_fn_t` returns F only.
- No iterator object and no alloc/free. One call drives the iteration to
  convergence; `xtol`, `ftol`, `max_iter` are optional.
- Termination is `|F|_inf <= ftol` or `|dx|_inf <= xtol*(|x|_inf + xtol)`.
  A singular Jacobian reports `FORTNUM_DOMAIN_ERROR`.

Derivative status:
- Solvers: policy `implicit_rule`. Now primal only; reserved `multiroot_jvp`,
  `multiroot_vjp`, `multiroot_grad`.
- `deriv_central`, `argsort`: policy `primal_only`.

---

## Numerical integration

| libneo / FGSL / QUADPACK | fortnum |
|--------------------------|---------|
| `gsl_integration_qag` | `integrate_qag` |
| `gsl_integration_qags` | `integrate_qags` |
| `gsl_integration_qagp` | `integrate_qagp` |
| `gsl_integration_qagiu` (semi-infinite) | `integrate_qagiu` |
| Direct QUADPACK `dqag` / `dqags` | `integrate_qag` / `integrate_qags` |
| Single-rule QK15/21/31/61 panel | `gk_apply` |

Interface differences:
- The integrand interface requires an optional `ctx` argument:
  `function f(x, ctx) result(fx)`. A parameter-free integrand ignores `ctx`.
  This replaces GSL's `gsl_function` struct and the FGSL function pointer.
- Workspace and epsilon-table types are caller-owned and explicit.
  `integrate_workspace_t` replaces the GSL workspace handle;
  `integrate_epstab_t` is the Wynn epsilon table. Pass default-initialized
  values on the first call; reuse them across calls to avoid reallocation.
- Status is a `fortnum_status_t` argument, not an integer error code.
  `FORTNUM_OK` maps to GSL's `GSL_SUCCESS = 0`; `FORTNUM_CONVERGENCE_ERROR`
  maps to `GSL_EMAXITER` or `GSL_EROUND`; `FORTNUM_DOMAIN_ERROR` covers
  non-smooth derivative cases (break-point singularity at an endpoint, or
  bad integrand behaviour).
- `integrate` is a flat single-call wrapper that owns its workspace; it is
  the closest equivalent to a GSL integration call without explicit workspace
  management.

Derivative status (policy `trace_rule`):
- Now: primal only.
- Issue #40: `integrate_qag_jvp`, `integrate_qags_jvp`, `integrate_qag_vjp`,
  `integrate_qags_vjp`, `integrate_qag_grad`. Derivatives are taken at the
  frozen accepted subdivision recorded in `integrate_result_t`.

---

## Fixed-rule Gauss quadrature

| libneo / Burkardt | fortnum |
|-------------------|---------|
| `cgqf(order, kind=5, alpha, beta, a, b, x, w)` (gen_laguerre_rule.f90) | `gauss_gen_laguerre(n, alpha, x, w)` |

`transport.f90` builds the generalized-Laguerre rule for the `1/nu` transport
coefficient by calling `cgqf` with `kind = 5` (weight `x^alpha exp(-x)`),
`alpha = 5/2` for `D11` and `7/2` for `D12`. `gauss_gen_laguerre` replaces that
call and the bundled `gen_laguerre_rule.f90`.

Interface differences:
- `gauss_gen_laguerre` covers only `kind = 5` with the unit shift `a = 0`,
  scale `b = 1`. `beta` does not apply. For other Burkardt kinds or a shifted
  interval, use `gauss_legendre_ab` or open an issue.
- Nodes come back ascending, so `calc_D_one_over_nu` sums them in increasing
  speed order. The Burkardt routine returns the same set in the same order.
- No LAPACK: the rule is the eigendecomposition of the Jacobi matrix by the
  in-tree implicit-shift QL solver, matching the `gauss_legendre` Newton path.

Derivative status (policy `transparent`):
- The coefficient `calc_D_one_over_nu` takes the weights and abscissas as
  `intent(in)` and differentiates only the integrand `1/coll_freq_tot`. The
  map `f -> I = sum_i w_i f_i` is linear, so `gauss_legendre_jvp`,
  `gauss_legendre_vjp`, and `gauss_legendre_grad` apply to the returned weight
  vector unchanged. Nodes, weights, `n`, and `alpha` are inactive.

---

## FFT

| libneo / FFTW | fortnum |
|---------------|---------|
| `fftw_plan_dft_1d(..., FFTW_FORWARD)` + `fftw_execute` | `fft_plan_init` + `fft_c2c(z, -1)` |
| `fftw_plan_dft_1d(..., FFTW_BACKWARD)` + `fftw_execute` | `fft_plan_init` + `fft_c2c(z, +1)` |
| `fftw_plan_dft_r2c_1d` + `fftw_execute` | `fft_plan_init` + `fft_r2c(x, c, plan)` |

Interface differences:
- Plans are `fortnum_fft_plan_t` values on the stack, not pointers. No
  `fftw_destroy_plan` call is needed.
- The plan is optional in `fft_r2c` and `fft_c2c`; omit it to let the routine
  build a temporary plan. Supply it to reuse twiddle tables across transforms
  of the same length.
- The normalization convention matches FFTW: forward is unnormalized,
  inverse is unnormalized. Divide by `n` after the inverse when the
  application requires it.

Derivative status (policy `transparent`):
- Now: primal only.
- Issue #40: `fft_c2c_jvp`, `fft_c2c_vjp`, `fft_r2c_jvp`, `fft_r2c_vjp`.
  The DFT is linear; Enzyme differentiates the implementation directly.

---

## ODE integration

| libneo | fortnum |
|--------|---------|
| `odeint` (GSL `gsl_odeiv2_*` or custom) | `ode_integrate` / `ode_solve` |
| `gsl_odeiv_step_rk8pd` (RK8(7), tight tolerances) | `ode_integrate_dop` / `ode_solve_dop` |

Interface differences:
- The RHS is an abstract interface: `subroutine rhs(t, y, dydt, ctx)` with
  an optional unlimited polymorphic `ctx`. Thread an arbitrary parameter
  struct through `ctx` instead of a `void *` or a module variable.
- `ode_solve` is the flat call for callers that do not need the workspace or
  the recorded trace.
- `ode_at` evaluates the solution at a specified set of output times, matching
  the common pattern of requesting output at prescribed points.
- The default method is Cash-Karp RK5(4) with PI step control. If libneo used a
  different order or tableau, check that the tolerances give the needed
  accuracy; the step controller uses the same `rtol` and `atol` parameters.
- `ode_integrate_dop` / `ode_solve_dop` are the Prince-Dormand RK8(7)-13M
  (DOP853) drop-ins for the very tight tolerances where the higher order pays
  off. Same carriers and same `trace_rule` policy as the Cash-Karp path.

Derivative status (policy `trace_rule`):
- Now: primal only.
- Issue #40: `ode_integrate_jvp`, `ode_integrate_vjp`. Sensitivities
  differentiate the frozen accepted-step mesh in `ode_solution_t`.

---

## B-spline basis

| libneo / GSL / FGSL (NEO-2) | fortnum |
|-----------------------------|---------|
| `gsl_bspline_alloc` + `gsl_bspline_knots` | `bspline_init` + `bspline_set_knots` |
| `gsl_bspline_eval` | `bspline_eval_basis(ws, x, values, status)` |
| `gsl_bspline_deriv_eval` | `bspline_eval_deriv(ws, x, nderiv, dvalues, status)` |
| span lookup | `bspline_span_index(ws, x)` |

Clean-room Cox-de Boor recursion (de Boor 2001; Piegl & Tiller); no GSL source
is used. The convention matches NEO-2 `gsl_bspline_routines_mod.f90`: order
`k = degree + 1`, clamped end knots, `ncoef = nbreak + k - 2`. The workspace
`bspline_workspace_t` is caller-owned, replacing the GSL workspace handle and
its free call.

Derivative status (policy `transparent` inside a fixed span): `bspline_eval_jvp`
/ `bspline_eval_vjp` ship now; the span index is `primal_only`. The Taylor
extrapolation above the last breakpoint (`collop_bspline_taylor`) stays in the
NEO-2 consumer.

---

## C ABI for C/C++ consumers

C and C++ callers (KAMEL/KiLCA, MEPHIT) link the static library and include the
hand-written header `include/fortnum.h`, installed to `<prefix>/include`. The
header declares the `bind(c)` entry points in `fortnum_capi` and
`fortnum_capi_bspline`: scalar specials (`fortnum_bessel_in`, `fortnum_dawson`,
`fortnum_gamma_lower`, `fortnum_erf`, `fortnum_erfc`), complex Bessel and 1F1
(`fortnum_bessel_*_complex`, `fortnum_hyperg_1f1`), quadrature and adaptive
integration (`fortnum_gauss_legendre*`, `fortnum_integrate_qag*`),
`fortnum_levin_u_accel`, root finding (`fortnum_root_brent`,
`fortnum_multiroot_hybrid`, `fortnum_deriv_central`, `fortnum_argsort`), the
DOP853 ODE drivers (`fortnum_ode_integrate_dop`, `fortnum_ode_solve_dop`), and a
B-spline handle API (`fortnum_bspline_create` / `_set_knots` / `_eval_basis` /
`_eval_deriv` / `_span_index` / `_destroy`). Status codes cross as plain
integers matching `fortnum_status_t%code`. The header is the single source of
truth for the C signatures.

---

## Polynomial interpolation and grid search

| libneo / math-kit | fortnum |
|-------------------|---------|
| `plag_coeff(nlag, nder, x, xp, coef)` | `lagrange_weights` + `lagrange_deriv_weights` |
| `binsrc(p, nmin, nmax, xi, i)` | `grid_search(p, nmin, nmax, xi, i)` |

Interface differences for `plag_coeff` -> `lagrange_weights` / `lagrange_deriv_weights`:
- `lagrange_weights(n, x, xp, coef)` computes value weights only.
- `lagrange_deriv_weights(n, x, xp, dcoef)` computes derivative weights only.
- Call both when you need both. The split keeps the derivative policy clean:
  value weights are `transparent`; derivative weights use `analytic_rule`.
- The node count `n` is the first argument (not embedded in `nlag`/`nder`).
- Nodes `xp(1:n)` must be distinct; no check is performed.

`binsrc` and `grid_search` are equivalent. The index convention is the same:
`p(nmin:nmax)` strictly increasing; result `i` in `[nmin+1, nmax]` such that
`p(i-1) < xi <= p(i)`.

Derivative status:
- `lagrange_weights` / `lagrange_deriv_weights`: policy `analytic_rule`.
  Now: primal only. Issue #40: `lagrange_weights_jvp`.
- `grid_search`: policy `primal_only` (integer output; no derivative ever).

---

## Random number generation

| libneo / GSL | fortnum |
|--------------|---------|
| `gsl_rng_alloc` + `gsl_rng_set(seed)` | `rng_seed(g, seed, status)` |
| `gsl_rng_uniform(rng)` | `rng_uniform(g, value)` |
| Standard normal via GSL | `rng_normal(g, value)` |
| `gsl_rng_free` | nothing (rng_t lives on the stack or in a derived type) |

Interface differences:
- `rng_t` is a plain derived type; no allocation or free needed.
- `rng_split(parent, stream, child, status)` derives independent per-thread
  substreams from one seeded parent without locking. This replaces a pattern
  of seeding separate generators with different seeds.
- The algorithm is Threefry-2x64-20, not the Mersenne Twister or other GSL
  generators. Output is deterministic: the same seed and the same stream index
  produce the same sequence on any compiler and platform that passes the
  built-in determinism gate.

Derivative status (policy `primal_only`):
- Now: primal only.
- Never: a pseudorandom draw is not a differentiable function of its seed.
  No derivative entry point exists or will exist. Gradients of estimators
  built on draws live in the caller.
