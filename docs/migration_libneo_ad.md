# Migration: libneo derivative products on fortnum

This document extends `docs/migration_libneo.md` with derivative availability.
For each mapping it states what is available NOW (issue #40 landed), what is
planned with a tracking reference, and what is excluded with a reason.

"Preserve OLD primal imports" means a libneo caller that only calls the primal
keeps compiling unchanged; derivative routines are new names added beside the
primal, never added arguments to it (`ad.md §6`).

---

## Special functions

### Modified Bessel functions

| libneo backend | fortnum primal | derivative product | status |
|---|---|---|---|
| modified Bessel `I_n(x)` | `bessel_in(n, x)` | `bessel_in_jvp(n, x, v, jv)` | NOW (`analytic_rule`) |
| modified Bessel `K_n(x)` | `bessel_kn(n, x)` | `bessel_kn_jvp(n, x, v, jv)` | NOW (`analytic_rule`) |
| modified Bessel `I_n(x)` array | `bessel_in_array(nmax, x, values)` | `bessel_in_array_jvp` / `bessel_in_array_vjp` | NOW (`analytic_rule`, recurrence over array) |

Active argument in both JVPs: `x` (the argument). Inactive: `n` (integer
order, selects branch; `ad.md §3`).

Analytic rules:

- `bessel_in_jvp`: `dI_n/dx = (I_{n-1}(x) + I_{n+1}(x)) / 2`  (DLMF 10.29.2)
- `bessel_kn_jvp`: `dK_n/dx = -(K_{n-1}(x) + K_{n+1}(x)) / 2`  (DLMF 10.29.4)

HVP deferred: second-order recurrences involve two further orders; no scalar
loss context in the module.

### Dawson integral

| libneo / math-kit | fortnum primal | derivative product | status |
|---|---|---|---|
| `dawson(x)` | `dawson(x)` | `dawson_jvp(x, v, jv)` | NOW (`analytic_rule`) |
| | | `dawson_grad(x, u, jtu)` | NOW (scalar adjoint) |

Active argument: `x`. Analytic rule: `F'(x) = 1 - 2 x F(x)`.

### Lower incomplete gamma

| libneo backend | fortnum primal | derivative product | status |
|---|---|---|---|
| lower incomplete gamma | `gamma_lower(a, x)` | `gamma_lower_jvp(x, v, jv)` | NOW d/dx (`analytic_rule`) |
| | | `gamma_lower_jvp_da(x, v, jv)` | NOW full [x, a] (`analytic_rule`) |
| regularized lower incomplete gamma P | `gamma_reg_p(a, x)` | `gamma_reg_p_jvp` / `gamma_reg_p_grad` | NOW (`analytic_rule`, digamma series) |

The packed input is `x = [x_val, a_val]`. `gamma_lower_jvp` treats only the
integration limit `x` as active (shape `a` inactive); `gamma_lower_jvp_da` and
the `gamma_reg_p_*` products treat both `x` and `a` as active.

Analytic rules:
- `d/dx gamma_lower(a, x) = x^{a-1} exp(-x)` (DLMF 8.8.1).
- `d/da gamma_lower(a, x) = x^a exp(-x) (ln x . S - T)` (the Gamma(a) factor
  cancels the digamma); `d/da P(a, x) = A[(ln x - psi(a)) S - T]` keeps the
  digamma `psi(a)` (DLMF 5.2.2), with `S`, `T` the regularized-P series sums.

---

## FFT

| libneo / FFTW | fortnum primal | derivative product | status |
|---|---|---|---|
| `fftw_plan_dft_1d` + forward execute | `fft_c2c(z, -1)` | `fft_c2c_jvp(dz, sign)` | NOW (`transparent`) |
| `fftw_plan_dft_1d` + backward execute | `fft_c2c(z, +1)` | `fft_c2c_vjp(u, sign)` | NOW |
| `fftw_plan_dft_r2c_1d` + execute | `fft_r2c(x, c, plan)` | `fft_r2c_jvp(dx, dc, plan)` | NOW |
| | | `fft_r2c_vjp(u, xbar, plan)` | NOW |

The DFT is C-linear; Enzyme differentiates the implementation directly
(`transparent` policy). The Jacobian is the transform itself.

HVP omitted: the DFT Hessian is zero (linear map).

Primal callers that still use `fft_c2c` and `fft_r2c` directly are unaffected;
the JVP/VJP routines are new names.

---

## ODE integration

| libneo | fortnum primal | derivative product | status |
|---|---|---|---|
| `odeint` (external adaptive integrator or custom) | `ode_integrate(problem, ws, sol, st)` | `ode_integrate_jvp` forward sensitivity | NOW (`trace_rule`) |
| | | `ode_integrate_vjp` discrete adjoint | NOW |
| | `ode_solve(rhs, t0, t1, y0, ...)` | primal convenience call | primal only |

Both products re-run the Cash-Karp stepper over the frozen `solution%t` /
`solution%h` schedule that the primal recorded. The adjoint walks the same
trace backward.

Policy `trace_rule`: the derivative is taken at the accepted-step mesh the
primal chose. Changing tolerances changes the trace and the derivative; keep
tolerances fixed when comparing derivative calls to finite difference.

HVP deferred: needs a caller-defined scalar loss on `y(t1)`.

---

## Quadrature / QUADPACK

| libneo backend / QUADPACK | fortnum primal | derivative product | status |
|---|---|---|---|
| adaptive quadrature QAG / `dqag` | `integrate_qag(...)` | `integrate_qag_jvp(dfdp, result, di_dp, st, ctx)` | NOW (`trace_rule`) |
| adaptive quadrature QAGS / `dqags` | `integrate_qags(...)` | `integrate_qags_jvp` | NOW (`trace_rule`) |
| adaptive quadrature QAGP | `integrate_qagp(...)` | `integrate_qagp_jvp` | NOW (`trace_rule`) |
| adaptive quadrature QAGIU | `integrate_qagiu(...)` | `integrate_qagiu_jvp` | NOW (`trace_rule`) |
| Single-rule QK panel | `gk_apply(...)` | primal only | |

`integrate_qag_jvp` differentiates at the frozen accepted subdivision recorded
in `integrate_result_t`. The scalar output makes the reverse product equal the
forward sensitivity; no separate VJP/grad is shipped.

The `integrate` flat convenience call discards the trace; callers that need
a derivative must use `integrate_qag` with explicit `integrate_result_t`.

Fixed-rule `gauss_legendre` is `transparent` in the integrand values (linear
map). Products NOW: `gauss_legendre_jvp`, `gauss_legendre_vjp`,
`gauss_legendre_grad`.

---

## Polynomial interpolation and grid search

| libneo / math-kit | fortnum primal | derivative product | status |
|---|---|---|---|
| `plag_coeff(nlag, nder, x, xp, coef)` | `lagrange_weights(n, x, xp, coef)` | `lagrange_weights_jvp(n, x, xp, f, vx, jv)` | NOW (d/dx active) |
| | | `lagrange_weights_vjp(n, x, xp, f, u, jtu)` | NOW |
| | `lagrange_deriv_weights(n, x, xp, dcoef)` | derivative weights are the analytic rule themselves | NOW (primal) |
| | | `lagrange_fval_jvp(n, x, xp, vf, jv)` | NOW (d/df active, transparent) |
| | | `lagrange_fval_vjp(n, x, xp, u, jtu)` | NOW |
| `binsrc(p, nmin, nmax, xi, i)` | `grid_search(p, nmin, nmax, xi, i)` | none | never (`primal_only`: integer output) |

`lagrange_weights_jvp` and `_vjp` treat `x` as active (d p(x)/dx). The
`lagrange_fval_*` pair treats the nodal values `f` as active (linear, transparent).

Derivatives are valid inside a fixed cell. The caller must hold the cell index
from `grid_search` fixed when differentiating with respect to the evaluation
point; crossing a cell boundary is a non-smooth event (`ad.md §4`).

---

## Random number generation

| libneo backend | fortnum primal | derivative product | status |
|---|---|---|---|
| RNG alloc + seed | `rng_seed(g, seed, status)` | none | never |
| uniform draw | `rng_uniform(g, value)` | none | never |
| standard normal draw | `rng_normal(g, value)` | none | never |
| RNG free | nothing (rng_t lives on stack) | | |

Policy `primal_only`. A pseudorandom draw is not a differentiable function of
its seed (`ad.md §4`). Gradients of estimators built on draws live in the
caller, not in fortnum_rng.
