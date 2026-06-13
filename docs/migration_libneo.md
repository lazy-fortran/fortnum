# Migration: libneo external backend and math-kit to fortnum

This document maps the math routines used in libneo's external numerical backend and
math-kit to their fortnum equivalents. For each mapping it states whether
derivative products are available now, planned, or excluded.

The fortnum derivative policy vocabulary is defined in `docs/design/ad.md`.
No derivative products ship yet; all are reserved for issue #40.

---

## Special functions

### Modified Bessel functions

| libneo backend | fortnum |
|--------------|---------|
| `modified Bessel `I_n(x)`` | `bessel_in(n, x)` |
| `modified Bessel `I_n(x)` array` | `bessel_in_array(nmax, x, values(0:nmax))` |
| `modified Bessel `K_n(x)`` | `bessel_kn(n, x)` |

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

| libneo backend | fortnum |
|--------------|---------|
| `lower incomplete gamma` (lower incomplete, unnormalized) | `gamma_lower(a, x)` |
| `regularized lower incomplete gamma P` (regularized P) | `gamma_reg_p(a, x)` |

Both routines are pure functions. Domain violations (`a <= 0` or `x < 0`) stop
with `error stop`; the primal signature carries no status argument by design.
A status path may be added later without changing the signature.

Derivative status (policy `analytic_rule`):
- Now: primal only.
- Issue #40: `gamma_lower_grad`, `gamma_reg_p_grad` using
  d/dx gamma(a,x) = x^{a-1} exp(-x).

---

## Numerical integration

| libneo backend / QUADPACK | fortnum |
|--------------------------|---------|
| `adaptive quadrature QAG` | `integrate_qag` |
| `adaptive quadrature QAGs` | `integrate_qags` |
| `adaptive quadrature QAGp` | `integrate_qagp` |
| `adaptive quadrature QAGIU (semi-infinite)` | `integrate_qagiu` |
| Direct QUADPACK `dqag` / `dqags` | `integrate_qag` / `integrate_qags` |
| Single-rule QK15/21/31/61 panel | `gk_apply` |

Interface differences:
- The integrand interface requires an optional `ctx` argument:
  `function f(x, ctx) result(fx)`. A parameter-free integrand ignores `ctx`.
  This replaces the reference's `C integrand struct and function pointer.
- Workspace and epsilon-table types are caller-owned and explicit.
  `integrate_workspace_t` replaces the reference workspace handle;
  `integrate_epstab_t` is the Wynn epsilon table. Pass default-initialized
  values on the first call; reuse them across calls to avoid reallocation.
- Status is a `fortnum_status_t` argument, not an integer error code.
  `FORTNUM_OK` maps to the reference's `success code 0`; `FORTNUM_CONVERGENCE_ERROR`
  maps to `max-iteration or round-off-limit failures`; `FORTNUM_DOMAIN_ERROR` covers
  non-smooth derivative cases (break-point singularity at an endpoint, or
  bad integrand behaviour).
- `integrate` is a flat single-call wrapper that owns its workspace; it is
  the closest equivalent to a the reference integration call without explicit workspace
  management.

Derivative status (policy `trace_rule`):
- Now: primal only.
- Issue #40: `integrate_qag_jvp`, `integrate_qags_jvp`, `integrate_qag_vjp`,
  `integrate_qags_vjp`, `integrate_qag_grad`. Derivatives are taken at the
  frozen accepted subdivision recorded in `integrate_result_t`.

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
| `odeint` (the reference `the adaptive evolve routines` or custom) | `ode_integrate` / `ode_solve` |

Interface differences:
- The RHS is an abstract interface: `subroutine rhs(t, y, dydt, ctx)` with
  an optional unlimited polymorphic `ctx`. Thread an arbitrary parameter
  struct through `ctx` instead of a `void *` or a module variable.
- `ode_solve` is the flat call for callers that do not need the workspace or
  the recorded trace.
- `ode_at` evaluates the solution at a specified set of output times, matching
  the common pattern of requesting output at prescribed points.
- The method is Cash-Karp RK5(4) with PI step control. If libneo used a
  different order or tableau, check that the tolerances give the needed
  accuracy; the step controller uses the same `rtol` and `atol` parameters.

Derivative status (policy `trace_rule`):
- Now: primal only.
- Issue #40: `ode_integrate_jvp`, `ode_integrate_vjp`. Sensitivities
  differentiate the frozen accepted-step mesh in `ode_solution_t`.

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

| libneo backend | fortnum |
|--------------|---------|
| `RNG alloc + seed` | `rng_seed(g, seed, status)` |
| `uniform draw` | `rng_uniform(g, value)` |
| standard normal draw | `rng_normal(g, value)` |
| `RNG free` | nothing (rng_t lives on the stack or in a derived type) |

Interface differences:
- `rng_t` is a plain derived type; no allocation or free needed.
- `rng_split(parent, stream, child, status)` derives independent per-thread
  substreams from one seeded parent without locking. This replaces a pattern
  of seeding separate generators with different seeds.
- The algorithm is Threefry-2x64-20, not the Mersenne Twister or another
  generators. Output is deterministic: the same seed and the same stream index
  produce the same sequence on any compiler and platform that passes the
  built-in determinism gate.

Derivative status (policy `primal_only`):
- Now: primal only.
- Never: a pseudorandom draw is not a differentiable function of its seed.
  No derivative entry point exists or will exist. Gradients of estimators
  built on draws live in the caller.
