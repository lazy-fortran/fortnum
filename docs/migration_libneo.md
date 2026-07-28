# libneo primal migration

Status: current mapping guide.

Replace external-backend and math-kit calls with the public `fortnum` modules
below. Compare numerical conventions before changing a call site.

## Special functions

| Operation | `fortnum` interface | Notes |
| --- | --- | --- |
| modified Bessel \(I_n(x)\) | `bessel_in(n, x)` | real order index |
| sequence \(I_0,\ldots,I_n\) | `bessel_in_array(nmax, x, values)` | caller-sized output |
| modified Bessel \(K_n(x)\) | `bessel_kn(n, x)` | positive real argument |
| Dawson integral | `dawson(x)` | real argument |
| lower incomplete gamma | `gamma_lower(a, x)` | unnormalized |
| regularized lower gamma | `gamma_reg_p(a, x)` | \(P(a,x)\) |
| complex Bessel \(J_\nu\) | `bessel_j_complex` | status return |
| complex modified Bessel \(I_\nu\) | `bessel_i_complex` | scaled option |
| complex modified Bessel \(K_\nu\) | `bessel_k_complex` | scaled option |
| confluent hypergeometric | `hyperg_1f1` | general real parameters |
| specialized \(a=1\) | `hyperg_1f1_a1` | optimized path |
| modified specialized form | `hyperg_1f1m_a1` | compatibility form |
| error function C ABI | `fortnum_erf`, `fortnum_erfc` | `bind(c)` |

Check scaling conventions for complex Bessel calls. A scaled result is not
interchangeable with the unscaled function.

## Quadrature and integration

| Existing need | `fortnum` interface |
| --- | --- |
| Gauss-Legendre nodes and weights | `gauss_legendre` |
| mapped Gauss-Legendre rule | `gauss_legendre_ab` |
| generalized Gauss-Laguerre | `gauss_gen_laguerre` |
| one Gauss-Kronrod panel | `gk_apply` |
| simple finite integral | `integrate` |
| adaptive finite rule selection | `integrate_qag` |
| singular or extrapolated finite integral | `integrate_qags` |
| known breakpoints | `integrate_qagp` |
| infinite interval | `integrate_qagiu` |
| CQUAD-style integration | `integrate_cquad` |
| sequence acceleration | `levin_u_accel` |

Stateful drivers use caller-owned workspace and result types. Preserve status
handling and tolerance meaning during migration. For QAGP, breakpoints must be
strictly interior and distinct.

## FFT

| Operation | `fortnum` interface |
| --- | --- |
| complex transform | `fft_c2c(z, sign)` |
| real-to-complex transform | `fft_r2c(x, c, plan)` |
| reusable plan | `fft_plan_init(plan, n)` |

Verify sign and normalization conventions with a direct-transform oracle at
the call site.

## ODEs

Choose the interface that matches the existing solver contract:

| Need | Interface |
| --- | --- |
| adaptive Cash-Karp trace | `ode_integrate` |
| allocating Cash-Karp call | `ode_solve` |
| requested output times | `ode_at` |
| DOP853 | `ode_integrate_dop`, `ode_solve_dop` |
| stateful RK8PD | `rk8pd_evolve_init`, `rk8pd_evolve_apply` |
| Adams DDEABM | `ddeabm_init`, `ddeabm_integrate_to` |
| VODE-style integration | `vode_init`, `vode_integrate_to` |

Map tolerances, maximum steps, event direction, and terminal-event behavior
explicitly. Do not assume dense-output or event semantics from another solver.

## Roots

| Operation | `fortnum` interface |
| --- | --- |
| bracketed scalar root | `root_bisect` or `root_brent` |
| scalar Newton with derivative callback | `root_newton` |
| dense vector root | `multiroot_hybrid`, `multiroot_hybrids` |
| reverse-communication vector root | `multiroot_rc_init`, `multiroot_step` |
| complex zeros in a rectangle | `complex_region_roots` |
| central derivative helper | `deriv_central` |
| index sorting helper | `argsort` |

Preserve bracket, tolerance, iteration, and singular-Jacobian behavior.

## Interpolation and splines

| Operation | `fortnum` interface |
| --- | --- |
| ordered-grid cell lookup | `grid_search` |
| Lagrange weights | `lagrange_weights` |
| derivative weights | `lagrange_deriv_weights` |
| B-spline workspace | `bspline_workspace_t` |
| knot setup | `bspline_set_knots` |
| basis values | `bspline_eval_basis` |
| basis derivatives | `bspline_eval_deriv` |

Grid-cell and knot-span selection are discrete. Preserve endpoint conventions
and caller array bounds.

## Linear algebra

Small fixed-size call sites can use `det2`, `det3`, `inv2`, `inv3`, and
`lu_factorization_t`. `LINALG_MAX_N` bounds the internal small dense solver.
Large or specialized systems should continue using their established BLAS and
LAPACK path.

## Random numbers

Replace implicit global RNG state with caller-owned `rng_t`:

```fortran
call rng_seed(generator, seed, status)
call rng_split(generator, stream_id, child, status)
call rng_uniform(child, value)
```

Assign stream IDs from stable logical work identifiers. A migration changes the
random stream unless the old code already used the same Threefry contract.

## C and C++ callers

Use `include/fortnum.h`. B-spline and RK8PD wrappers return opaque handles.
Destroy each handle with its matching release function. Check every integer
status code before consuming outputs.

## Acceptance

For each migrated call site:

1. match mathematical conventions and status behavior
2. compare values with an independent oracle
3. cover boundary and failure cases
4. benchmark the complete downstream workload
5. migrate derivative products separately with
   [migration_libneo_ad.md](migration_libneo_ad.md)
