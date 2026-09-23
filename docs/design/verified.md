# Verified computing layer

`fortnum` provides a generic certified-computing runtime: arithmetic whose
every result is a mathematical enclosure of the exact real or complex result,
plus series, transforms, and linear algebra built on it. Two applications
motivate the layer. Both carried private copies of the same machinery:

- `kinetic-compression`: certified neoclassical transport (complex balls,
  ball Fourier series with l1 tails, FFT convolution with a rigorous error
  bound, Wiener-algebra reciprocals, verified linear algebra, quad intervals,
  complementary variational bounds, Lean package `KineticCert`).
- `gc-loss-certificate`: certified alpha-particle loss (real intervals,
  rectangular complex intervals, interval dual numbers, Lohner/Taylor
  validated ODE enclosures, analytic Bernstein field models, strip decay,
  Lean package `LossCert`).

This document records the inventory of both projects, the semantics adopted
by `fortnum`, the module layout, the split between `fortnum`, `fortsym`,
`libneo`, and a shared Lean package, and the migration order.

## Inventory

Line counts are from the project sources at the time of the inventory.
"Generic" means the code has no physics or file-format dependence.

| Area | kinetic-compression | gc-loss-certificate | Overlap |
| --- | --- | --- | --- |
| Directed rounding primitives | 15 private `up`/`dn` copies (about 90 lines), `nearest` except `kc_ball` and `kc_fseries` (successor formula) | 2 in `interval.f90`, `nearest` | same function, two methods |
| Real interval | `ivl_t` in `kc_collop` (about 90 lines) | `interval.f90` (322) | full |
| Quad-precision interval | `kc_qivl` (138) plus a private copy in `kc_linlandau` | none | internal duplicate |
| Complex enclosure | `kc_ball` (212, midpoint-radius disc) | `cint_t` in `cinterval.f90` (about 200, rectangle) | same role, two representations |
| First-order interval AD | none | `idual.f90` (237), `cdual_t` in `cinterval.f90` (about 165) | generic, one user |
| Taylor-series arithmetic | none | `taylor_c.f90` (80) | generic, one user |
| Fourier series with tails | `kc_fseries` (288), `kc_fseries2d` (448), `kc_fourier` (65) | aliasing bounds in Lean only | generic, one user |
| Rigorous FFT convolution | `kc_cfft` (378) | none | generic, one user |
| Wiener-algebra reciprocal | inside `kc_neo_geom` and `kc_neo_geom3d` (about 60 lines each) | none | internal duplicate |
| Verified linear algebra | `kc_verified_la` (229), Cholesky checks in `kc_fullop_mat` | `mm`, `mv`, `qr_frame`, `inverse_enclosure`, `step_jacobian` in `lohner7.f90`, repeated for n = 8 in `lohner8.f90` (about 110 each) | partial, plus internal duplicate |
| Validated ODE | none | `lohner7`, `lohner8`, `lohner_taylor`, `lohner_time`, `flow_enclosure` (1252 total, about half generic) | generic core, one user |
| Polynomial field models | none | Bernstein evaluation, `bezier_hull`, local re-expansion in `cheb_field.f90` (342, about 150 generic) | generic core |
| Boozer input | `kc_boozer_bc` (425, `.bc` reader) | `boozer_field.f90` (379), `booz_*` | both belong in `libneo` |
| Lean | `DiscLemma`, `StieltjesConvex`, `VariationalBounds`, `WienerInverse` (487) | `StripDecay`, `Aliasing`, `Lohner` (1332) | same toolchain and Mathlib revision |

About 2900 lines of generic Fortran and 1800 lines of generic Lean live in the
two projects. Roughly 800 Fortran lines are direct duplicates (rounding
primitives, real intervals, complex enclosure arithmetic, the n = 7/n = 8
Lohner helpers, the two Wiener reciprocals).

## Semantics

### Directed rounding without libm

All enclosures run in IEEE binary64 round-to-nearest, without changing the
rounding mode. After each correctly rounded operation the result is moved
outward by the successor formula of Rump, Zimmermann, Boldo, and Melquiond
("Computing predecessor and successor in rounding to nearest", BIT 49, 2009,
Algorithm 2):

```text
up(x) = fl(x + fl(phi |x| + eta)),  phi = u (1 + 2u),  u = 2^-53,
eta = 2^-1074
```

For every finite `x` without flush-to-zero, `up(x)` is `succ(x)` or
`succ(succ(x))`, hence `up(x) >= nearest(x, +1)`. `dn(x) = -up(-x)`.
`kinetic-compression` uses this formula. `gc-loss-certificate` uses the
`nearest` intrinsic, which is exact (one ulp) but branches on exceptional
values and does not vectorize. `fortnum` adopts the successor formula as the
default because it is branch-free and needs no library call; both are
rigorous. Infinite inputs stay infinite and NaN propagates, so an overflowed
enclosure degrades to the whole line instead of becoming wrong.

The build must keep IEEE semantics: no `-ffast-math`, no
`-funsafe-math-optimizations`, no flush-to-zero. Fused multiply-add
contraction does not break the enclosures: the proof of `up(x) >= x` only
needs `fl(phi |x| + eta) > ulp(x)/2`, which holds whether the inner
expression is rounded once (fused) or twice, and a fused final addition
`fl(a b + c)` with `c > ulp(a b)/2` is still at least `a b`.

### Trust assumptions

`fortnum` verified modules use only IEEE `+`, `-`, `*`, `/`, and `sqrt`,
which are correctly rounded. Transcendental functions (`exp`, `log`, `sin`,
`cos`, `sinh`, `cosh`) and the constant pi are enclosed by Taylor series with
explicit Lagrange remainders evaluated in interval arithmetic, after exact
argument reduction against two-word enclosures of `ln 2` and `pi`. No libm
accuracy claim enters.

The projects currently trust libm in these places, which the migration
removes:

- `gc-loss-certificate`: `exp`, `cos`, `sin`, `sinh`, `cosh`, `acos(-1)` to
  within one ulp (widened by two `nearest` steps).
- `kinetic-compression`: `cos`/`sin` twiddle factors to within four units of
  roundoff in `kc_cfft`, `exp` for the sigma weights in `kc_fseries2d`, and
  `hypot` (`abs` of a complex number) to within one ulp outside the safe
  exponent range in `kc_ball`.

### Findings from the review

These gaps do not invalidate current project results in practice (margins
dominate), but each one is an unproved step. The `fortnum` versions close
them.

1. `kc_cfft::rconv_pow2` bounds an l1 norm by `up(sum(aa))`. A floating sum of
   `n` terms can undershoot the exact sum by `(n - 1) u` relative, which one
   upward step does not cover. `fortnum` accumulates with an upward step per
   addition.
2. `kc_cfft::compose_conv_err` forms `rho = up(a)/dn(b)` without rounding the
   quotient upward.
3. `kc_fseries::fs_mul_fast` feeds `abs(c)` (round-to-nearest `hypot`) into
   the nonnegative radius convolution, which needs an upper bound of `|c|`.
4. `kc_verified_la::verified_inverse` forms `R Z - I` with an unanalysed
   floating subtraction on the diagonal; `norm2_tight` factors
   `fl(tau + marg - h_ii)` and certifies a matrix that differs from the
   intended one by that rounding.
5. `gc-loss-certificate` `lohner7::inverse_enclosure` and `step_jacobian`
   round `eps` and `x` with fixed factors `1 + 1e-10` and `1 + 1e-12` instead
   of directed rounding.
6. `kc_fseries2d::wgt` rounds a libm `exp` upward once; tail estimates that
   divide by the weight need a lower bound instead.

## Representations

| Type | Module | Set | Use |
| --- | --- | --- | --- |
| `interval_t` | `fortnum_interval` | `[lo, hi]` real | scalar enclosures, interval Newton, a priori boxes |
| `cinterval_t` | `fortnum_interval` | `re + i im`, both intervals | complex boxes for strip bounds (`gc-loss-certificate` `cint_t`) |
| `ball_t` | `fortnum_ball` | disc `abs(z - c) <= r` | Fourier coefficients, complex linear algebra (`kc_ball`) |
| `idual_t` | `fortnum_idual` | interval value plus interval gradient | box Jacobians for Lohner steps (`gc-loss-certificate` `idual_t`) |
| `fseries_t` | `fortnum_fseries` | 1D ball Fourier series plus l1 tail | `kc_fseries` |
| `fseries2d_t` | `fortnum_fseries` | 2D series with sigma-weighted l1 tail | `kc_fseries2d` |

Rectangles and discs are both kept. A disc is the natural object for
products and reciprocals (exact disc image under inversion, `DiscLemma`);
a rectangle is tighter for separable functions such as `sin(x + iy)`.
Conversion in both directions is provided.

`idual_t` carries a fixed-capacity gradient (`idual_max_vars = 8`) and an
active count `n`. The capacity avoids allocation in hot loops;
`gc-loss-certificate` uses four and eight variables. A complex dual type
(`cdual_t`) follows the same pattern and is scheduled with the validated ODE
module, which is its only consumer.

## Module layout

All verified sources live under `src/verified/`.

| Module | Content | Stage |
| --- | --- | --- |
| `fortnum_rounding` | `round_up`, `round_down` (successor formula), `sum_up`, `gamma_up` (Higham `gamma_k`) | step 2 |
| `fortnum_interval` | `interval_t`, `cinterval_t`; `+ - * /`, `sqr`, `sqrt`, `pow`, `exp`, `log`, `sin`, `cos`, `sinh`, `cosh`, `pi`; hull, midpoint, width, magnitude, containment | step 2 |
| `fortnum_ball` | `ball_t`; `+ - * /`, real scaling, exact-disc reciprocal, modulus bounds, conversions | step 2 |
| `fortnum_idual` | `idual_t` interval forward AD | step 2 |
| `fortnum_cfft_rigorous` | radix-2 FFT with certified twiddles, Higham error factor, 1D/2D convolution error bounds, upper-bound convolution of nonnegative sequences | step 2 |
| `fortnum_fseries` | 1D/2D ball Fourier series with tails: add, scale, product (direct and FFT), derivative, weighted inner products, Wiener reciprocal | step 2 |
| `fortnum_verified_linalg` | matrix balls, rigorous norm bounds, product error bounds, verified inverse, Cholesky-certified eigenvalue lower bounds | step 3 |
| `fortnum_taylor_series` | order-by-order Taylor arithmetic on intervals and complex duals | planned with `validated_ode` |
| `fortnum_validated_ode` | Lohner QR enclosure and high-order Taylor-Lohner step with an abstract right-hand side (reverse communication or deferred binding), Picard a priori boxes, Gronwall variational bounds, interval Newton event crossings | planned |
| `fortnum_bernstein` | interval Bernstein evaluation, convex-hull bounds, local re-expansion | planned |
| `fortnum_stieltjes` | Stieltjes/Pick bounds: Pade-type two-sided bounds for `c^T (A + z B)^-1 c` from moments, convexity and monotonicity certificates | planned after `kinetic-compression` settles the algorithm |
| `fortnum_interval_qp` | real128 intervals for closed-form Gaussian integrals (`kc_qivl`) | planned |

### Rigorous twiddle factors

`kc_cfft` trusts libm `cos`/`sin` to four units of roundoff. `fortnum`
removes the assumption. For `n = 2^p`, each `k` is reduced exactly by the
eighth-turn symmetries to `0 <= k' <= n/8`, so `theta = 2 pi k'/n` lies in
`[0, pi/4]`. `cos theta` and `sin theta` are enclosed by their Taylor series
to degree 26 with a Lagrange remainder bound, in interval arithmetic, with
`pi` enclosed by `[fl(pi), succ(fl(pi))]`. The stored twiddle is the midpoint
of the enclosure and the per-twiddle error `mu` is computed as the maximum of
`abs(w_stored - w_exact)` over all twiddles, bounded from the enclosure
radii. The FFT error factor then uses this certified `mu` (about `u`)
instead of the assumed `4u`.

### FFT convolution bound

The non-asymptotic bound of Higham, *Accuracy and Stability of Numerical
Algorithms*, 2nd ed., Theorem 24.2, for the radix-2 transform with
per-twiddle error `mu`:

```text
||fl(F x) - F x||_2 <= kappa ||F x||_2,
kappa = p eta / (1 - p eta),  eta = mu + gamma_4 (sqrt 2 + mu),  p = log2 n
```

composed through two forward transforms, the pointwise product (complex
multiplication error `sqrt 5 u`, Brent, Percival, Zimmermann 2007), and the
inverse transform, bounds every entry of a computed linear convolution by
`Gamma(n) ||a||_1 ||b||_1`. The 2D bound composes row and column stages:
`kappa2 = kappa_m + kappa_n + kappa_m kappa_n`. The derivation is in the
source header of `fortnum_cfft_rigorous` and is the one from `kc_cfft`, with
the rounding gaps listed above closed.

## Library split

| Capability | Owner | Notes |
| --- | --- | --- |
| Interval, ball, dual arithmetic, rounding | `fortnum` | runtime interface targeted by `fortsym` emission |
| Series with tails, rigorous FFT, Wiener reciprocal | `fortnum` | generic Banach-algebra operations |
| Verified linear algebra, eigenvalue certificates | `fortnum` | |
| Validated ODE, Taylor arithmetic, Bernstein models | `fortnum` | abstract right-hand side; the physics system stays in the project |
| Emission of interval and ball kernels | `fortsym` | the `interval-emit` target emits code against `fortnum_interval` and `fortnum_ball` |
| Symbolic bracket mechanism | `fortsym` | from an operator specification derive the primal and complementary functionals, their constraints, and the proof obligations; emit the evaluation kernels |
| Rigorous Boozer field enclosures | `libneo` | interval/ball evaluation of Boozer spectra, spline profiles with enclosure |
| `.bc` and `booz_xform` readers | `libneo` | with round-trip tests (read, write, read, compare bitwise) |
| Guiding-centre and drift-kinetic systems, collision operators, certificates | projects | |

`fortsym` owns symbolic derivation and code emission only. A particular
kernel (for example the guiding-centre right-hand side on a Bernstein field)
is generated in the project against a `fortsym` build and linked against
`fortnum`.

### Runtime interface for emitted code

Emitted rigorous kernels need, for each of `interval_t` and `ball_t`: a
constructor from a point value, elemental `+ - * /` with mixed real and
integer operands, integer powers, `sqrt`, `exp`, `log`, `sin`, `cos`, and
the constant pi. `fortnum_interval` extends the intrinsic generic names
`sqrt`, `exp`, `log`, `sin`, `cos`, `sinh`, `cosh`, `abs` so emitted code that
calls them on interval arguments resolves without renaming. `fortsym` is
defining its runtime-interface document on its `interval-emit` branch; the
interface names above are the `fortnum` side and any mismatch is resolved by
a thin adapter module in `fortnum`, never by weakening the rounding
semantics.

## Shared Lean package

A `CertCore` Lake library, required by both `KineticCert` and `LossCert`
(same toolchain and Mathlib revision today), holds the generic theorems:

| File | From | Statement |
| --- | --- | --- |
| `StripDecay` | `LossCert` | exponential Fourier decay of functions analytic in a strip, sup bounds |
| `Aliasing` | `LossCert` | DFT aliasing error of sampled series with exponential decay |
| `WienerInverse` | `KineticCert` | Neumann-series inverse in a Banach algebra |
| `DiscLemma` | `KineticCert` | disc images under affine maps and inversion |
| `VariationalBounds` | `KineticCert` | Ritz lower and complementary upper bounds, polarisation |
| `StieltjesConvex` | `KineticCert` | monotonicity and convexity of Stieltjes resolvent forms |
| `Lohner` | `LossCert` | step confinement, Taylor-Lagrange step, mean-value update, Gronwall |
| `FFTBound` | new | Higham Theorem 24.2 composition used by `fortnum_cfft_rigorous` |
| `Rounding` | new | successor-formula lemma as an axiom over a stated IEEE model |

Project-specific files (`Barrier`, `Tube`, `Grazing`, `Collisions`,
`Hermite`, and the certificates) stay in the projects.

## Migration

Each project migrates module by module. Every step keeps the project's own
certificate tests green and compares the certificate numbers before and after
(new bounds must contain the old enclosures' exact targets and be no wider
than the old bounds plus the documented rounding fixes).

`kinetic-compression`:

1. `kc_ball` -> `fortnum_ball`. Guard: `test_fseries_fast`,
   `test_fseries2d`, `test_cfft_bound`.
2. Private `up`/`dn` copies and `ivl_t` -> `fortnum_rounding`,
   `fortnum_interval`. Guard: `test_collop`, `test_weights`, full `fo test`.
3. `kc_cfft` -> `fortnum_cfft_rigorous` (drops the libm twiddle assumption).
   Guard: `test_cfft_bound`, `test_fseries_fast`, `test_fseries2d_fast`.
4. `kc_fseries`, `kc_fseries2d`, Wiener parts of `kc_neo_geom*` ->
   `fortnum_fseries`. Guard: `test_neo_geom3d`, `test_neo_bz3d`,
   `test_neo_full`.
5. `kc_verified_la` -> `fortnum_verified_linalg`. Guard: `test_realeq`,
   `test_ntv_realeq`, `test_linlandau`.
6. `kc_qivl` and the `kc_linlandau` copy -> `fortnum_interval_qp` once
   implemented. Guard: `test_linlandau`, `test_landau_capped`.
7. `kc_boozer_bc` -> `libneo` reader. Guard: `test_boozer_bc`,
   `test_neo_bz3d_bc`.

`gc-loss-certificate` (after its current work lands):

1. `interval` -> `fortnum_interval`. Guard: `test_model`, `test_boozer`.
2. `cinterval` (`cint_t`) -> `fortnum_interval` (`cinterval_t`); `idual` ->
   `fortnum_idual`. Guard: `test_taylor`, `test_corrected`.
3. `lohner7`/`lohner8` helpers -> `fortnum_verified_linalg` and, once
   implemented, `fortnum_validated_ode`. Guard: `test_enclosure_slow`,
   `test_taylor_time_slow`.
4. Bernstein core of `cheb_field` -> `fortnum_bernstein`; Boozer parts ->
   `libneo`. Guard: `test_bern`, `test_boozer`.
