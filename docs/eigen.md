# Matrix-free eigenproblems

`fortnum_eigen` exports kind-generic `krylov_schur`, `lanczos`, `eigen_bounds`,
and `inertia_count`. FP64 is the reference default through `fortnum_kinds:dp`;
FP32 is a separate overload. Both use the same algorithm include. The
kind-independent options and counters live in `fortnum_eigen_types`.

## Operators and solves

```fortran
use fortnum_kinds, only: dp
use fortnum_eigen
! start(n), values(nev), vectors(n,nev), residuals(nev) are caller-owned.
options = eigen_options_t(max_basis=32, retain=16, max_restarts=100)
call lanczos(apply_k, start, values, vectors, residuals, 1.e-9_dp, &
    options, info, stats, mass=apply_m, solve_mass=solve_m)
```

An operator is `subroutine apply(x,y)` with assumed-shape rank-one arrays,
`intent(in)` and `intent(out)`. It may evaluate a matrix-free JVP at a fixed
state. Solves have an additional `integer, intent(out) :: info`; zero means
success. Callbacks must be linear, dimension preserving, and fixed for a solve.
No finite differences or differentiation through iterations are performed.
The public callback interfaces have names `eigen_real_operator_dp`,
`eigen_complex_operator_dp`, `eigen_real_solve_dp`, and
`eigen_complex_solve_dp`; replace `_dp` by `_sp` for FP32.

Without `mass`, the problem is K x = lambda x. With `mass`, it is
K x = lambda M x with Hermitian positive-definite M, and `solve_mass` applies
M inverse. The basis is orthonormal in the M inner product. M actions are
recomputed during reorthogonalization to avoid cancellation in recursively
updated metric images. Positivity checks detect some invalid metrics; they
cannot certify that an arbitrary callback is globally SPD.

For shift-invert, provide `solve_shift` and `sigma` together. The callback
solves (K - sigma M)y = rhs; FortNum applies M before this solve. The operator
is (K - sigma M) inverse M, and physical eigenvalues are sigma + 1/theta.
`solve_mass` is unnecessary on this path. `lanczos` requires a real shift.
A factorization may be captured by the callback without any sparse-library
dependency in FortNum. Factors, pivoting and inertia belong to
[FortSparse](https://github.com/lazy-fortran/fortsparse).

Complex `krylov_schur` also accepts `transform(x,y,info)` and
`map(theta)` together, for a caller-defined spectral transformation and its
inverse eigenvalue map. These replace the default transformed action; `apply`
and optional `mass` still define the original residual. A noninjective map or
an inaccurate solve can prevent convergence. Nonfinite actions and maps, or
nonzero solve status, return `EIGEN_OPERATOR_FAILED`.

Real Arnoldi callbacks return complex values/vectors, so conjugate pairs are
representable. Real Lanczos returns real values/vectors. Complex Lanczos
returns real values and complex vectors. Real callbacks are extended to complex
vectors by separate real/imaginary actions; this favors correctness over cost.

## Restarting, locking and status

`krylov_schur` uses twice-reorthogonalized Arnoldi and an ordered complex Schur
basis. `lanczos` uses its Hermitian specialization, with full
reorthogonalization and thick Ritz restart. This deliberately stores full
operator images rather than relying on a fragile three-term recurrence.
Only the projected problem uses LAPACK (GEES/TREXC/TREVC or HEEV).
Schur vectors, rather than nonorthogonal eigenvectors, are retained on restart.
Converged directions are soft locked by retaining them and their images;
there is no separate hard-lock tolerance or removal from final residual checks.
At happy breakdown, deterministic coordinate injections expose missing
invariant subspaces and repeated eigenvalues. A partial subspace alone cannot
prove that all requested multiplicities or extremal eigenvalues were found.

`eigen_options_t` supplies `max_basis`, `retain` (zero chooses a default),
`max_restarts`, and `which`: `EIGEN_SMALLEST`/`EIGEN_LARGEST` order physical
real parts, `EIGEN_LARGEST_MAGNITUDE` orders physical magnitudes, and
`EIGEN_NEAREST` orders distance to sigma (zero without shift-invert).
A restarted basis must have `nev <= retain < max_basis <= n`; a full n-vector
basis needs no restart. `stats` records restarts, actions, converged pairs,
and the leading soft-locked subspace size. Callback action counts refer to
complex actions; a real bridge can invoke two underlying real actions.

Acceptance recomputes r = Kx - lambda Mx with ||x||_M = 1, and requires
||r||_2 <= tolerance (k_scale + |lambda| m_scale) ||x||_2.
This is the normwise backward-error test for perturbations of K and M when
`k_scale = ||K||_2` and `m_scale = ||M||_2`. Known norm bounds can be supplied
as the nonnegative FP64 options `operator_scale` and `mass_scale`, including
for FP32 solves (they must be representable in the working kind). Negative
or nonfinite scales are invalid. Each zero option independently selects an
automatic estimate: the largest observed ||Kq||_2/||q||_2 or
||Mq||_2/||q||_2 over basis vectors and checked Ritz vectors; M=I uses one.
These are lower estimates, not certified matrix norms, so automatic stopping
can be stricter than stopping with exact norms. Unseen directions need not
be represented. No absolute unit floor is imposed, preserving scaling for
small operators as well as highly cancelling Kx and lambda Mx.

Original K images are reused on the ordinary path. Shift-invert and custom
transforms evaluate K once more on each new basis vector when
`operator_scale=0`; these actions are included in `stats%applications`.
Supplying `operator_scale` avoids these extra actions. Workspace is reused;
no matrix entries or assembled operators are needed.
`residuals` always reports the original Euclidean residual, including after
shift-invert. It is a floating-point estimate, without an outward rounding
certificate. The returned vector and residual kind follows the input kind;
FP32 output must be independently rechecked in FP64 for a precision-policy
acceptance decision. No reduced-precision policy is silently selected here.

`EIGEN_OK` means every requested pair passes this residual test.
`EIGEN_NOT_CONVERGED` returns the last candidates at the iteration limit;
`EIGEN_INVALID` rejects invalid dimensions, options or callback combinations;
`EIGEN_METRIC_FAILED` rejects observed metric/self-adjointness failures;
`EIGEN_DENSE_FAILED` reports a projected LAPACK failure. Outputs after a
negative status must not be accepted as eigenpairs.

Workspace is allocated once per call: eight n-by-max_basis complex arrays,
O(n) vector scratch and O(max_basis squared) projected storage. There is no
allocation in expansion/restart loops. BLAS performs basis rotations; OpenMP
parallelizes large reorthogonalization updates. Small problems remain serial.
Mixed storage, CUDA, block starts and distributed operators are future work.

## Conditional error bounds and counts

For r = Kx - theta Mx, a supplied positive lower bound m0 for lambda_min(M)
and `metric_norm = ||x||_M` give rho = ||r||_2 / (sqrt(m0) metric_norm).
For a Hermitian pencil, an eigenvalue lies within rho of theta. This does not
identify the eigenvalue or establish extremal ordering. For a diagonalizable
non-Hermitian pencil, supply `eigenvector_condition >= 1`, a trusted upper
bound for the eigenbasis condition number of M inverse-square-root K M
inverse-square-root; the Bauer-Fike radius is that number times rho. A
defective pencil has no such finite eigenbasis bound.

For a normalized Hermitian Rayleigh quotient theta and a trusted interval
(a,b) containing exactly one eigenvalue and no other spectrum, Kato-Temple
provides theta - rho squared/(b-theta) and theta + rho squared/(theta-a).
Supply `rayleigh`, `gap_lower`, `gap_upper`, and optionally `lower`, `upper`.
The routine checks a < theta < b and rho squared < (theta-a)(b-theta).
It does not establish the spectral-gap premise. Bauer-Fike and Temple modes
cannot be combined. Nonfinite/unsupported inputs return `EIGEN_INVALID`.
These are conditional floating-point estimates, not validated enclosures.

`inertia_count(factorization,n,left,right,count,info)` requests negative,
zero and positive inertia for K - s M at both endpoints. It checks statuses,
nonnegative counts, dimension sums and monotonicity. The difference of negative
counts gives the number of eigenvalues in [left,right), including
multiplicity. This applies only to a Hermitian pencil with SPD M. A supplied
floating-point factorization needs its own pivot/sign certification before the
count constitutes a rigorous proof of spectral ordering.

## Evidence and regeneration

The committed FortSym generator proves the Toeplitz exponential-mode
recurrence, the interval Laplacian identity, the axisymmetric disc Bessel
series recurrence, generalized shift-invert cancellation, and the bound
identities with symbolic zero tests. It emits radius and Temple leaves in both
kinds, plus contracted Temple JVP/VJP leaves. Generation is not a proof of the
spectral theorems or a floating-point rounding certificate.

Native tests cover restarted Toeplitz sine spectra, congruent SPD pencils,
real/complex nonnormal triangular spectra, real/complex Hermitian cases,
shift-invert ordering, custom maps, multiplicities, callback failures,
nonconvergence, inertia endpoints, both kinds, and a radial disc finite-volume
limit against roots of the exact J0 power series. The n=32768 disc regression
uses backward-error tolerance 1e-10; the separate eigenvalue comparison allows
1% spatial discretization error. A scaled Toeplitz regression checks that a
small unresolved operator still fails, and verifies supplied-scale acceptance
against the analytic norm bound and sine spectrum. Run `ctest -R fortnum_eigen --output-on-failure` or
`fo test 'test_fortnum_eigen*'`. Regenerate with
`cd tools/codegen && fo exec gen_eigen_bounds`; the registered generated-source
check byte-compares its output.

Algorithm references: [Netlib implicit restart templates](https://www.netlib.org/utk/people/JackDongarra/etemplates/node222.html),
[LAPACK GEES](https://www.netlib.org/lapack/explore-html/d5/d38/group__gees_gac426d99ddc8ac986de2a17e7ad3f2512.html),
and [LAPACK TREVC](https://www.netlib.org/lapack/explore-html/da/d6a/group__trevc_ga2d7a7abf901cb7ff1a8ade8d6f11ac80.html).
