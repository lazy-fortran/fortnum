# Error-control primitives

The immutable-matrix workhorse and coordinate sampler support a KIN6D
consumer without changing its Hamiltonian, residual theorem or probability
semantics. Both allocate preparation state once; matrix actions allocate
no scratch or second matrix.

## Point matrix acting on interval data

`prepare_matvec_bound` encloses each sum `s_i = sum_j abs(A_ij)` in
`matvec_bound_t`. The plan applies only to the immutable binary64 matrix
used at preparation. The caller owns that identity contract; changing matrix
bytes requires a new plan. KIN6D stores its matrix privately.

For interval inputs choose their exactly represented lower endpoints `l_j`,
and enclose `w = max_j (x_j.hi-x_j.lo)` and `m = max_j abs(l_j)`.
`enclose_matvec` evaluates a sequential point dot product and uses

```
error_i = s_i * (gamma_(2n) * m + w) + underflow
underflow = 2n * tiny * (1 + gamma_(2n))
```

Every budget operation and final endpoint is outward. There are at most
`2n` rounded multiply/add operations on a summand path; the usual product of
relative errors gives `gamma_(2n)`. The input-width contribution follows
from the triangle inequality. The absolute error of each underflowed
operation is bounded by `tiny` (smallest normal), deliberately larger than
a subnormal rounding unit, and amplified by the same product bound.
Signs ±1 are exact; optional accumulation adds intervals outward.

Premises: IEEE binary64, round to nearest, gradual underflow, no unsafe
reassociation or contraction. Compile with `-fno-fast-math -ffp-contract=off`.
Unsupported dimensions, invalid intervals and nonfinite preparation/results
fail explicitly. A plan prepared under these premises does not authorize
changing the process rounding/flush mode during use.

Independent tests use binary128 box vertices, signed cancellation,
gradual underflow and overflow rejection. CMake registers `verified_matvec`.

## Randomized norm of a fixed deterministic vector

`sampled_norm_upper` requests `m` interval coordinate evaluations through
a callback. The caller supplies a proved global envelope `M >= |v_i|`
for every coordinate. A sampled maximum cannot establish this premise.
Sampling `I_j` independently and uniformly from `1..N` gives
`X_j = (v_(I_j)/M)^2` in `[0,1]`. Hoeffding's inequality gives

```
P(E X > mean(X) + sqrt(log(1/delta)/(2m))) <= delta
||v||_2 <= M * sqrt(N * min(1, mean(X) + tail))
```

Interval callback values bound each sampled square above. Log, square root,
reduction and final multiplication are outward, preserving this one-sided
confidence statement. The returned sample center is an indicator; `upper`
is a probability bound conditional on the sampling premises.

The law is algorithmic randomness, independent of physical/input
uncertainty. Threefry supplies replayable pseudorandom draws with unbiased
integer rejection; confidence assumes an ideal independent uniform bit
stream. A fixed seed is a realization, not a proof of PRNG independence.
Repeated adaptive requests need conditional independence and a total
failure budget (for example a union bound), not repeated use of one delta.

The primitive bounds a named vector norm. State/observable error additionally
needs stability and residual integration; a norm at one time proves neither.
Coarse envelopes can make this method useless for sparse or cancellation-
dominated residuals. No speed advantage is promised.

Source: [Hoeffding (1963), Theorem 2](https://doi.org/10.1080/01621459.1963.10500830).
The independent `sampled_norm` test checks constant, heterogeneous and sparse
vectors against binary128 norms over 2000 streams per family. At delta=0.05,
observed failure rates were 0, 0.0005 and 0; a five-sigma binomial margin
guards the empirical test. Empirical coverage supplements the theorem.

## Development evidence

GNU Release CMake/CTest focused gates pass 2/2, and the original KIN6D
consumer passes its independent evolution and derivative gates. The native
CMake route is authoritative. An attempted `fpm test --target
test_verified_matvec --profile release` links the pre-existing C smoke
`main` into the Fortran test, failing with duplicate `_main`; this is an
existing FPM mixed-test discovery limitation. Reproducer and next action:
repeat that command and isolate/register the C test in the owning FPM test
discovery path. No consumer workaround is installed.

Changes remain on local branch `error-strategy`, based on `e97a493`, with
no push. KIN6D uses an explicit FetchContent source override until promotion.

Chris&AI
