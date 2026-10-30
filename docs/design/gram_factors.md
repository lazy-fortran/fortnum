# Gram and semidefinite square roots

Keep a covariance in factored form when rounding its entries would lose
positive semidefiniteness. `fortnum_gram_factor` compresses a rectangular factor
by applying LAPACK `dgeqrf` to its transpose. Transposing the resulting upper
trapezoid and padding with zero columns gives the same row Gram product to
floating-point accuracy. This avoids forming the covariance itself. The
Householder calculation uses LAPACK's scaled norm arithmetic, including for
subnormal factor entries.

`fortnum_psd_factor` initializes a factor from an existing matrix with LAPACK
`dpstrf`. Zero pivot tolerance preserves positive pivots without imposing a
global rank cutoff. After restoring the original row order, reconstruction
must satisfy the relative maximum-entry bound `128*n*epsilon`, where `n` is
the matrix dimension and `epsilon` is real64 machine precision. The same bound
checks symmetry. A nonpositive pivot therefore cannot by itself establish
validity. Factorization uses the original scale to preserve widely separated
positive diagonal entries; the reconstruction check scales the factors by the
square root of the maximum input magnitude. Success means a semidefinite
approximation within the stated bound. It does not certify exact rank or exact
semidefiniteness, and it introduces no jitter or spectral projection.

Independent binary128 sums of products check reconstruction after scale
normalization. Cases include rectangular and rank-deficient factors, an exactly
representable rank-one matrix, positive diagonals separated by 32 orders of
magnitude, covariance scales from `1e-300` to `1e300`, and factor entries near
`1e-320`. Representable subnormal spacing limits relative accuracy in the last
case, whose test tolerance is `2e-3`. Invalid dimensions, nonfinite values,
asymmetry and a strongly indefinite matrix must fail.

`bench_gram_factor` compares a reused workspace with a fresh one for a synthetic
matrix with seven rows and sixteen columns. Five CPU-time samples include the
fresh workspace's allocation and LAPACK query costs. Application performance
needs a separate measurement.
