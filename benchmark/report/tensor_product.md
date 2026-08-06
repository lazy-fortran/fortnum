# Tensor-product operator evidence

The new `fortnum_tensor_product` module represents

\[
  A_3 \otimes A_2 \otimes A_1
\]

without assembling the full matrix. Factor 1 is the innermost, fastest
varying grid dimension. A matrix-vector or multi-right-hand-side product
applies one factor to each mode in turn. The diagonal product uses the same
flattening convention.

The behavioral test uses independent hand-written Kronecker assembly for three
factors of sizes 2, 3, and 2. It checks vector, three-right-hand-side matrix,
and diagonal products, plus invalid-factor and invalid-lifecycle status cases.
The implementation is therefore not validated against its own contraction
routine.

The benchmark uses three tridiagonal factors, four right-hand sides, one pinned
CPU core, and float64 arithmetic. The small dense comparison assembles the
512-by-512 Kronecker matrix. The 4096-element structured case deliberately
does not assemble its 4096-by-4096 dense matrix. Each reported time is the
mean of the repeated product calls after one warm-up. Peak RSS is from
`/usr/bin/time -v` around the complete benchmark process. Code size is the
linked executable's `size` total. The raw rows are in
[`tensor_product.csv`](data/tensor_product.csv).

| compiler | N | tensor product (ms) | dense (ms) | speedup | peak RSS (KiB) |
| --- | ---: | ---: | ---: | ---: | ---: |
| gfortran `-O3` | 512 | 0.0270 | 0.4505 | 16.7x | 5968 |
| gfortran `-O3` | 4096 | 0.9908 | not measured | not measured | 5968 |
| nvfortran 26.5 `-O3` | 512 | 0.0366 | 1.8773 | 51.3x | 10188 |
| nvfortran 26.5 `-O3` | 4096 | 0.8049 | not measured | not measured | 10188 |

The build times for the four source files are 0.48 seconds with gfortran and
0.64 seconds with nvfortran. The linked code sizes are 31,716 and 36,843
bytes, respectively. This is CPU evidence for the structured contraction, not
a GPU-offload claim. The next accelerator slice must supply a device-resident
mode contraction and compare transfer-inclusive and resident timings.

The implementation is in `fortnum` and is intentionally independent of
`fortml`: regular-grid GP operators can consume it without moving numerical
primitives across repository boundaries.
