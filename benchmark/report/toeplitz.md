# FFT Toeplitz operator evidence

The `fortnum_toeplitz` module represents a one-dimensional Toeplitz grid
operator from its first column and optional first row. It caches a circulant
embedding spectrum at initialization, then applies matrix-free FFT products.
The symmetric benchmark case uses the same exponential covariance column for
both the lower and upper diagonals.

The behavioral test assembles an asymmetric 5-by-5 Toeplitz matrix by hand and
checks vector, multi-right-hand-side, diagonal, symmetric-default-row, and
invalid-shape behavior. This dense assembly is an independent oracle; it does
not call the implementation's FFT or embedding routines.

The benchmark uses float64 arithmetic, four right-hand sides, one pinned CPU
core, one warm-up product, and repeated products for each size. Dense timing is
reported only at 512 samples, where the explicit matrix is small enough for a
matched comparison. The 128-to-4096 rows expose scaling without constructing
the larger dense matrices. Peak RSS is from `/usr/bin/time -v` around the full
benchmark process, build time is a direct clean compile of the benchmark
sources, and code size is the linked executable's `size` total. Raw data are in
[`toeplitz.csv`](data/toeplitz.csv).

| compiler | samples | Toeplitz (ms) | dense (ms) | speedup |
| --- | ---: | ---: | ---: | ---: |
| gfortran `-O3` | 512 | 0.0400 | 0.4515 | 11.29x |
| nvfortran 26.5 `-O3` | 512 | 0.0352 | 1.8576 | 52.71x |

The measured 128-to-4096 time ratio is 51.3x for gfortran and 54.1x for
nvfortran while the sample count grows by 32x, consistent with an
approximately (O(n\log n)) product rather than the (O(n^2)) dense path.
The gfortran and nvfortran peak RSS values are 6,636 KiB and 10,264 KiB;
build times are 1.00 s and 1.52 s; linked sizes are 71,423 and 92,988 bytes.
These are host-resident CPU measurements. The GeForce RTX 5060 Ti is recorded
as available hardware, not as an accelerator result: the current FFT backend
does not keep the embedding spectrum and work arrays device-resident.

The implementation is in `fortnum` and is a reusable primitive for structured
GP consumers. Multilevel tensor-grid embeddings and accelerator-resident FFT
products remain open roadmap items.
