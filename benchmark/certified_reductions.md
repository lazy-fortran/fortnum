# Certified reduction experiment

This CPU benchmark compares the existing strict interval primitives with
experimental floating-point centers and aggregate error radii. The candidate
module is benchmark-local: it changes no production API or selection. A separate
row measures existing interval `exp` on point arguments in `[-2,0]`, representative
of bounded Gaussian-weight evaluation costs. No approximate `exp` is introduced.

## Enclosure argument

Inputs denote their exact binary64 values. Let `u=2^-53` and `eta=2^-1074`.
Under IEEE round-to-nearest with gradual underflow, each finite rounded addition
or multiplication admits

```
fl(z) = z*(1+delta) + xi,  |delta| <= u,  |xi| <= eta.
```

The absolute underflow allowance is conservative (one smallest subnormal).
For a controlled serial sum of `n` inputs use `k=n`; for a serial dot product
with separately rounded products use `k=2n`. A term in either computation has
at most `k` rounding factors. For `k*u<1`, their deviation from one is at most
`gamma_k=k*u/(1-k*u)`, and their magnitude is at most `1/(1-k*u)`.
There are at most `k` absolute underflow errors. Consequently, with exact
absolute sum `A=sum(abs(x))`, or `A=sum(abs(x*y))` for the dot product,

```
|center - exact_result| <= gamma_k*A + b,
b = k*eta/(1-k*u).
```

The candidate accumulates a second ordinary floating-point sum `m` of input
magnitudes, or magnitudes of the same rounded products. This has the same error
bound, so `m >= (1-gamma_k)*A-b`. The implementation uses outward arithmetic to
compute `g >= gamma_k`, `B >= b`, and

```
A_upper = (m+B)/(1-g),
radius = g*A_upper+B.
```

`B` is evaluated as `k*eta*(1+g)`. The denominator is rounded downward; all
other budget operations are rounded upward. Require `g<1` and finite center,
magnitude and radius. This also bounds complete cancellation and products that
underflow to zero. Empty reductions return exactly zero with radius zero.
The proof uses the documented `fortnum_rounding` enclosure contract; it is
not a newly machine-checked theorem.

## Validity and failure

- Preserve the declared serial operation order: no unsafe reassociation,
  fast-math, flush-to-zero or denormals-are-zero. GNU builds disable fast-math
  and FMA contraction for the benchmark core. Compiler flags remain part of
  the experiment's evidence.
- The executable checks binary64 metadata, round-to-nearest and basic gradual
  underflow behavior before running. These checks cannot prove compiler
  correctness or detect every unsafe transformation.
- `status=0` encloses the exact sum/dot in the mathematical set
  `[center-radius,center+radius]`. To construct binary64 endpoints, round their
  subtraction/addition outward; an ordinary rounded endpoint is insufficient.
- `status=1` means invalid input (shape mismatch or nonfinite input).
  `status=2` means range/budget failure. Failed results carry infinite radius
  and must not be accepted as finite certificates. Intermediate overflow is
  rejected, including a magnitude sum overflowing despite an exact canceled
  result. The strict baseline also reports unbounded endpoints as range failure.
- A caller may retry with strict intervals or higher precision. There is no
  automatic fallback hidden in the timing or an assertion that fallback always
  produces a finite enclosure.

Independent fixtures use exact dyadic identities, integer-accumulated workload
oracles, lost-unit cancellation around `2^53`, smallest-subnormal additions,
rounded-to-zero products, and overflow/nonfinite/shape/environment rejection.
Their exact expected values are represented in binary128 without evaluating the
tested reduction to construct the reference. Existing general transcendental
tests remain the enclosure oracle for the interval `exp` timing row.

## Timing protocol

Configure the existing standalone benchmark build in Release mode, then build
only `bench_certified_reductions`. Run the registered validation before timing:

```sh
cmake -S benchmark -B build-certified-bench -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-certified-bench --target bench_certified_reductions -j 1
ctest --test-dir build-certified-bench -R '^certified_reductions_validation$' -V
build-certified-bench/bin/bench_certified_reductions aggregate sum positive 4096 100
build-certified-bench/bin/bench_certified_reductions interval sum positive 4096 100
build-certified-bench/bin/bench_certified_reductions aggregate dot cancel 65536 20
build-certified-bench/bin/bench_certified_reductions interval dot cancel 65536 20
build-certified-bench/bin/bench_certified_reductions interval exp negative 4096 10
```

Each process validates before timing, uses three warmups and emits fifteen raw
wall-clock samples as CSV with an observable sink. Post-sample validation is
outside timing. Candidates return the same center/radius semantics, but their
radii can differ; compare both time and enclosure width. Sum and dot rows
include input scanning and budget construction. The `exp` row includes point
interval construction and stores all output intervals; it is not a reduction
or an end-to-end physics solver comparison. Record source/patch identity,
compiler, flags, host, affinity, lengths and repetitions with raw output.

Caller tolerances apply after stability or sensitivity amplification. These
arithmetic enclosures alone do not certify a continuous equation, quadrature
remainder, unresolved tail, observable accuracy or derivative of an implicit
solution. No production speedup or universal hardware/backend support is claimed.
