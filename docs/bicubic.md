# Fixed bicubic coefficient jet

`fortnum_bicubic` evaluates the fixed polynomial

\[
P(x,y)=\sum_{i=0}^3\sum_{j=0}^3 c_{i+1,j+1}x^i y^j.
\]

The first Fortran index is the power of `x`. Coordinates are displacements in
the coefficient basis; no cell lookup, normalization or grid-spacing factor is
applied. `jet(1:3)` is `[P, dP/dx, dP/dy]`.

```fortran
use fortnum_kinds, only: dp
use fortnum_bicubic, only: bicubic_jet, bicubic_jet_jvp, bicubic_jet_vjp
real(dp) :: c(4,4), x, y, jet(3), dc(4,4), dx, dy, djet(3)
real(dp) :: bar_jet(3), bar_c(4,4), bar_x, bar_y
! Initialize coefficients, coordinates, directions and cotangents first.
call bicubic_jet(c, x, y, jet)
call bicubic_jet_jvp(c, x, y, dc, dx, dy, jet, djet)
call bicubic_jet_vjp(c, x, y, bar_jet, jet, bar_c, bar_x, bar_y)
```

All sixteen coefficients and both coordinates are active. The analytical JVP
returns the derivative of the entire three-component jet, including the
coordinate dependence of both first partials. The analytical VJP uses the
ordinary Euclidean pairing and overwrites its input cotangents rather than
accumulating them. Inputs must not alias output arguments. The routines are
pure, have fixed scalar scratch, allocate no memory and own no mutable state.

One retained polynomial definition lives in
`tools/codegen/app/gen_bicubic.f90`. It proves the nested Horner value and
partials against the monomial definition and proves the contracted JVP/VJP
identity with FortSym's native polynomial engine. Ordinary builds use committed
`src/generated/fortnum_bicubic_kernel.f90`; no CAS is needed at runtime or
during an ordinary CMake/fpm build. The generator's tested revision is recorded
in `tools/codegen/fortsym-bicubic.lock`, independently of older generators.
Regenerate from `tools/codegen` with `fo exec gen_bicubic`; the existing
`check_generated.sh` includes the new generator and compares emitted bytes.

## Floating-point and enclosure limits

These APIs return floating-point centers, not public certificates. Finite
coefficients/coordinates alone do not prevent intermediate overflow or loss of
relative accuracy through cancellation. There is no clipping, hidden scaling,
finite-input rejection, or fallback: exceptional arithmetic follows the
caller's IEEE environment. Infinities/NaNs and overflowed intermediate values
must be rejected by a caller that requires a finite physical result. Gradual
underflow can lose the usual relative-error model; flush-to-zero/fast-math
builds are outside the following bound. Coordinate/cell reconstruction error
and coefficient uncertainty are separate from arithmetic error.

For round-to-nearest IEEE binary64, with every nonzero intermediate normal,
no overflow, and no unsupported reassociation, let `u=epsilon(1.0_dp)/2` and
`gamma_n=n*u/(1-n*u)`. For each exact jet component define its absolute
coefficient-weight sum:

\[
S_0=\sum |c_{ij}| |x|^{i-1}|y|^{j-1},\quad
S_x=\sum_{i\ge2}(i-1)|c_{ij}| |x|^{i-2}|y|^{j-1},\quad
S_y=\sum_{j\ge2}(j-1)|c_{ij}| |x|^{i-1}|y|^{j-2}.
\]

The emitted primal has 71 scalar additions/multiplications in total. Expanding
its arithmetic graph into coefficient monomials associates a product of at
most 71 rounding factors with each contribution, including shared temporaries;
the factors need not be independent. Triangle inequality therefore gives
`abs(computed_jet(k)-exact_jet(k)) <= gamma_71*S_k`. The regression uses the
more conservative `gamma_128*S_k`, computed with outward-rounded sums and
bounds, and an independent strict-interval monomial oracle. A small explicit
512-minimum-subnormal allowance covers the oracle's outward padding at exact
zeros. This allowance is test infrastructure, not a certificate for arbitrary
production underflow. The argument concerns the primal jet only; no public
rigorous JVP/VJP enclosure API is claimed.

The independent test covers all sixteen monomial positions, zero and signed
coordinates, mixed/cancelling coefficients, coefficient scales from `1e-100`
to `1e100`, independent first/second partials, all coefficient directions,
coordinate directions, a combined direction, central differences and VJP
pairing. These executable checks support the implementation; test agreement
does not establish an enclosure outside the stated arithmetic assumptions.

## Performance status

This is an analytical candidate for a concrete fixed-cell interpolation
consumer. No complete-workload speedup, mechanism winner or memory improvement
has been selected by this implementation. Whole-consumer CPU/wall/peak-memory
comparisons require unchanged physical gates and a controlled measurement
window; microkernel operation counts alone do not establish an improvement.
