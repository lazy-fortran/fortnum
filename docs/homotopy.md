# Polynomial homotopy continuation

`fortnum_homotopy` finds the isolated solutions of square polynomial systems
in n complex unknowns. It fills the solution set for random coefficients by
monodromy, then moves that set to the coefficients of interest along a
straight parameter homotopy. Tracking runs in double precision and stays
accurate for solutions whose coordinates span many orders of magnitude.

## Systems

Equation i of F(x; c) = 0 is a sum of terms, and each term multiplies one
entry of the coefficient vector c by a monomial:

    F_i(x; c) = sum over terms k of equation i of c(coef(k)) * prod_j x(vars_k(j))

The coefficient vector is the homotopy parameter: the family is "same
supports, arbitrary coefficients". A term usually owns its coefficient index,
but indices may be shared.

```fortran
use fortnum_kinds, only: dp
use fortnum_homotopy
type(poly_system_t) :: sys
! x^2 + y^2 - 5 = 0, x y - 2 = 0 with coefficients c(1..5)
call poly_system_init(sys, 2)
call poly_system_add_term(sys, 1, [1, 1], 1)
call poly_system_add_term(sys, 1, [2, 2], 2)
call poly_system_add_term(sys, 1, [integer ::], 3)
call poly_system_add_term(sys, 2, [1, 2], 4)
call poly_system_add_term(sys, 2, [integer ::], 5)
call poly_system_finalize(sys)
```

`vars` lists variable indices with multiplicity, so `[1, 1, 2]` is x1^2 x2 and
an empty list is a constant. A term may contain at most 32 distinct
variables. `sys%ncoef` is the largest coefficient index used. Every equation
needs at least one term.

`poly_system_eval(sys, x, c, dc, f, jac, ft, rowscale)` returns F(x; c), the
dense Jacobian dF/dx, F(x; dc) (the t-derivative along c(t) = ca + t dc) and,
for each equation, the largest term magnitude |c_k m_k(x)|. It is pure and
does not allocate.

## Solving

```fortran
type(homotopy_options_t) :: opts
type(fortnum_status_t) :: status
complex(dp) :: x0(2), c0(5), ct(5)
complex(dp), allocatable :: sols(:, :), out(:, :)
logical, allocatable :: ok(:)
real(dp), allocatable :: acc(:)
integer :: nsol

call homotopy_plant(sys, x0, c0, seed=1)
opts%target = 4                 ! 0 when the generic count is unknown
call homotopy_monodromy(sys, x0, c0, opts, sols, nsol, status)
ct = [1, 1, -5, 1, -2]
allocate (out(2, nsol), ok(nsol), acc(nsol))
call homotopy_move(sys, sols, nsol, c0, ct, opts, out, ok, acc)
```

- `homotopy_plant(sys, x0, c0, seed)` draws x0 on the unit torus and complex
  Gaussian coefficients, then solves the coefficient of the first term of
  each equation so that F(x0; c0) = 0. That coefficient index must not occur
  in another equation.
- `homotopy_monodromy(sys, x0, c0, opts, sols, nsol, status, stats)` starts
  from one solution x0 of F(x; c0) = 0. Each loop draws fresh random c1, c2
  and tracks every known solution along c0 -> c1 -> c2 -> c0 in parallel.
  New endpoints are added. It stops when `opts%target` solutions are known or
  after `opts%stagnant_loops` loops without a new one, also when the target is
  known, so set `stagnant_loops` high in that case. `sols(n, nsol)` is
  allocated on return.
- `homotopy_move(sys, sols, nsol, ca, cb, opts, out, ok, acc, stats)` tracks
  every column of `sols` from ca to cb in parallel. `ok(i)` reports whether
  path i reached cb; `acc(i)` is its final relative Newton correction.
- `homotopy_track(sys, x, ca, cb, opts, ok, nsteps, acc)` tracks one solution
  along the segment.

`status` from `homotopy_monodromy` is `FORTNUM_OK` when the target count was
reached exactly or, with an unknown target, when the search stagnated. It is
`FORTNUM_CONVERGENCE_ERROR` when the target was missed or exceeded or
`max_loops` ran out while new solutions still appeared, and
`FORTNUM_DOMAIN_ERROR` when x0 is not a regular solution. The optional
`homotopy_stats_t` counts loops, segment tracks, failed paths, predictor-
corrector steps, LU factorizations and path jumps.

Monodromy and `homotopy_move` give identical results for every thread count:
paths are independent and deduplication is serial. Random draws come from
`fortnum_rng` streams derived from `opts%seed`.

## Method

Every linear solve is equilibrated. Rows are scaled by the largest term
magnitude of their equation and columns by the current |x_j|, so corrections
are relative, and convergence is measured in max_j |dx_j| / |x_j|. This keeps
solutions whose coordinates span many orders of magnitude, with unscaled
Jacobian condition numbers far beyond 1/eps, trackable in double precision.

The predictor is a coordinatewise Pade [7/1] approximant of the order-8
Taylor polynomial of the path at the accepted point. The Taylor coefficients
come from truncated power-series evaluation of the terms, which keeps partial
products per term so that each order costs O(j) per factor, and from back
substitutions with the LU factors of the last corrector iteration. When the
Pade pole lies within twice the step, the Taylor polynomial is used. A step
therefore costs about two factorizations, one per Newton iteration.

The corrector is Newton with a hand-written column-oriented LU with partial
pivoting. A step is accepted only when the first relative correction is below
`max_first_correction` (1e-3), which keeps path jumping rare, and Newton then
either converges to `tol`, predicts its next correction below `tol`, or stops
contracting once the previous correction was below 1e-7 (roundoff). The step
size aims at a first correction of `step_target`, assuming an error of order
h^9. A path fails when the step falls below `h_min` times the distance in t
to the nearer end of the segment, or after `max_steps` steps. So a start fibre
at a very different coefficient scale can be left with tiny steps, and an
endpoint can be approached over many orders of magnitude: the tracker keeps
t in the first half of a segment and the remaining distance 1 - t in the
second, and forms the coefficients from the nearer end. At the end of each
segment Newton refines the point to `final_tol`, and the last relative
correction is reported as its accuracy.

Near an ill-conditioned point the residual formed in double precision is
dominated by roundoff, so Newton stalls above roundoff level or, worse, sees a
residual that has rounded to zero. When Newton stalls, the tracker therefore
forms every further residual of that path with all products and sums in quad
precision (`real128`) and keeps the double Jacobian and its LU factors. This
is iterative refinement, so the point reaches double accuracy as long as the
scaled condition number stays well below 1/epsilon. The final refinement of
every segment uses the quad-precision residual, so the reported accuracy is
not an artefact of a residual that rounded to zero.

Where the scaled condition number approaches 1/epsilon, even the double Jacobian
and LU factors are too inaccurate. `homotopy_move` then resumes a path that
failed in double precision from its last accepted point with a tracker whose
every operation is in quad precision (RK4 predictor, scaled Newton corrector);
`quad_fallback = .false.` turns this off. The fallback is slower but runs only
on the paths that need it.

A straight coefficient segment maps the start solutions one to one onto the
end solutions, so two coinciding endpoints mean that a path jumped to another
branch. `homotopy_move` detects coinciding endpoints by sorting (as the
monodromy deduplication does) and tracks every member of such a group again
with a step target 100 times smaller, for up to three rounds. Endpoints that
still coincide are counted in `stats%jumps`.

Monodromy accepts a loop endpoint only if its final correction is below
`dup_tol`. Two points are one solution when every coordinate agrees to the
relative distance `dup_tol`. Deduplication sorts by a real key, a fixed
random combination of log|x_j| and the phases x_j/|x_j|, which moves by at
most a known multiple of `dup_tol` for equal solutions, and compares
coordinates only within that key window. Its cost is O(N log N) per loop.
Within one loop the endpoints must be distinct, because a loop permutes the
fibre; coinciding endpoints are counted as path jumps in `stats%jumps`.

## Limits

- Solutions are computed in double precision only.
- Monodromy finds the solutions in the complex torus of the generic member of
  the family. Solutions with a coordinate fixed at zero for all coefficients
  are not reached.
- Singular target solutions (multiple roots, roots at infinity, or exact
  zeros of a coordinate, where relative coordinates break down) end with a
  large `acc` or `ok = .false.`; there is no endgame.
