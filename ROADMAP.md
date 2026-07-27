# Differentiation implementation plan

This file tracks the remaining work for differentiable programming in
`fortnum`. Public terminology is defined in `docs/design/ad.md`: `autodiff`,
`analytical`, and `hybrid`.

## Hard execution rules

Work through this file one checkbox at a time.

1. Select exactly one unchecked item.
2. Keep the implementation as small as possible while fully satisfying that
   item.
3. Avoid overengineering. Prefer simple, direct, effective code over general
   frameworks, speculative abstractions, and unused flexibility.
4. Add an independent behavioral oracle and run the relevant test matrices.
5. Measure and report performance for every item. Report at least runtime and
   peak memory, with before and after numbers or competing-candidate numbers.
   Record the workload, compiler, hardware, and measurement method. Wall-clock
   time for the complete workload is the primary selection metric.
6. Measure scaling with active-input count, output count, and derivative
   directions wherever those dimensions can change the forward-versus-reverse
   verdict. Report CPU work, peak memory, and cache behavior in addition to
   wall-clock time when the host exposes reliable counters.
7. Update committed benchmark evidence when the item affects differentiation
   performance or candidate selection.
8. Check off the completed item in this file.
9. Commit the implementation, tests, evidence, and checklist update together.
10. Push the commit to `main`.
11. Only after the push succeeds, select the next unchecked item.

Do not combine several checklist items into one implementation cycle. If an item
turns out to contain several independent changes, split it into smaller
checkboxes before writing code.

A mechanism name never selects a winner. Production selection requires
independent validation plus measured application runtime and peak memory.

## Completed foundation

- [x] Define `autodiff`, `analytical`, and `hybrid` as candidate mechanisms.
- [x] Implement derivative provenance and quality merging.
- [x] Implement deterministic candidate metadata and pre-loop selection.
- [x] Build a real Flang and Enzyme custom-rule path.
- [x] Generate and consume a fused Dawson value/JVP kernel with `fortsym`.
- [x] Benchmark analytical, autodiff, and hybrid Dawson candidates.
- [x] Implement analytical implicit linear-solve JVP and VJP products.
- [x] Implement combined-active Lagrange and B-spline JVP/VJP products.

## Core infrastructure

- [x] Add factorization reuse to the analytical linear-solve JVP.
- [x] Add transpose-factorization reuse to the analytical linear-solve VJP.
- [x] Add reusable preconditioner hooks for implicit products.
- [x] Connect committed benchmark records to build-time candidate selection.
- [x] Generate a static per-workload selection registry.
- [x] Add conditioning diagnostics to implicit derivative products.
- [x] Add derivative-reliability status for ill-conditioned implicit products.
- [x] Measure candidate-specific peak memory instead of shared-process RSS.
- [x] Audit the linear-solve VJP benchmark target and refresh its committed
  Release and candidate-memory evidence.

## Roots and residual equations

- [x] Add a generic analytical implicit tangent boundary for scalar roots.
- [x] Add a generic analytical implicit adjoint boundary for scalar roots.
- [x] Add a generic analytical implicit tangent boundary for vector roots.
- [x] Add a generic analytical implicit adjoint boundary for vector roots.
- [x] Add analytical implicit tangent products for fixed points.
- [x] Add analytical implicit adjoint products for fixed points.
- [x] Add a hybrid scalar-root JVP using forward-mode autodiff residual products.
- [x] Add a hybrid scalar-root VJP using reverse-mode autodiff residual products.
- [x] Add a hybrid vector-root JVP using forward-mode autodiff residual products.
- [x] Add a hybrid vector-root VJP using reverse-mode autodiff residual products.
- [x] Return and reuse converged Jacobians from analytical vector-root solves.
- [x] Reuse converged root Jacobian factorizations for JVPs.
- [x] Reuse converged transposed root Jacobian factorizations for VJPs.
- [x] Report JVP reliability near singular root Jacobians.
- [x] Report VJP reliability near singular root Jacobians.
- [x] Benchmark analytical implicit, hybrid, autodiff-through-iterations, and
  finite-difference scalar-root candidates.
- [x] Benchmark analytical implicit, hybrid, autodiff-through-iterations, and
  finite-difference vector-root candidates.

## Integration and quadrature

- [x] Add analytical differentiation under a fixed-bound integral.
- [x] Add analytical moving-lower-bound terms.
- [x] Add analytical moving-upper-bound terms.
- [x] Compose autodiff integrand JVPs with analytical fixed quadrature.
- [x] Compose autodiff integrand VJPs with analytical fixed quadrature.
- [x] Add an analytical frozen-trace adaptive-integration candidate.
- [x] Add an autodiff frozen-trace adaptive-integration candidate.
- [x] Add a hybrid adaptive-integration candidate.
- [x] Benchmark fixed quadrature candidates.
- [x] Benchmark smooth adaptive-integration candidates.
- [x] Benchmark singular adaptive-integration candidates.
- [x] Benchmark batched integration candidates.

## ODE solvers and events

- [x] Add hybrid forward sensitivities using autodiff RHS JVPs.
- [x] Add an analytical discrete adjoint for one explicit ODE method.
- [ ] Accumulate parameter VJPs in the discrete adjoint.
- [ ] Add checkpointed reverse differentiation as a candidate.
- [ ] Add recomputation-based reverse differentiation as a candidate.
- [ ] Add analytical implicit-stage tangent products.
- [ ] Add analytical implicit-stage adjoint products.
- [ ] Differentiate event times through the event residual equation.
- [ ] Document and test the continuous sensitivity contract.
- [ ] Document and test the discrete sensitivity contract.
- [ ] Benchmark a short nonstiff trajectory.
- [ ] Benchmark a long nonstiff trajectory.
- [ ] Benchmark a stiff trajectory.
- [ ] Benchmark a many-parameter trajectory.
- [ ] Benchmark an event-driven trajectory.

## Linear algebra

- [ ] Add an analytical determinant JVP.
- [ ] Add an analytical determinant VJP.
- [ ] Add an analytical inverse JVP for callers that require an inverse.
- [ ] Add an analytical inverse VJP for callers that require an inverse.
- [ ] Add a reusable LU factorization object.
- [ ] Add multiple-right-hand-side tangent solves.
- [ ] Add multiple-right-hand-side adjoint solves.
- [ ] Add a forward-mode autodiff direct-solver JVP comparator.
- [ ] Add a reverse-mode autodiff direct-solver VJP comparator.
- [ ] Add an autodiff iterative-solver comparator where applicable.
- [ ] Add a hybrid BLAS or LAPACK custom-rule candidate where measurement
  justifies it.
- [ ] Benchmark all implemented linear-algebra candidates.

## Interpolation and splines

- [ ] Add analytical products for active support-node locations.
- [ ] Add analytical products for active knot locations where defined.
- [ ] Add analytical implicit differentiation of fitted spline coefficients.
- [ ] Define derivative status at interpolation-cell crossings.
- [ ] Define derivative status at knot crossings.
- [ ] Benchmark separate and fused combined-active Lagrange products.
- [ ] Benchmark separate and fused combined-active B-spline products.
- [ ] Compare analytical basis products with fixed-span autodiff candidates.

## Special functions

- [ ] Run an analytical, autodiff, and hybrid Bessel tournament.
- [ ] Run an analytical, autodiff, and hybrid gamma-function tournament.
- [ ] Run an analytical, autodiff, and hybrid error-function tournament.
- [ ] Run an analytical, autodiff, and hybrid hypergeometric tournament.
- [ ] Generate region-specific analytical candidates with `fortsym`.
- [ ] Benchmark fused and separate value/JVP special-function candidates.
- [ ] Benchmark representative small-argument regions.
- [ ] Benchmark representative transition regions.
- [ ] Benchmark representative asymptotic regions.

## FFT and transforms

- [ ] Benchmark the analytical FFT JVP against autodiff of the implementation.
- [ ] Benchmark the analytical FFT VJP against autodiff of the implementation.
- [ ] Verify complex-adjoint and normalization conventions independently.
- [ ] Add a hybrid custom rule for an external FFT library.
- [ ] Benchmark scalar FFT workloads.
- [ ] Benchmark batched FFT workloads.
- [ ] Benchmark an FFT-based spectral objective.

## `fortsym`

- [ ] Normalize generated algebraic signs and constant products.
- [ ] Generate contracted JVP expressions directly.
- [ ] Generate contracted VJP expressions directly.
- [ ] Generate and retain several algebraic variants.
- [ ] Rank generated variants by structural operation count before compilation.
- [ ] Add stability-preserving transformations.
- [ ] Add target-aware factoring only when benchmarks justify it.
- [ ] Record generated native code size.
- [ ] Add generated-source drift checking to CI.
- [ ] Record the exact `fortsym` revision for every generated artifact.

## Module tournaments

- [ ] Complete the special-functions tournament and commit its evidence.
- [ ] Complete the interpolation tournament and commit its evidence.
- [ ] Complete the FFT tournament and commit its evidence.
- [ ] Complete the quadrature tournament and commit its evidence.
- [ ] Complete the roots tournament and commit its evidence.
- [ ] Complete the linear-algebra tournament and commit its evidence.
- [ ] Complete the ODE tournament and commit its evidence.

Every tournament record must include:

- [ ] Independent validation error.
- [ ] Median runtime and dispersion.
- [ ] Candidate-specific peak memory.
- [ ] Generated or native code size where relevant.
- [ ] Hardware and compiler identity.
- [ ] Reusable primal state.
- [ ] The selected candidate and deterministic selection rationale.
- [ ] Scaling over representative active-input and output counts.
- [ ] Forward-mode versus reverse-mode crossover evidence.
- [ ] CPU work and cache counters where supported.

## Application-level benchmarks

- [ ] Benchmark a root-constrained scalar objective.
- [ ] Benchmark a parameterized integral.
- [ ] Benchmark spline fitting followed by a downstream gradient.
- [ ] Benchmark an FFT-based spectral objective.
- [ ] Benchmark nonlinear least squares.
- [ ] Benchmark stiff ODE sensitivity.
- [ ] Benchmark nonstiff ODE sensitivity.
- [ ] Benchmark an event-driven trajectory.
- [ ] Benchmark a small PDE or residual application.
- [ ] Benchmark a representative `itpplasma` workload.
- [ ] Measure end-to-end forward/reverse scaling over active-input counts.
- [ ] Measure end-to-end forward/reverse scaling over output counts.
- [ ] Record wall-clock, peak-memory, and cache-performance crossover curves.

## Second order

- [ ] Identify a concrete consumer that requires an HVP.
- [ ] Add an autodiff forward-over-reverse HVP candidate.
- [ ] Add an autodiff reverse-over-forward HVP candidate.
- [ ] Add an analytical or hybrid implicit HVP where justified.
- [ ] Validate HVPs with an independent second-directional or high-precision
  oracle.
- [ ] Measure HVP runtime and peak memory.
- [ ] Add a full Hessian only when a measured consumer requires one.

## Repository quality

- [ ] Resolve or explicitly baseline the `fo lint` array-temporary warnings.
- [ ] Repair the fpm and `fo test` oracle-argument mismatch.
- [ ] Add CI coverage for the real Flang and Enzyme hybrid pipeline.
- [ ] Add scheduled benchmark runs on the reference Ryzen 9 5950X.
- [ ] Re-run candidate selection after material compiler changes.
- [ ] Re-run candidate selection after material Enzyme changes.
- [ ] Re-run candidate selection after material `fortsym` changes.
- [ ] Re-run candidate selection after material primal-code changes.
- [ ] Re-run candidate selection on new target hardware.
