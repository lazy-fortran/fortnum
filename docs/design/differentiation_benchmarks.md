# Differentiation evidence catalog

Status: index of committed machine-readable records.

Exact measurements live under `benchmark/reference/`. This document identifies
their scope, oracle, and selection question. It does not duplicate append-only
run logs.

## Comparison contract

Two rows support a speed comparison only when they return the same
caller-visible values for the same:

- mathematical operator
- derivative product
- active arguments
- workload dimensions
- reusable primal state
- target and residency

“Faster” means lower median wall clock for that complete workload. Dispersion
and the selection tie threshold remain part of the decision. Peak memory is
candidate-specific.

Finite differences are labeled diagnostics. Their time includes the perturbed
primal evaluations required by the stated oracle.

## CPU record families

| Family | Record pattern | Products and comparisons | Independent validation |
| --- | --- | --- | --- |
| generated algebra | `ryzen9_5950x_dawson_generated_family.json`, determinant and inverse records | fused/separate, generated/diagnostic, JVP/VJP | formulas, finite differences, matrix identities, adjoint identities |
| special functions | `ryzen9_5950x_bessel_*.json`, gamma, erf, and hypergeometric tournament records | analytical, autodiff, hybrid across numerical regions | complete-objective finite differences and custom-rule provenance |
| Enzyme infrastructure | `ryzen9_5950x_enzyme_*.json`, `*_fixture_migration.json` | wrapper ABI, shared support, code size, migration regression | real Enzyme formulas, adjoint identities, negative repository guards |
| interpolation | `ryzen9_5950x_lagrange_*.json`, `ryzen9_5950x_bspline_*.json` | active points, values, nodes, knots, coefficients, fused products | independent cubic formulas, finite differences, adjoint identities, crossing status |
| quadrature and integration | `ryzen9_5950x_integrate_*.json`, fixed-quadrature records | fixed, moving-bound, frozen-trace, autodiff, hybrid, batched | exact integrals, Leibniz terms, trace replay, finite differences |
| roots and fixed points | `ryzen9_5950x_scalar_root_*.json`, `ryzen9_5950x_vector_root_*.json`, `ryzen9_5950x_multiroot_*.json`, fixed-point records | implicit, hybrid, iteration autodiff, diagnostic, factor reuse | closed-form roots, complete-solve finite differences, residual equations, adjoint identities |
| linear algebra | `ryzen9_5950x_linear_*.json`, direct/iterative solver records, LU record | factor reuse, multiple RHS, forward/reverse Enzyme, diagnostic | solve residuals, matrix products, finite differences, adjoint identities |
| ODE | `ryzen9_5950x_ode_*.json` | continuous tangent, discrete tangent/adjoint, parameters, checkpoint/recompute, events, implicit stages | closed forms, refinement, frozen-map finite differences, matrix exponentials, adjoint identities |
| selection | build and static-selection records | deterministic registry and CMake consumption | exact workload lookup and negative malformed-record tests |

CPU records use separate processes when peak RSS must be attributed to one
candidate. Low-latency reference runs record CPU affinity. Cache and work
counters appear when they explain a wall-clock difference.

## GPU record families

| Family | Record pattern | Selection question |
| --- | --- | --- |
| Dawson | `rtx5060ti_dawson_*.json` | algebraic form, fusion, residency, JVP/VJP |
| multi-input scalar | `rtx5060ti_multi_input_*.json` | input/direction scaling, layout, profiling |
| backend selection | `rtx5060ti_gpu_backend_selection.json` | OpenACC versus OpenMP target |
| interpolation | `rtx5060ti_lagrange4_products.json` | batch crossover and returned-adjoint transfer |
| quadrature | `rtx5060ti_fixed_quadrature_products.json` | arbitrary-order contraction |
| linear algebra | `rtx5060ti_linalg3_products.json` | determinant and inverse products |
| FFT | `rtx5060ti_fft8_products.json` | transform and adjoint transform |
| implicit root | `rtx5060ti_implicit_root_jvp.json` | generated residual plus analytical solve |
| ODE | `rtx5060ti_ode_trace2_products.json`, `rtx5060ti_ode_terminal_objective.json` | trace products and complete value-plus-gradient application |

GPU records separate resident from transfer-inclusive time and prove
non-host execution. They include peak device allocation and released profiler
counters where available.

## Scaling evidence

The committed records cover:

- 1, 4, and 16 JVP directions
- 1, 4, and 16 VJP cotangents
- 2, 4, 8, and 16 active scalar inputs
- scalar and vector outputs
- fixed and increasing linear-system right-hand sides
- ODE parameter, trajectory, and checkpoint dimensions
- CPU cache behavior
- launch-bound, transition, and large GPU batches
- resident and transfer-inclusive GPU execution

These samples establish local crossover evidence. They do not define a
universal forward/reverse threshold.

The gamma tournament keeps the shape fixed and activates the integration
limit. Enzyme 22 crashes in type analysis for active shape differentiation
through Flang's `log_gamma` lowering, so that unreliable candidate is neither
selected nor normalized. Full shape derivatives remain available through the
independently tested analytical products.

The hypergeometric tournament activates real `z` at fixed real parameters.
Because Enzyme 22 cannot differentiate Flang's complex `cexp`, its candidate
uses a real specialization of the same Kummer, Taylor, and asymptotic
algorithms. Primal equivalence to the production complex path is a mandatory
validation gate.

## Current normalized set

`benchmark/report/data/mechanism_tournaments.csv` includes only tournaments
with competing mechanism timings and a selection decision. Each row points to
its source JSON.

The normalized set excludes:

- comparisons solely between variants of one mechanism
- reliability-only experiments
- incomplete candidate sets
- primal-only performance records

See [differentiation_report.md](differentiation_report.md) for aggregate
statistics and figures.

## Reproduction

Benchmark commands and record requirements are in
[benchmark/README.md](../../benchmark/README.md). The optional CPU Enzyme path
is described in [enzyme_toolchain.md](enzyme_toolchain.md). GPU gates are in
[gpu.md](gpu.md).
