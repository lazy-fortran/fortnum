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
| FFT | `ryzen9_5950x_fft8_{jvp,vjp}_tournament.json`, `ryzen9_5950x_fftw8_custom_rule.json`, `ryzen9_5950x_fft_scalar_products.json` | analytical transform products versus Enzyme, external-library custom rule, scalar scaling | direct DFT, direct adjoint DFT, dot identity, finite difference, and NumPy transform oracle |
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

For the generated Dawson outer objective, one fused value/JVP call takes
9.67 ns versus 15.28 ns for separate value and JVP calls. The fused value/VJP
takes 9.75 ns versus 15.81 ns. Fusion is therefore 1.58× and 1.62× faster,
respectively, and executes 22.4% fewer instructions without a cache penalty.
All four scalar candidates use fixed storage; their roughly 2.8 to 2.9 MB process
RSS does not resolve derivative workspace differences.

Root and fixed-point records support the same operator-level conclusion.
For the two-state fixed-point fixture, the analytical implicit JVP takes
97.80 ns versus 1,224.40 ns for two complete perturbed solves; its VJP takes
132.40 ns versus 2,495.14 ns for four solves. The speedups are 12.52× and
18.85×, with maximum absolute errors of \(1.18\times10^{-11}\) and
\(1.01\times10^{-10}\). Reusing a converged 16-state root factorization makes
the JVP 2.96× faster and the VJP 2.98× faster than refactorization. The fixed
point and factor-reuse choices change candidate-process RSS by at most 184,320
bytes, so wall clock decides these workloads. Scalar and two-state root
tournaments compare analytical, hybrid, fixed-iteration autodiff, and
diagnostic implementations separately for JVP and VJP; no mechanism is assumed
to win both products.

The linear-algebra tournament covers 12 complete workloads. Analytical
candidates win all 12 on this CPU, but for different reasons. Generated
contracted determinant and inverse products are 1.28× to 8.19× faster than
complete finite-difference diagnostics. Reusing a 16-by-16 LU factorization is
2.82× faster for a JVP and reusing its transpose factorization is 2.93× faster
for a VJP; current complete-solve finite-difference errors are
\(3.72\times10^{-11}\) and \(4.45\times10^{-12}\). Repeated scalar products
beat the available batched interfaces by 1.76× for 16 JVP directions and 1.08×
for 16 VJP cotangents. Direct-solver analytical products beat Enzyme by 3.81×
for 16 JVP directions and 1.53× for 16 VJP cotangents. These are local
wall-clock results, not a mechanism-wide rule.

ODE choices change with the requested product and scale. On the short
frozen-trace fixture, a hybrid generated Enzyme RHS JVP is 1.01× faster than
the analytical RHS JVP in the Enzyme build; the normal build selects the
analytical implementation. For a full 16-by-16 initial-state Jacobian,
analytical forward sweeps are 1.23× faster than analytical reverse sweeps.
For a scalar objective, analytical reverse is 1.74× faster on the short
trajectory and 1.71× faster on the 400-step trajectory. Parameter scaling
crosses over between one and four parameters: forward wins at one, while
reverse wins by 2.54× at four and 8.51× at sixteen. Implicit event products
reach 31.37× over repeated event relocation. Checkpointing reduces retained
trace bytes but does not reduce complete-workload peak RSS in the measured
implementation, so full-trace reverse remains faster. The stiff record covers
differentiation around an explicit primal; it does not claim a production
stiff-solver selection.

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

The fixed length-eight batch record covers 256, 65,536, and 1,048,576
transforms. At the largest batch, resident GPU execution is 26.8 times faster
than CPU for JVP and 29.6 times faster for VJP. Transfer-inclusive GPU
execution remains slower than CPU for both products. The shared CPU/GPU batch
leaf has a current CPU CTest gate against direct DFT and adjoint oracles.

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

Small-argument measurements do not admit one mechanism-wide winner. At 16
products, analytical Bessel takes 818 ns for both JVP and VJP, while the
closest hybrid JVP takes 825 ns and autodiff VJP takes 3,202 ns. Gamma selects
autodiff for its 1,524 ns JVP but analytical for its 2,654 ns VJP. Erf is a
practical tie: raw autodiff takes 298 ns for JVP and 305 ns for VJP, versus
313 and 307 ns analytically; the normal-build, device-portable analytical
candidate remains selected. The corresponding record includes peak RSS and
cache counters.

Transition regions remain product-specific. At 16 products, Bessel selects a
3,531 ns autodiff JVP and a 5,566 ns analytical VJP. Gamma similarly selects a
2,819 ns autodiff JVP and a 4,245 ns analytical VJP. Erf is again a practical
tie: autodiff is raw-fastest for JVP, while analytical wins VJP and remains
the portable normal-build selection. The cache counters track the same
wall-clock verdicts: the clear winners execute fewer instructions and incur
fewer misses.

Tail and asymptotic measurements reach the same conclusion. Bessel selects a
2,212 ns autodiff JVP and a 2,798 ns analytical VJP; gamma selects a 2,820 ns
autodiff JVP and a 5,710 ns analytical VJP. Erf remains a practical tie.
Hypergeometric is the important exception enabled by code generation: its
fortsym-generated asymptotic products take 59 ns for JVP and 57 ns for VJP,
3.89× and 7.67× faster than Enzyme. The selected generated kernels also use
slightly less process RSS and fewer measured instructions.

The completed special-function records report the observed finite-difference
error, timing MAD, candidate-attributed peak RSS, native code size, reusable
state, product-count scaling, and separate JVP/VJP results. Maximum observed
relative errors are \(2.90\times10^{-10}\) for Bessel,
\(4.03\times10^{-12}\) for gamma, \(2.38\times10^{-11}\) for erf, and
\(4.19\times10^{-10}\) for hypergeometric. Median candidate RSS at 16 products
ranges from 2.72 to 3.37 MB across these fixtures.

The completed interpolation tournament covers active Lagrange nodes, active
B-spline knots, implicit spline fitting, and fixed-span analytical versus
autodiff products. At the largest measured sizes, analytical JVP/VJP speedups
over the finite-difference diagnostics are 2.31×/16.61× for Lagrange nodes,
2.04×/2.09× for active knots, and 4.90×/1,224.35× for spline fitting.
Generated-wrapper fixed-span analytical products take 47.46 ns per 16 JVPs and
52.36 ns per 16 VJPs, versus 59.40 and 58.31 ns for Enzyme. Independent
finite-difference errors at the largest sizes range from
\(8.71\times10^{-12}\) to \(7.20\times10^{-9}\). Each record identifies
reused spline or factorization state, native code size, candidate memory,
direction scaling, and cache counters.

For the length-eight complex FFT JVP, analytical and Enzyme candidates share
the same production radix execution leaf. Analytical takes 251, 655, and
2,464 ns for 1, 4, and 16 directions; Enzyme takes 274, 696, and 2,360 ns.
The 16-direction raw gap is 4.4%, but repeated hardware counters favor
analytical, which also has lower RSS, fewer instructions and misses, and is
available in normal and device builds. Plan construction is correctly
inactive: Enzyme 22 cannot differentiate Flang's `_FortranAAssign`, and plan
derivatives are not part of the mathematical FFT contract.

For the corresponding VJP, the analytical adjoint transform takes 272, 691,
and 2,420 ns for 1, 4, and 16 cotangents. A valid autodiff construction uses
16 Enzyme forward-mode basis products per cotangent and takes 859, 2,920, and
11,865 ns. The analytical VJP is 3.16 to 4.90 times faster and uses 2.81 MB
peak RSS versus 2.90 MB. Direct reverse Enzyme through the in-place complex
radix kernel produced NaNs and failed the direct-adjoint oracle, so it is a
rejected candidate rather than a timing result.

The external-library pilot wraps an FFTW 3.3.11 plan in a scalar objective.
Enzyme differentiates the surrounding Fortran expression and invokes an
analytical forward rule at the FFTW boundary. A generated forward-only
fixed-array wrapper prevents an unsupported reverse request from entering the
Enzyme module. Analytical and hybrid JVPs differ by at most 3.3% over 1, 4,
and 16 directions. At 16 directions they take 4,745 and 4,897 ns, with peak
RSS of 4.94 and 4.99 MB. The explicit analytical composition remains selected.

A complete frequency-weighted spectral objective maps 16 real components to
one scalar and differentiates the composition, not an isolated transform. One
primal FFT followed by the analytical adjoint takes 196 ns. Contracting the
same gradient from 16 Enzyme forward-mode basis products takes 776 ns, so the
analytical composition is 3.96 times faster. Peak RSS is 2.92 versus 3.00 MB.
The analytical path also executes 4.13 times fewer instructions and incurs
5.43 times fewer measured cache misses. An independent direct DFT and adjoint,
plus a central finite difference of the complete objective, validate both
candidates.

The completed FFT tournament records observed oracle errors, timing MAD,
candidate RSS, native code size, reusable plan state, product-count scaling,
and CPU counters. The shared native FFT fixture's maximum scaled error is
\(5.63\times10^{-11}\), including the complete-objective finite difference.
The FFTW hybrid fixture's maximum relative error is
\(6.13\times10^{-11}\). Production selections remain the analytical radix JVP
and adjoint VJP, explicit analytical composition around FFTW, and the
analytical spectral-objective gradient.

The completed quadrature tournament covers fixed JVP/VJP products, smooth and
singular frozen adaptive traces, and batched full Jacobians. For four fixed
products, analytical JVP takes 1,223 ns versus 1,225 ns hybrid, 1,702 ns
autodiff, and 2,212 ns finite-difference diagnostic. One analytical VJP returns
the four-parameter gradient in 301 ns, 4.07× faster than four analytical JVPs.
At batch size 16, the complete analytical reverse Jacobian takes 4.75
microseconds, versus 19.54 microseconds for analytical forward products and
6.98 microseconds for reverse autodiff. Smooth and singular adaptive winners
remain within 1.1% of their nearest validated competitors. Maximum observed
fixed-rule errors are \(4.11\times10^{-11}\) for JVP and
\(7.00\times10^{-11}\) for VJP. The singular frozen-trace candidates agree
within \(4.94\times10^{-11}\); their 0.0116 difference from the continuous
integral is the documented frozen-trace bias.

Single-transform analytical products were measured at lengths 8, 64, 256,
and 1024. JVP wall clock grows from 178 ns to 27.0 microseconds. VJP grows
from 170 ns to 26.3 microseconds and stays within 6.3% of JVP at every length.
Length-1024 peak RSS is 2.86 MB for JVP and 2.91 MB for VJP. These complete
public calls include plan construction for lengths other than 8. The scalar
record therefore describes the current caller-visible cost, including its
lack of reusable complex-plan state.

The hypergeometric tournament activates real `z` at fixed real parameters.
Because Enzyme 22 cannot differentiate Flang's complex `cexp`, its candidate
uses a real specialization of the same Kummer, Taylor, and asymptotic
algorithms. Primal equivalence to the production complex path is a mandatory
validation gate. In the `x > 60` region, fortsym also emits pure-elemental
analytical JVP and VJP kernels from the asymptotic definition. At 16 products
they take 59.10 and 56.55 ns, versus Enzyme's 229.71 and 433.71 ns. Peak RSS
differs by less than 2%, so wall clock selects the generated region kernels.

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
