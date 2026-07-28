# Differentiation evidence report

Status: normalized report for the committed tournament CSV.

Post-generation-migration refresh: 2026-07-28. Revalidation retained the
existing selections and the gamma tournament added six product-specific
results. The eight temporary figures generated in 0.72 s wall clock with
29.4 MB peak RSS; generated PNGs were inspected but are not repository
artifacts.

## Scope

`benchmark/report/data/mechanism_tournaments.csv` contains 38 derivative
workloads with competing mechanism timings. Every row names its source JSON
record.

`selected_mechanism` is the admissible production decision.
`fastest_mechanism` is the lowest raw median. They differ when the raw winner
requires an optional toolchain that is unavailable in the normal build.

Finite differences are diagnostics and win no current production selection.

## Aggregate result

| Mechanism | Selected workloads | Share | Raw fastest workloads |
| --- | ---: | ---: | ---: |
| `analytical` | 29 | 76.3% | 26 |
| `autodiff` | 7 | 18.4% | 7 |
| `hybrid` | 2 | 5.3% | 5 |
| diagnostic | 0 | 0.0% | 0 |

The second-fastest/fastest wall-clock ratio ranges from 1.001 to 1,224.355.
Its median is 1.465 and geometric mean is 2.005. 13 workloads have a
runner-up within 20%.

These counts describe the measured workload set. They do not select new
kernels.

The migration replaced eligible transcribed algebra and Enzyme scaffolding
with generated sources. Its same-mechanism raw/simplified/factored and
generated/manual records remain outside this table, so they cannot inflate the
cross-mechanism totals. No measured production selection changed in this
refresh.

## Current verdict by family

| Family | Current measured verdict | Reason |
| --- | --- | --- |
| compact generated algebra | `analytical` | direct contracted and fused products avoid repeated primals |
| Bessel | mixed | raw forward autodiff reuses recurrence work in two JVP regions; analytical wins the measured VJPs and series JVP |
| regularized gamma, fixed shape | product-specific | forward autodiff reuses the primal iteration and wins JVPs; the analytical scalar adjoint wins VJPs |
| fixed quadrature | `analytical` | fixed linear contractions beat Enzyme and complete finite differences |
| adaptive integration | mixed | whole-trace autodiff narrowly wins two measured JVP workloads; normal-build analytical remains for an optional hybrid raw win |
| scalar/vector roots | product-specific | hybrid residual JVP can win; analytical implicit VJP and factor reuse remain strongest |
| interpolation and spline fits | `analytical` | basis products and implicit coefficient solves exploit structure |
| direct and iterative solves | `analytical` | factor-reusing tangent/adjoint equations beat iteration autodiff |
| ODE | `analytical` for current selections | frozen recurrences and adjoints reuse the accepted schedule; optional hybrid RHS JVP is a narrow raw win |

Forward and reverse mode require separate tournaments. In the fixed four-by-four
direct solve, 16 analytical JVPs take 0.643 microseconds and 16 analytical VJPs
take 0.637 microseconds. The corresponding forward and reverse Enzyme times are
2.449 and 0.974 microseconds. This one workload shows reverse Enzyme approaching
the analytical adjoint more closely. It does not establish a general crossover.

## Figures

The `fortplot` generator creates:

- selected-mechanism counts
- a histogram of runner-up slowdown
- direct-solver JVP and VJP scaling
- fixed-span B-spline JVP and VJP scaling
- Bessel-region JVP and VJP comparisons

Generate them:

```bash
output_dir=$(mktemp -d)
cd benchmark/report
fo exec plot_differentiation_report \
  data/mechanism_tournaments.csv "${output_dir}"
```

The repository commits the generator and normalized CSV. PNG outputs remain
outside version control.

## Interpretation

Wall clock decides among validated candidates for the exact complete workload.
Peak memory can reject a time winner or break a practical timing tie. Cache and
work counters explain results but do not replace wall clock.

Re-run selection after material changes to primal code, compiler, Enzyme,
`fortsym`, workload shape, or target hardware.
