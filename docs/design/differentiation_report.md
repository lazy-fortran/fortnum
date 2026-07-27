# Differentiation evidence report

This report summarizes the committed candidate tournaments without replacing
their workload-specific decisions. Complete workload wall clock is the primary
metric. Peak memory, validation error, scaling, and cache counters remain in
the source records under `benchmark/reference/`.

## Scope

The cumulative data set contains 26 measured derivative-product tournaments
for which competing mechanism timings are committed. It excludes comparisons
between two implementations of the same mechanism, incomplete selections, and
reliability-only experiments. Finite differences remain a validation
`diagnostic`, not a production mechanism.

`benchmark/report/data/mechanism_tournaments.csv` is the normalized plotting
input. Every row names its source JSON record. `selected_mechanism` records the
evidence-based production decision; `fastest_mechanism` records the lowest raw
median. They differ for Dawson, where the 2% raw `hybrid` advantage is within
measurement noise and the smaller analytical implementation remains selected.

## Overall result

| Mechanism | Selected workloads | Share | Raw fastest workloads |
|---|---:|---:|---:|
| `analytical` | 19 | 73.1% | 18 |
| `autodiff` | 3 | 11.5% | 3 |
| `hybrid` | 4 | 15.4% | 5 |
| finite-difference diagnostic | 0 | 0.0% | 0 |

Across the 26 workloads, the second-fastest candidate ranges from 1.004x to
1,224.355x the fastest wall clock. The median ratio is 1.179x and the
geometric mean is 2.144x. Sixteen workloads are tightly grouped near the
fastest candidate; full finite-difference VJPs of fitted coefficients provide
the largest separation.
Consequently, mechanism counts describe only the current measured workload
set and are not a rule for selecting a new kernel.

## Forward and reverse direct-solver evidence

The fixed 4-by-4 direct solver provides a controlled comparison between
forward JVP and reverse VJP products. Times below are complete product wall
clock on the AMD Ryzen 9 5950X reference host.

| Product count | Analytical JVP | Autodiff JVP | Diagnostic JVP | Analytical VJP | Autodiff VJP | Diagnostic VJP |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.055 us | 0.189 us | 0.072 us | 0.059 us | 0.074 us | 1.025 us |
| 4 | 0.208 us | 0.721 us | 0.257 us | 0.227 us | 0.308 us | 3.972 us |
| 16 | 0.774 us | 2.719 us | 0.982 us | 0.853 us | 1.090 us | 15.455 us |

Both forward and reverse products scale approximately linearly with the number
of requested directions or cotangents in this small fixed-size workload.
Analytical is fastest in both cases. Reverse autodiff is much closer to
analytical than forward autodiff: at 16 products analytical is 1.278x faster
than reverse autodiff, while it is 3.514x faster than forward autodiff. This
does not establish a universal forward/reverse crossover; it shows why the two
products need separate tournaments as input and output dimensions change.

The fixed-span cubic B-spline comparison provides a smaller interpolation
kernel where neither mechanism dominates. At 16 products, forward Enzyme JVP
is 2.8% faster than analytical, while analytical VJP is 6.6% faster than
reverse Enzyme. The selected mechanism therefore changes with both derivative
product and workload shape.

## Figures

The generator creates six independent PNGs:

- selected-mechanism workload counts
- a histogram of
  `log10(second-fastest wall clock / fastest wall clock)`
- complete-wall-clock scaling with forward JVP directions
- complete-wall-clock scaling with reverse VJP cotangents
- fixed-span B-spline JVP scaling
- fixed-span B-spline VJP scaling

The scaling figures use microseconds and encode mechanisms with an
Okabe-Ito-derived color-safe palette plus distinct line styles and markers.
The checked-in repository contains the generator and data, never generated
PNGs.

Reproduce the figures from the repository root:

```bash
mkdir -p /tmp/fortnum-differentiation-report
cd benchmark/report
fo exec plot_differentiation_report \
  data/mechanism_tournaments.csv /tmp/fortnum-differentiation-report
```

The report package pins the exact `fortplot` revision. The generated color and
grayscale figures were inspected for readable labels, units, legend placement,
and redundant series distinction.
