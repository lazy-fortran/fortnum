# Differentiation evidence report

This report summarizes the committed candidate tournaments without replacing
their workload-specific decisions. Complete workload wall clock is the primary
metric. Peak memory, validation error, scaling, and cache counters remain in
the source records under `benchmark/reference/`.

## Scope

The cumulative data set contains 32 measured derivative-product tournaments
for which competing mechanism timings are committed. It excludes comparisons
between two implementations of the same mechanism, incomplete selections, and
reliability-only experiments. Finite differences remain a validation
`diagnostic`, not a production mechanism.

`benchmark/report/data/mechanism_tournaments.csv` is the normalized plotting
input. Every row names its source JSON record. `selected_mechanism` records the
evidence-based production decision; `fastest_mechanism` records the lowest raw
median. They differ for Dawson, smooth adaptive integration, and ODE forward
sensitivity, where raw `hybrid` advantages are available only through the
optional Enzyme pipeline while the selected analytical implementations are
the fastest candidates admissible in a normal fortnum build.

## Overall result

| Mechanism | Selected workloads | Share | Raw fastest workloads |
|---|---:|---:|---:|
| `analytical` | 26 | 81.2% | 23 |
| `autodiff` | 4 | 12.5% | 4 |
| `hybrid` | 2 | 6.2% | 5 |
| finite-difference diagnostic | 0 | 0.0% | 0 |

Across the 32 workloads, the second-fastest candidate ranges from 1.001x to
1,224.360x the fastest wall clock. The median ratio is 1.275x and the
geometric mean is 2.076x. Thirteen workloads have a runner-up within 20% of the
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
| 1 | 0.041 us | 0.152 us | 0.049 us | 0.046 us | 0.081 us | 0.868 us |
| 4 | 0.163 us | 0.612 us | 0.200 us | 0.164 us | 0.241 us | 3.523 us |
| 16 | 0.643 us | 2.449 us | 0.816 us | 0.637 us | 0.974 us | 14.037 us |

Both forward and reverse products scale approximately linearly with the number
of requested directions or cotangents in this small fixed-size workload.
Analytical is fastest in both cases. Reverse autodiff is much closer to
analytical than forward autodiff: at 16 products analytical is 1.529x faster
than reverse autodiff, while it is 3.812x faster than forward autodiff. This
does not establish a universal forward/reverse crossover; it shows why the two
products need separate tournaments as input and output dimensions change.

The fixed-span cubic B-spline comparison provides a smaller interpolation
kernel. After migration to generated wrappers and pre-loop candidate dispatch,
analytical is 20.1% faster than forward Enzyme for 16 JVPs and 10.2% faster
than reverse Enzyme for 16 VJPs. The selected complete-workload wall clocks
also improve by 66.4% and 60.4%, respectively, relative to the old fixture.

The two-parameter scalar root shows that forward and reverse decisions can
still differ at the same interface. Generated Enzyme residual products make
the hybrid implicit JVP 1.215x faster than analytical, while the analytical
implicit VJP is 1.381x faster than the reverse hybrid. Autodiff through the
fixed Newton iterations loses both tournaments.

The coupled two-state root makes the same point with a matrix implicit solve.
The generated-wrapper migration leaves the hybrid JVP narrowly fastest at
78.27 ns versus 78.31 ns analytical, but changes the VJP selection to
analytical at 78.80 ns versus 79.82 ns hybrid. Autodiff through 12 Newton
iterations is 2.79x slower for JVP and 2.62x slower for VJP than the selected
implicit candidates. Cache misses are nearly identical, so the extra executed
work—not cache behavior—explains the larger wall clock.

Fixed quadrature now uses generated Enzyme wrappers for both the five-input
integrand and four-input complete operator. For its four-input, one-output
shape, four analytical forward JVPs take 1.223 µs while one analytical reverse
VJP takes 0.301 µs. The corresponding autodiff times are 1.702 µs and
0.440 µs. Analytical remains selected for both products; peak RSS and cache
misses do not materially distinguish the leading candidates.

Adaptive integration generates only the smooth/singular integrand and compact
frozen-trace Enzyme boundaries. The accepted-trace construction and identity
guards remain explicit. On the smooth workload, compact hybrid is the
Enzyme-enabled raw winner at 1.260 µs, followed by whole-trace autodiff at
1.264 µs and compact analytical at 1.269 µs; analytical remains the normal
build selection. The singular workload selects whole-trace autodiff at
2.807 µs over compact analytical at 2.836 µs.

ODE forward sensitivity similarly generates only the local three-input RHS
JVP boundary. The Enzyme-enabled hybrid takes 9.071 µs for the complete primal
trajectory plus sensitivity, versus 9.169 µs analytical and 14.846 µs for
complete-solve finite differences. Analytical remains selected for normal
builds; cache events are too dispersed to override wall clock.

The modified-Bessel outer objective strengthens that conclusion across primal
regions. At 16 products, analytical wins series JVP/VJP, recurrence VJP, and
asymptotic VJP; raw Enzyme wins recurrence and asymptotic JVP. The forward
hybrid custom rule is validated but does not win because its `I0` and `I1`
evaluations duplicate work that raw forward Enzyme reuses. The recurrence VJP
selection changed from autodiff to analytical after the shared-fixture
migration; the controlled winner is 0.7% faster than the pre-migration winner.

## Figures

The generator creates eight independent PNGs:

- selected-mechanism workload counts
- a histogram of
  `log10(second-fastest wall clock / fastest wall clock)`
- complete-wall-clock scaling with forward JVP directions
- complete-wall-clock scaling with reverse VJP cotangents
- fixed-span B-spline JVP scaling
- fixed-span B-spline VJP scaling
- Bessel JVP mechanisms across three primal regions
- Bessel VJP mechanisms across three primal regions

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
