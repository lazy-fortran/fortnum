# Derivatives in fortnum

Many fortnum procedures come with derivative products: forward products (JVPs),
reverse products (VJPs), gradients and, where useful, Hessian-vector products.
For example `fortnum_special` provides `dawson_jvp` next to `dawson`.

## Terminology

- `autodiff`: compiler or source-transformation differentiation.
- `analytical`: explicit formulas, stable recurrences, linear operators,
  differentiation of a recorded computation trace, and implicit tangent or
  adjoint solves.
- `hybrid`: autodiff composed with at least one analytical rule.

A procedure can have several implementations of the same product. Each is
checked against independent reference values, and the one used by default is
chosen from measured wall-clock time and peak memory on representative workloads.
Finite differences are used for validation, not as production derivatives.
Quantities defined by a residual equation always have an implicit
(adjoint-solve) implementation.

## CPU and GPU support

| Target | `autodiff` | `analytical` | `hybrid` |
| --- | --- | --- | --- |
| CPU | selected kernels | supported | selected kernels |
| GPU | unsupported | listed kernels | unsupported |

## Generated code

Symbolic algebra and code generation are done at build time by
[fortsym](https://github.com/lazy-fortran/fortsym). Generated production
kernels are committed under `src/generated`; temporary variants and Enzyme
wrappers stay in the build tree.

## Benchmarks

Measurements are committed under `benchmark/reference`, with cross-kernel
statistics in `benchmark/report/data/mechanism_tournaments.csv`. Each record
includes the validation result, workload shape, compiler and hardware, wall-clock
time and peak memory. See [benchmark/README.md](../benchmark/README.md) for how
to reproduce the reports.

## Further reading

- [Derivative contract](design/ad.md)
- [Differentiation plan](design/differentiation_plan.md)
- [Kernel ownership inventory](design/derivative_kernel_inventory.csv)
- [CPU/GPU contract](design/gpu.md)
- [Evidence report](design/differentiation_report.md)
- [Choosing the differentiation backend](api.md#choosing-the-differentiation-backend)
