# Documentation map

The production source is authoritative for exact interfaces and arithmetic.
The documents below define contracts, explain use, track work, or summarize
measured evidence.

## User guides

| Document | Scope |
| --- | --- |
| [README](../README.md) | Build, package scope, derivative model, and entry points |
| [API guide](api.md) | Public module families and common call patterns |
| [libneo migration](migration_libneo.md) | Primal API mapping |
| [libneo derivative migration](migration_libneo_ad.md) | Derivative-product mapping |

## Design contracts

| Document | Scope |
| --- | --- |
| [Architecture](design.md) | Module boundaries and repository invariants |
| [Derivative contract](design/ad.md) | Terminology, products, candidate selection, and generated-source ownership |
| [Optimizer API](design/optimizer_api.md) | Flat active vectors and backend-independent callbacks |
| [Downstream active kernels](design/downstream_ad.md) | Packing and composing downstream differentiable kernels |
| [Enzyme toolchain](design/enzyme_toolchain.md) | Supported CPU Flang/Enzyme boundary |
| [CPU/GPU contract](design/gpu.md) | Device leaves, supported mechanisms, and offload gates |
| [Adaptive integration](design/integrate.md) | Integration state, status, traces, and derivative semantics |
| [ODE](design/ode.md) | Solver state, events, traces, and sensitivity semantics |
| [RNG](design/rng.md) | Explicit state, stream splitting, and reproducibility |

## Plans and evidence

| Document | Scope |
| --- | --- |
| [ROADMAP](../ROADMAP.md) | Authoritative checklist and execution rules |
| [Differentiation plan](design/differentiation_plan.md) | Current differentiation architecture and implementation sequence |
| [Kernel inventory](design/derivative_kernel_inventory.csv) | Machine-checked derivative source ownership |
| [Evidence catalog](design/differentiation_benchmarks.md) | Benchmark records and validation coverage |
| [Evidence report](design/differentiation_report.md) | Aggregate mechanism statistics and generated figures |
| [Strategy paper](performance_optimal_differentiation.md) | Rationale for benchmark-selected differentiation |

Machine-readable benchmark records under `benchmark/reference/` own exact
measurements. Documentation cites those records and avoids duplicating
append-only run logs.

## Maintenance

Update the relevant user guide and design contract with every behavior change.
Update benchmark documents only when validation or selection evidence changes.
Do not hard-code dependency revisions in prose. Link the lock file or benchmark
record that owns the revision.

CTest runs `scripts/check_documentation.py`. The check is quiet on success and
rejects broken local links, missing document-map entries, undocumented
production modules, absent API symbols, legacy terminology, stale issue-era
claims, hard-coded revisions, generator-lock drift, stale aggregate benchmark
statistics, unindexed derivative ownership, and committed report PNGs.

CTest also rejects generated JVP or VJP paths that call a Jacobian constructor,
use the wrong product construction, or lose their generator ownership.
