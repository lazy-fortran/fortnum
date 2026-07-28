# Architecture

Status: current repository contract.

## Package boundaries

`fortnum` is a static Fortran 2018 numerical library. Numerical areas live in
separate modules:

```text
special       special functions
quadrature    fixed rules, adaptive integration, sequence acceleration
fft           real and complex transforms
ode           ODE methods, traces, events, and sensitivities
roots         scalar, vector, reverse-communication, and complex roots
interp        polynomial interpolation and B-splines
linalg        fixed-size primitives and small LU solves
rng           explicit-state random generation
ad            derivative interfaces, active-vector layout, and selection
bindings      C ABI adapters
generated     committed fortsym numerical leaves
testing       oracle-table support
```

The exact source list is in `src/CMakeLists.txt`. LAPACK and BLAS are the only
runtime numerical dependencies.

## State ownership

The library has no mutable global numerical state. Work arrays, solver traces,
plans, factorizations, RNG counters, and callback context belong to the caller
or an explicit derived type.

Callbacks accept optional caller context where needed. Internal host-associated
wrappers may adapt a context-bearing callback to a simpler local interface.
They exist only for the duration of the call.

This ownership model supports:

- independent concurrent calls with separate state
- deterministic tests
- explicit memory accounting
- derivative replay of accepted traces
- reuse of primal factors and plans

## Kinds and status

Public numerical code uses `real64` through `dp`. The status module defines
success, domain, convergence, and unimplemented codes. A routine that can fail
returns `fortnum_status_t` or its documented compatibility status.

Failure paths initialize outputs before returning. Numerical routines do not
print, stop the process, or update global error state.

## Numerical semantics

IEEE behavior is part of adaptive integration, root, and special-function
contracts. The CMake target disables inherited fast-math modes that would
remove NaN checks or reassociate validated sums.

Primal code is `pure` when its algorithm and callback interface allow it.
Hot loops avoid allocation. Fixed-capacity local algorithms declare their
limits in the public interface.

## Differentiation

The stable public products are value, JVP, VJP, gradient, and HVP interfaces.
Each product may have `autodiff`, `analytical`, and `hybrid` candidates.
Selection uses independent validation plus complete-workload wall clock and
peak memory.

Mathematical operator boundaries determine custom rules:

- linear transforms apply the operator or its adjoint
- roots and solves use implicit tangent or adjoint equations
- adaptive methods replay a documented frozen trace
- stable special-function recurrences remain intact
- explicit local algebra can be generated

See [design/ad.md](design/ad.md) for the normative derivative contract.

## Generated sources

`fortsym` runs at development and build time. Symbolic definitions generate
local explicit products and Enzyme wrappers. Selected production kernels are
committed under `src/generated` with provenance and regeneration commands.
Temporary candidates and routine Enzyme wrappers stay in the build tree.

Stable recurrences, region selection, solver orchestration, traces, scheduling,
and independent oracles remain hand-written.

## CPU and GPU

CPU builds support selected `autodiff`, `analytical`, and `hybrid` candidates.
GPU builds support individually validated `analytical` device leaves.
OpenACC and OpenMP target wrappers use the same pure numerical body. Backend
selection occurs outside launch loops, and host fallback is an error.

See [design/gpu.md](design/gpu.md) for the support matrix and acceptance gates.

## Verification

Unit tests check local behavior. Oracle tests use exact results, identities,
manufactured problems, or independent reference tables. Derivative tests add
directional checks, adjoint identities, and residual linearizations.

Benchmark records under `benchmark/reference` own exact performance
measurements. Figure generators and normalized report data live under
`benchmark/report`.
