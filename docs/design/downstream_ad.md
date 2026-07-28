# Downstream active kernels

Status: integration guide.

Downstream applications should expose a small differentiable kernel around one
mathematical objective or residual. Solver setup, file I/O, allocation,
logging, and optimizer control stay outside that kernel.

## Boundary

A useful active kernel has:

- explicit real-valued active inputs
- explicit outputs
- immutable configuration in a context object
- no mutable global state
- stable primal calls into `fortnum`
- value, JVP, VJP, gradient, or HVP callbacks required by the consumer

Discrete mode flags, dimensions, tolerances, mesh topology, and random seeds
remain inactive.

## Packing

Build one `fortnum_active_layout_t` for the optimizer variables. Write and read
named blocks at the boundary:

```text
optimizer x(:)
    -> read named geometry, profile, and control blocks
    -> evaluate objective or residual
    -> pack gradients or cotangents into the same layout
```

The layout order is part of the downstream API. Store it with checkpoint or
restart metadata when vectors persist across runs.

## Composition

Call the derivative product matching the outer algorithm:

- forward JVPs for a small number of parameter directions
- reverse VJPs for a small number of output cotangents
- analytical implicit products for equilibria and solves
- frozen-trace products for a specified discrete adaptive map

Do not build a full Jacobian only to contract it once.

A typical equilibrium objective uses:

1. a primal equilibrium solve
2. local residual JVP/VJP products
3. an analytical implicit tangent or adjoint solve
4. an objective JVP/VJP
5. packing into the optimizer vector

If local products use autodiff and the equilibrium solve uses an analytical
implicit rule, the complete candidate is `hybrid`.

## Context and concurrency

The callback context owns grids, plans, traces, factors, and fixed physics
data. Give each concurrent evaluation separate mutable workspace. Shared
read-only tables are safe when their lifetime exceeds every callback.

## Validation

Validate the complete downstream kernel, not only its local library calls:

- objective gradients against directional differences
- residual JVP/VJP dot-product identities
- implicit residual linearizations
- invariance under `pack_block` and `unpack_block` round trips
- frozen-trace identity checks
- representative application wall clock and peak memory

Microbenchmark winners can lose after batching, cache pressure, factor reuse,
or optimizer iteration count. Select with the complete downstream workload.
