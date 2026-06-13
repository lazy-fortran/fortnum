---
name: Module port
about: Port an existing algorithm from a published reference
labels: port
---

## Routine

Name and domain (e.g. `bessel_j0`, special functions).

## Reference implementation

Link or citation: library, file, function name, version or commit.

Note any license restrictions. fortnum is MIT; the port must be a
clean-room reimplementation if the source is GPL or LGPL.

## Derivative policy

Select one per procedure or mark `TBD`. See `docs/design/ad.md` for
definitions and annotation syntax.

- [ ] `transparent`
- [ ] `analytic_rule`
- [ ] `implicit_rule`
- [ ] `trace_rule`
- [ ] `primal_only` (justify below)

Justification for `primal_only` (if applicable):

## Active arguments

List which dummy arguments carry derivative information and which are inactive
(indices, tolerances, mode flags).

## Oracle reference

Where will the oracle test get its reference values? (scipy, DLMF table,
analytic result, …)

## Milestone

Which milestone does this belong to? (M1 kernels, M2 ODE, M3 roots/interp/rng,
M4 adaptive integration)
