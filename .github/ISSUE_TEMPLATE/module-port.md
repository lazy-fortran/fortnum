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

## Derivative candidates

Select every admissible candidate per derivative product or mark `TBD`. See
`docs/design/ad.md`.

- [ ] `autodiff`
- [ ] `analytical`
- [ ] `hybrid`
- [ ] `finite_difference_reference`
- [ ] `primal_only` (justify below)

Justification for `primal_only` (if applicable):

Requested products and candidate details:

Validation oracle:

Representative benchmark workloads:

## Active arguments

List which dummy arguments carry derivative information and which are inactive
(indices, tolerances, mode flags).

## Oracle reference

Where will the oracle test get its reference values? (scipy, DLMF table,
analytic result, …)

## Milestone

Which milestone does this belong to? (M1 kernels, M2 ODE, M3 roots/interp/rng,
M4 adaptive integration)
