---
name: Clean-room implementation
about: Implement a routine from primary sources with no reference code
labels: new-routine
---

## Routine

Name and domain (e.g. `dop853`, ODE solvers).

## Primary sources

List the algorithms, papers, or standards (DLMF, textbooks) the implementation
will follow. No reference to GPL/LGPL source code.

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
analytic result, manufactured solution, …)

## Public interface sketch

Paste a draft Fortran interface or subroutine signature. Enough to agree on
argument names, intents, and types before implementation starts.

## Milestone

Which milestone does this belong to? (M1 kernels, M2 ODE, M3 roots/interp/rng,
M4 adaptive integration)
