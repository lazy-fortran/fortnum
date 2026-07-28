---
name: Clean-room implementation
about: Implement a numerical operator from primary sources
labels: new-routine
---

## Operator

Name, mathematical definition, domain, and expected public module.

## Primary sources

List equations, algorithms, standards, and test data. Do not consult
license-incompatible implementation source.

## Interface

Sketch the Fortran signature. Identify active inputs, inactive controls,
outputs, status behavior, numerical regimes, and stability requirements.

## Derivative products

| Product | `autodiff` candidates | `analytical` candidates | `hybrid` candidates |
| --- | --- | --- | --- |
| JVP | | | |
| VJP | | | |
| gradient or HVP | | | |

Explain any unsupported product. Include an analytical implicit candidate when
the output is defined by a residual.

## Validation

Name an independent oracle and include boundary, regime-transition, singular,
and failure-status cases.

## Performance

List representative dimensions, batches, active-input/output counts,
directions, reusable primal state, and target hardware.
