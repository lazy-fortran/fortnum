---
name: Module port
about: Reimplement a published numerical algorithm
labels: port
---

## Operator

Name, mathematical definition, domain, and expected public module.

## Sources and license

List the paper, standard, textbook, or reference implementation with its exact
version. State the source license and the clean-room boundary required for an
MIT implementation.

## Interface

Sketch the Fortran signature. Identify active inputs, inactive controls,
outputs, status behavior, and numerical regimes.

## Derivative products

| Product | `autodiff` candidates | `analytical` candidates | `hybrid` candidates |
| --- | --- | --- | --- |
| JVP | | | |
| VJP | | | |
| gradient or HVP | | | |

Explain any unsupported product. Include an analytical implicit candidate when
the output is defined by a residual.

## Validation

Name an independent oracle and its boundary or singular test cases.

## Performance

List representative dimensions, batches, active-input/output counts,
directions, reusable primal state, and target hardware.
