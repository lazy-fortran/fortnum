## Scope

Closes #

What behavior changes, and what remains outside this pull request?

## Numerical contract

| Procedure or kernel | Active arguments | Product | Candidates | Selected candidate |
| --- | --- | --- | --- | --- |
| | | | | |

For residual-defined outputs, include the analytical implicit candidate.

## Validation

What independent oracle checks the behavior? Include boundary, singular, or
trace cases relevant to the change.

## Performance

Give complete-workload median wall clock, dispersion, and candidate-specific
peak memory. State the comparison baseline, workload dimensions, compiler,
hardware, and affinity. Add scaling and cache or device counters when they
affect selection.

## Generated sources and documentation

- [ ] Generated files reproduce byte-for-byte, or this change has none
- [ ] README, API, design, migration, and benchmark docs are current
- [ ] No generated PNG files are committed

## Verification

Paste the commands and result counts:

```text
fo
ctest --test-dir build --output-on-failure
```
