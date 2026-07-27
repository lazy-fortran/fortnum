## Issue

Closes #

## Summary

What does this add or change, and why?

## Derivative candidates

For each new public procedure and product, list the admissible `autodiff`,
`analytical`, and `hybrid` candidates. Identify `primal_only` operations with a
justification. See `docs/design/ad.md`.

| Procedure | Product | Candidates | Active args | Validation oracle |
|---|---|---|---|---|
| | | | | |

If this PR selects a production winner, link runtime and peak-memory evidence.

## Verification

Paste real `ctest` output below. Do not paste partial output or summarize it.

```
ctest --test-dir build --output-on-failure
```

All tests must pass. Do not weaken or skip failing tests.
