# Array-temporary baseline

Status: 74 compiler warnings explicitly baselined on 2026-07-28.

`fo lint` reports array temporaries in 27 files. Of these, 24 warnings are in
production `src/`; the remaining 50 are in benchmark fixtures, Enzyme probes,
examples, and tests. A warning is not evidence of a performance regression by
itself. Complete-workload wall clock decides whether replacing a temporary is
worthwhile.

The stable baseline records counts by file rather than line number:

```bash
tools/check-array-temporaries.sh
```

The checker runs the real compiler lint and rejects new, removed, or relocated
warnings. When a warning changes, measure the affected workload before
updating `tools/array_temporary_baseline.csv`. Remove a baseline entry as soon
as its temporary is eliminated.
