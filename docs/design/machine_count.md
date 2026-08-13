# Machine-code operation counts (N_machine)

Status: implemented.

`N_machine` is the third level of the four-count cost chain:

| Level | Meaning | Source |
| --- | --- | --- |
| `N_sym` | symbolic DAG as first built | `fortsym` |
| `N_emit` | what the generator wrote (post-CSE roots) | `fortsym` `count_operations` |
| `N_machine` | what the compiler produced | this collector, from disassembly |
| measured time | end-to-end wall clock and peak memory | benchmark harness |

`N_emit` is what the generator wrote; `N_machine` is what the compiler
produced. The gap is the compiler's doing and moves in both directions:

- **down**, when FMA fusion turns two floating-point operations into one
  instruction;
- **up**, when register pressure forces spills, adding loads and stores that
  never appeared in the source.

That second case is the observable consequence of over-aggressive CSE
(lazy-fortran/fortsym#64), and it is invisible at every level above the
machine code.

## Tool

`benchmark/collect_machine_count.py` counts floating-point instructions in
the compiled artifact, per target:

- **CPU**: `objdump -d` on the object file
- **CUDA**: `cuobjdump -sass` / `nvdisasm` for SASS. SASS, not PTX: PTX is not
  the machine code and ptxas changes it substantially.

```bash
python3 benchmark/collect_machine_count.py \
  --object /path/to/kernel.o \
  --source src/generated/fortnum_det2_jvp_kernel.f90 \
  --target cpu \
  --output benchmark/reference/xeon_e5_2630v4_det2_machine_count.json
```

`--source` supplies `N_emit` from the emitted Fortran source so the side-by-side
comparison is one command and does not depend on a separate `fortsym` run.
`N_emit` and `N_machine` are kept in separate manifest blocks and never
conflated.

## Reported numbers

Reported separately, never conflated:

- **instruction count** and **FLOP count**. One FMA is one instruction and two
  FLOPs. Conflating them is the standard error in this measurement and makes
  every derived number wrong.
- **spill loads and stores**, counted apart from ordinary memory traffic. On
  x86-64 a spill is a floating-point load or store through the stack frame
  (`%rsp`/`%rbp`); in SASS it is an `LDL`/`STL` local-memory access.
- **register count**, where the tool reports it. `nvdisasm` reports it;
  `objdump` and `cuobjdump -sass` do not, so the field is absent on those
  targets.

## Opt-in

The binary utilities are not present everywhere. When a tool is absent the
level is simply not measured: the manifest field is absent, never `"skipped"`,
never zero. Absence is the honest encoding of "no evidence" and keeps the
record from asserting more than it knows. CI stays portable; local runs and
the benchmark machine get full attribution.

## Acceptance demonstration

For the generated determinant-2 JVP kernel
(`src/generated/fortnum_det2_jvp_kernel.f90`), compiled with
`gfortran -O3 -march=x86-64-v3`:

```
jvp = a*vd - b*vc - c*vb + d*va
```

- `N_emit`: 4 multiplies + 3 add/subtract = 7 instructions, 7 FLOPs
- `N_machine`: 1 multiply + 3 FMA = 4 instructions, 7 FLOPs

The FMA fusion is visible as a reduction (7 instructions to 4) and the FLOP
count is preserved (7 to 7). The committed record
`benchmark/reference/xeon_e5_2630v4_det2_machine_count.json` captures this.

For kernels with division or transcendental calls the compiler additionally
reassociates and folds constants, so the machine FLOP count may diverge from
the emitted source count even without spills; the collector reports the
measured `flop_count_preserved` verdict rather than asserting one. Spill
traffic is the separate signal that register pressure is the bottleneck.
