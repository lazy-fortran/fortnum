#!/usr/bin/env bash
# Independent behavioral oracle for the N_machine collector.
#
# Compiles a kernel with a known pure-arithmetic expression and asserts the
# disassembly-derived counts match values derived by hand from that
# expression (not by re-running the collector):
#
#     jvp = a*vd - b*vc - c*vb + d*va
#
# N_emit:  4 multiplies + 3 add/subtract = 7 instructions, 7 FLOPs.
# N_machine with FMA (x86-64-v3) fuses each mul with the following add/sub
# into one FMA, so 7 instructions collapse to 4 while the FLOP count is
# preserved at 7.  That is the acceptance criterion for #77.
#
# The test is opt-in: when objdump or a Fortran compiler is absent it exits
# cleanly, because absence of the binary utility is "level not measured" and
# must never be reported as a skipped or zeroed measurement.
set -euo pipefail

source_root=${1:?repository root is required}
collector="${source_root}/benchmark/collect_machine_count.py"
kernel="${source_root}/src/generated/fortnum_det2_jvp_kernel.f90"

command -v gfortran >/dev/null 2>&1 || { echo "gfortran absent; N_machine not measured"; exit 0; }
command -v objdump >/dev/null 2>&1 || { echo "objdump absent; N_machine not measured"; exit 0; }

temp_dir=$(mktemp -d)
trap 'rm -rf "${temp_dir}"' EXIT
object="${temp_dir}/det2.o"

# FMA fusion on x86-64-v3 is the reduction under test; it also needs the
# compiler to fuse mul+add (gfortran enables contraction by default).
gfortran -O3 -march=x86-64-v3 -c -o "${object}" "${kernel}"

record=$(python3 "${collector}" \
    --object "${object}" --source "${kernel}" --target cpu)
if ! command -v python3 >/dev/null 2>&1; then
    echo "python3 absent; N_machine not measured" >&2
    exit 1
fi

# Independent oracle values, derived from the expression itself.
read_json() {
    python3 -c "import json,sys; d=json.loads(sys.argv[1]); print(d$1)" "$record"
}

emit_inst=$(read_json "['n_emit']['instruction_count']")
emit_flop=$(read_json "['n_emit']['flop_count']")
mach_inst=$(read_json "['n_machine']['instruction_count']")
mach_flop=$(read_json "['n_machine']['flop_count']")
fma=$(read_json "['n_machine']['fma_count']")
loads=$(read_json "['n_machine']['spill_loads']")
stores=$(read_json "['n_machine']['spill_stores']")

fail=0
[ "$emit_inst" = "7" ]  || { echo "N_emit instructions: got $emit_inst, expected 7"; fail=1; }
[ "$emit_flop" = "7" ]  || { echo "N_emit flops: got $emit_flop, expected 7"; fail=1; }
[ "$mach_inst" = "4" ]  || { echo "N_machine instructions: got $mach_inst, expected 4 (FMA fusion)"; fail=1; }
[ "$mach_flop" = "7" ]  || { echo "N_machine flops: got $mach_flop, expected 7 (preserved)"; fail=1; }
[ "$fma" = "3" ]        || { echo "FMA count: got $fma, expected 3"; fail=1; }
# det2 has no register pressure, so there are no spills.
[ "$loads" = "0" ]      || { echo "spill loads: got $loads, expected 0"; fail=1; }
[ "$stores" = "0" ]     || { echo "spill stores: got $stores, expected 0"; fail=1; }

if [ "$fail" -ne 0 ]; then
    echo "N_machine collector failed its independent oracle" >&2
    echo "$record" >&2
    exit 1
fi
echo "N_machine collector oracle passed (N_emit 7 -> N_machine 4, FLOPs 7 preserved)"
