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
# The same expression in real64 CUDA, when nvcc and cuobjdump are present,
# must count 1 DMUL + 3 DFMA: 4 instructions, 7 FLOPs, 3 FMAs.
#
# Each part is opt-in: when its compiler or binary utility is absent, or the
# host is not x86-64 for the x86-64-v3 CPU oracle, it is not run, because an
# absent level is "not measured" and must never be reported as a skipped or
# zeroed measurement.
set -euo pipefail

source_root=${1:?repository root is required}
collector="${source_root}/benchmark/collect_machine_count.py"
kernel="${source_root}/src/generated/fortnum_det2_jvp_kernel.f90"

command -v python3 >/dev/null 2>&1 || { echo "python3 absent; N_machine not measured"; exit 0; }

temp_dir=$(mktemp -d)
trap 'rm -rf "${temp_dir}"' EXIT
fail=0
measured=0

field() {
    python3 -c "import json,sys; d=json.loads(sys.argv[1]); print(d$2)" "$1"
}

expect() {
    local label=$1 got=$2 want=$3
    [ "$got" = "$want" ] || { echo "$label: got $got, expected $want"; fail=1; }
}

# Both oracles share the det2 JVP: 4 instructions, 7 FLOPs, 3 FMAs and no
# spills, because det2 has no register pressure.
expect_machine() {
    local record=$1
    expect "N_machine instructions" "$(field "$record" "['n_machine']['instruction_count']")" 4
    expect "N_machine flops" "$(field "$record" "['n_machine']['flop_count']")" 7
    expect "FMA count" "$(field "$record" "['n_machine']['fma_count']")" 3
    expect "spill loads" "$(field "$record" "['n_machine']['spill_loads']")" 0
    expect "spill stores" "$(field "$record" "['n_machine']['spill_stores']")" 0
}

case "$(uname -m)" in
    x86_64 | amd64) cpu_host=1 ;;
    *) cpu_host=0 ;;
esac

if [ "$cpu_host" -ne 1 ]; then
    echo "host is not x86-64; x86-64-v3 CPU N_machine not measured"
elif ! command -v gfortran >/dev/null 2>&1; then
    echo "gfortran absent; CPU N_machine not measured"
elif ! command -v objdump >/dev/null 2>&1; then
    echo "objdump absent; CPU N_machine not measured"
else
    # FMA fusion on x86-64-v3 is the reduction under test; it also needs the
    # compiler to fuse mul+add (gfortran enables contraction by default).
    gfortran -O3 -march=x86-64-v3 -c -o "${temp_dir}/det2.o" "${kernel}"
    record=$(python3 "${collector}" \
        --object "${temp_dir}/det2.o" --source "${kernel}" --target cpu)
    expect "N_emit instructions" "$(field "$record" "['n_emit']['instruction_count']")" 7
    expect "N_emit flops" "$(field "$record" "['n_emit']['flop_count']")" 7
    expect_machine "$record"
    [ "$fail" -eq 0 ] || echo "CPU record: $record" >&2
    measured=1
fi

if ! command -v nvcc >/dev/null 2>&1; then
    echo "nvcc absent; CUDA N_machine not measured"
elif ! command -v cuobjdump >/dev/null 2>&1; then
    echo "cuobjdump absent; CUDA N_machine not measured"
else
    cat > "${temp_dir}/det2.cu" <<'CUDA'
extern "C" __global__ void det2_jvp(const double *a, const double *b,
    const double *c, const double *d, const double *va, const double *vb,
    const double *vc, const double *vd, double *jvp, int n)
{
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) {
        jvp[i] = a[i] * vd[i] - b[i] * vc[i] - c[i] * vb[i] + d[i] * va[i];
    }
}
CUDA
    nvcc -O3 -c -o "${temp_dir}/det2_cuda.o" "${temp_dir}/det2.cu"
    record=$(python3 "${collector}" \
        --object "${temp_dir}/det2_cuda.o" --target cuda)
    expect_machine "$record"
    [ "$fail" -eq 0 ] || echo "CUDA record: $record" >&2
    measured=1
fi

if [ "$fail" -ne 0 ]; then
    echo "N_machine collector failed its independent oracle" >&2
    exit 1
fi
if [ "$measured" -eq 1 ]; then
    echo "N_machine collector oracle passed (N_emit 7 -> N_machine 4, FLOPs 7 preserved)"
fi
