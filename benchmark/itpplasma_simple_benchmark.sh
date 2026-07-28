#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -ne 4 ]; then
    echo "usage: $0 CURRENT_SIMPLE PINNED_SIMPLE SIMPLE_IN OUTPUT_DIR" >&2
    exit 2
fi

current_exe=$1
pinned_exe=$2
config=$3
output=$4
samples=${FORTNUM_SIMPLE_SAMPLES:-15}
cpu=${FORTNUM_SIMPLE_CPU:-4}

mkdir -p "$output"

run_candidate() {
    label=$1
    exe=$2
    run_dir="$output/$label"
    times="$output/${label}_times_ms.txt"

    mkdir -p "$run_dir"
    cp "$config" "$run_dir/simple.in"
    : > "$times"

    (
        cd "$run_dir"
        for _ in 1 2 3; do
            OMP_NUM_THREADS=1 taskset -c "$cpu" "$exe" >/dev/null 2>&1
        done
        for sample in $(seq 1 "$samples"); do
            start_ns=$(date +%s%N)
            OMP_NUM_THREADS=1 taskset -c "$cpu" "$exe" \
                > "sample_${sample}.log" 2>&1
            end_ns=$(date +%s%N)
            awk -v first="$start_ns" -v last="$end_ns" \
                'BEGIN {printf "%.6f\n", (last-first)/1e6}' >> "$times"
        done
        OMP_NUM_THREADS=1 /usr/bin/time -v \
            -o "$output/${label}_memory.txt" \
            taskset -c "$cpu" "$exe" >/dev/null 2>&1
        OMP_NUM_THREADS=1 perf stat \
            -o "$output/${label}_cache.txt" \
            -e cycles,instructions,cache-references,cache-misses \
            taskset -c "$cpu" "$exe" >/dev/null 2>&1
        sha256sum class_parts.dat confined_fraction.dat \
            orbit_exit_code.dat times_lost.dat \
            | sort -k2 > "$output/${label}_hashes.txt"
    )
}

run_candidate current "$current_exe"
run_candidate pinned "$pinned_exe"
diff -u "$output/pinned_hashes.txt" "$output/current_hashes.txt"

printf 'candidate,median_ms,mad_ms,peak_rss_kib\n' > "$output/summary.csv"
for label in current pinned; do
    median=$(sort -n "$output/${label}_times_ms.txt" \
        | awk -v middle=$(((samples + 1)/2)) 'NR==middle {print $1}')
    mad=$(awk -v center="$median" \
        '{delta=$1-center; if(delta<0)delta=-delta; print delta}' \
        "$output/${label}_times_ms.txt" | sort -n \
        | awk -v middle=$(((samples + 1)/2)) 'NR==middle {print $1}')
    rss=$(awk '/Maximum resident set size/ {print $6}' \
        "$output/${label}_memory.txt")
    printf '%s,%s,%s,%s\n' "$label" "$median" "$mad" "$rss" \
        >> "$output/summary.csv"
done

cat "$output/summary.csv"
