# Benchmarks

`benchmark/` contains three kinds of performance evidence:

| Location | Purpose |
| --- | --- |
| `benchmark/reference/*.json` | Machine-readable measurements and selections |
| `benchmark/report/data/*.csv` | Normalized data for cross-kernel reports |
| `benchmark/report/app/*.f90` | `fortplot` figure generators |

Generated figures stay outside the repository.

## Primal microbenchmarks

Build the standalone harness in Release mode:

```bash
cmake -S benchmark -B build-bench -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-bench
./build-bench/bin/bench_main
```

JSON output can be checked against the committed baseline:

```bash
./build-bench/bin/bench_main --json | python3 benchmark/gate.py
python3 benchmark/gate.py --run run.json --factor 1.5
```

Refresh `benchmark/baseline.json` only for an intentional performance change
measured on a quiet reference host:

```bash
./build-bench/bin/bench_main --json > benchmark/baseline.json
```

Commit the new baseline with the implementation that changed it.

The weekly `reference benchmarks` workflow runs on the self-hosted runner
labeled `ryzen9-5950x`. It verifies the CPU model, pins execution to CPU 4,
checks the committed baseline, and retains raw JSON, peak-RSS, cache-counter,
host, compiler, and source-revision evidence as a workflow artifact. It never
rewrites or commits the baseline automatically.

The hosted CI runner also re-selects the fused and separate Dawson JVP and VJP
candidates. Reproduce its machine-readable wall-clock, peak-RSS, validation,
and cache-counter record with:

```bash
python3 benchmark/collect_dawson_family.py \
  --executable build-bench/bin/bench_dawson_generated_family \
  --output dawson-family.json
```

For the downstream `itpplasma/SIMPLE` comparison, build the same SIMPLE
revision once with this checkout through
`FETCHCONTENT_SOURCE_DIR_FORTNUM` and once with SIMPLE's pinned fortnum
revision. Then run:

```bash
benchmark/itpplasma_simple_benchmark.sh \
  /path/to/current/simple.x \
  /path/to/pinned/simple.x \
  /path/to/SIMPLE/test/golden_record/test_tokamak_classifier/simple.in \
  /tmp/simple-fortnum-benchmark
```

The script pins one CPU, runs three warmups and 15 samples, verifies
byte-identical deterministic application outputs, and records complete-process
wall clock, peak RSS, and cache counters.

## Derivative tournaments

Derivative benchmarks compare implementations of the same mathematical
product and workload. Use these public labels:

- `analytical`
- `autodiff`
- `hybrid`
- `finite_difference_reference` for diagnostics

Each record must state:

- operator, product, active inputs, outputs, and derivative directions
- compiler, flags, source revision, hardware, and affinity
- workload repetitions, warmups, samples, median, and dispersion
- validation error and oracle
- candidate-specific peak memory
- reusable primal values, traces, factors, or preconditioners
- selected candidate and deterministic reason

Report complete-workload wall clock first. Add CPU work and cache counters,
native code size, transfer-inclusive GPU time, resident GPU time, device
memory, bandwidth, occupancy, registers, and spills when available and
relevant.

“Faster” always compares rows with the same operator, product, workload, and
returned values. Never compare an isolated derivative kernel with a
value-plus-derivative candidate without labeling the difference.

The scalar FFT benchmark runs through CTest so validation and timing use the
same executable:

```bash
cmake -S benchmark -B build-bench -G Ninja -DCMAKE_BUILD_TYPE=Release
cmake --build build-bench --target bench_fft_scalar
FORTNUM_FFT_SCALAR_ACTION=benchmark \
FORTNUM_FFT_SCALAR_PRODUCT=jvp \
FORTNUM_FFT_SCALAR_LENGTH=1024 \
ctest --test-dir build-bench -R '^fft_scalar_validation$' -V
```

## Reports

Generate the cumulative CPU report:

```bash
output_dir=$(mktemp -d)
cd benchmark/report
fo exec plot_differentiation_report \
  data/mechanism_tournaments.csv \
  data/ode_parameter_crossover.csv \
  "${output_dir}"
```

Generate the GPU report:

```bash
output_dir=$(mktemp -d)
cd benchmark/report
fo exec plot_gpu_report \
  data/gpu_batch_scaling.csv \
  data/gpu_active_scaling.csv \
  data/gpu_profile.csv \
  "${output_dir}"
```

The report programs use the pinned `fortplot` dependency from
`benchmark/report/fpm.toml`.

## Add a benchmark

Keep the benchmark workload observable so the optimizer cannot delete it.
Measure candidates in separate processes when peak RSS must be attributed to
one candidate. Pin the process and record the affinity for low-latency CPU
workloads.

For a new derivative tournament:

1. implement an independent oracle
2. validate every candidate before timing
3. choose representative workload classes
4. record raw or machine-readable results
5. normalize report rows when the tournament meets the report inclusion rule
6. update the evidence index and selected-candidate registry where applicable

See [docs/design/differentiation_benchmarks.md](../docs/design/differentiation_benchmarks.md)
for the evidence catalog and
[docs/design/differentiation_report.md](../docs/design/differentiation_report.md)
for aggregate statistics.
