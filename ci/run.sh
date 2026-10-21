#!/bin/bash
# The same native test selection runs in a build container or a heavy Slurm tier.
set -euo pipefail
src=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$src"
if [ -n "${SLURM_JOB_ID:-}" ]; then
    root=$(dirname "$src")
    tier=$(cat "$root/tier.txt")
    . "$root/task.env"
    . /home/ert/opt/kin6d-env.sh
    export PATH=/home/ert/.local/bin:$PATH
else
    tier=${1:?usage: ci/run.sh TIER}
    root=$(mktemp -d "${RUNNER_TEMP:-/tmp}/native-ci.XXXXXX")
    export FC=/usr/bin/gfortran CC=/usr/bin/gcc CXX=/usr/bin/g++
    export OPENBLAS_NUM_THREADS=1 MPLBACKEND=Agg
fi
unset DISPLAY
export CMAKE_TOOLCHAIN_FILE=$src/ci/cluster-toolchain.cmake
if [ -d "$root/mirrors" ]; then
    export GIT_CONFIG_GLOBAL=$root/gitconfig
    : > "$GIT_CONFIG_GLOBAL"
    git config --global protocol.file.allow always
    shopt -s nullglob
    for mirror in "$root"/mirrors/*/*.git; do
        slug=${mirror#"$root/mirrors/"}; slug=${slug%.git}
        git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug.git"
        git config --global --add "url.file://$mirror.insteadOf" "https://github.com/$slug"
        git config --global --add "url.file://$mirror.insteadOf" "git@github.com:$slug.git"
    done
fi
export GITHUB_WORKSPACE=$src
export GITHUB_ENV=$root/github-env GITHUB_PATH=$root/github-path GITHUB_OUTPUT=$root/github-output
: > "$GITHUB_ENV"; : > "$GITHUB_PATH"; : > "$GITHUB_OUTPUT"
export SCI_CPUS=${SCI_CPUS:-8}
export CMAKE_BUILD_PARALLEL_LEVEL=$SCI_CPUS CTEST_PARALLEL_LEVEL=$SCI_CPUS
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export RUNNER_TEMP=$root/tmp
mkdir -p "$RUNNER_TEMP" ci-out
{ hostname; /usr/bin/gfortran --version | head -1; git rev-parse HEAD; } | tee ci-out/toolchain.txt
load_job_environment() {
    while IFS='=' read -r key value; do
        [[ $key =~ ^[A-Za-z_][A-Za-z0-9_]*$ ]] && export "$key=$value"
    done < "$GITHUB_ENV"
    while IFS= read -r item; do [ -z "$item" ] || export PATH="$item:$PATH"; done < "$GITHUB_PATH"
}
case "$tier" in
  ci-gfortran)
    echo Configure
    (
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Debug \
  -DFORTNUM_BUILD_TESTING=ON


    )
    load_job_environment
    echo Build
    (
cmake --build build -j "${SCI_CPUS:-8}"

    )
    load_job_environment
    echo Test
    (
ctest --test-dir build --output-on-failure

    )
    load_job_environment
    ;;
  ci-generated)
    echo 'Checkout locked fortsym revision'
    (
set -euo pipefail
git clone https://github.com/lazy-fortran/fortsym.git \
  ../fortsym
tools/codegen/checkout_fortsym.sh ../fortsym

    )
    load_job_environment
    echo 'Register and run CTest drift check'
    (
set -euo pipefail
export PATH="$HOME/.local/bin:$PATH"
# A drift failure should arrive with evidence. The generators log
# which stages each kernel takes and why, which is what distinguishes
# a real change from an environment-dependent classification.
export FORTNUM_CODEGEN_VERBOSE=1
cmake -S . -B build-generated -G Ninja \
  -DFORTNUM_BUILD_TESTING=OFF \
  -DFORTNUM_BUILD_EXAMPLES=OFF \
  -DFORTNUM_CHECK_GENERATED=ON
ctest --test-dir build-generated --output-on-failure \
  -R generated_sources_current

    )
    load_job_environment
    ;;
  ci-oracle)
    python3 -m venv "$root/venv"
    . "$root/venv/bin/activate"
    echo 'Install Python dependencies'
    (
set -euo pipefail
python -m pip install --upgrade pip
python -m pip install numpy scipy

    )
    load_job_environment
    echo 'Generate oracle tables'
    (
set -euo pipefail
# The generator is introduced by #4. Until then, skip cleanly so this
# job is informational only. Search the conventional locations.
gen=""
for cand in \
  tools/oracle/generate.py \
  tools/generate_oracle.py \
  test/oracle/generate.py \
  scripts/generate_oracle.py; do
  if [ -f "$cand" ]; then
    gen="$cand"
    break
  fi
done
if [ -z "$gen" ]; then
  echo "No oracle generator found yet (issue #4); skipping."
  exit 0
fi
echo "Running oracle generator: $gen"
python "$gen"

    )
    load_job_environment
    ;;
  ci-benchmark)
    echo 'Configure benchmark'
    (
cmake -S benchmark -B build-bench -G Ninja \
  -DCMAKE_BUILD_TYPE=Release

    )
    load_job_environment
    echo 'Build benchmark'
    (
cmake --build build-bench -j "${SCI_CPUS:-8}"
    )
    load_job_environment
    echo 'Run benchmark'
    (
./build-bench/bin/bench_main --json | tee bench-run.json
    )
    load_job_environment
    echo 'Gate against baseline'
    (
python3 benchmark/gate.py \
  --baseline benchmark/baseline.json \
  --run bench-run.json \
  --factor 4.0

    )
    load_job_environment
    echo 'Re-select Dawson candidates on hosted hardware'
    (
set -euo pipefail
lscpu > new-target-host.txt
python3 benchmark/collect_dawson_family.py \
  --executable build-bench/bin/bench_dawson_generated_family \
  --output new-target-dawson-family.json

    )
    load_job_environment
    ;;
  *) echo "unknown CI tier: $tier" >&2; exit 2 ;;
esac
