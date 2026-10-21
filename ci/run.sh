#!/bin/bash
# Execute the repository's original native CI checks in isolated Slurm scratch.
set -euo pipefail
: "${SLURM_JOB_ID:?ci/run.sh must run in Slurm}"
tier=$(cat "$(dirname "$(dirname "$(realpath "$0")")")/../tier.txt")
src=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
root=$(dirname "$src")
cd "$src"
. "$root/task.env"
. /home/ert/opt/kin6d-env.sh
export PATH=/home/ert/.local/bin:$PATH
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
export GITHUB_WORKSPACE=$src
export GITHUB_ENV=$root/github-env GITHUB_PATH=$root/github-path GITHUB_OUTPUT=$root/github-output
: > "$GITHUB_ENV"; : > "$GITHUB_PATH"; : > "$GITHUB_OUTPUT"
export CMAKE_BUILD_PARALLEL_LEVEL=${SCI_CPUS:-4} CTEST_PARALLEL_LEVEL=${SCI_CPUS:-4}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
mkdir -p ci-out
{ hostname; /usr/bin/gfortran --version | head -1; cat "$root/revision.txt"; } | tee ci-out/toolchain.txt
export CMAKE_TOOLCHAIN_FILE=$src/ci/cluster-toolchain.cmake
export RUNNER_TEMP=$root/tmp
mkdir -p "$RUNNER_TEMP"
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
  *) echo "unknown CI tier: $tier" >&2; exit 2 ;;
esac
