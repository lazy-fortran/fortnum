#!/bin/bash
# Bundle the exact checkout and public Git dependencies before Slurm dispatch.
set -euo pipefail
out=$(realpath -m "${1:?usage: ci/bundle.sh OUT.tar.zst}")
src=$(git rev-parse --show-toplevel)
cache=${SCI_CACHE:-${HOME}/.cache/sci-bundle}
stage=$(mktemp -d "${TMPDIR:-/tmp}/sci-bundle.XXXXXX")
trap 'rm -rf "$stage"' EXIT
mkdir -p "$cache/mirrors" "$stage/mirrors"
while read -r slug; do
    [[ -z $slug || $slug == \#* ]] && continue
    mkdir -p "$cache/locks"
    exec {lock}> "$cache/locks/${slug//\//_}"
    flock "$lock"
    mirror=$cache/mirrors/$slug.git
    if [ ! -d "$mirror" ]; then
        mkdir -p "$(dirname "$mirror")"
        git clone -q --mirror "https://github.com/$slug.git" "$mirror"
    else
        git -C "$mirror" fetch -q --prune origin '+refs/heads/*:refs/heads/*' '+refs/tags/*:refs/tags/*'
    fi
    mkdir -p "$stage/mirrors/${slug%/*}"
    cp -a "$mirror" "$stage/mirrors/$slug.git"
    exec {lock}>&-
done < "$src/ci/dependencies.txt"
git clone -q --no-hardlinks --dissociate "$src" "$stage/src"
git -C "$stage/src" remote set-url origin "https://github.com/${GITHUB_REPOSITORY:?GITHUB_REPOSITORY missing}.git"
git -C "$src" rev-parse HEAD > "$stage/revision.txt"
printf '%s\n' "${SCI_TIER:?SCI_TIER missing}" > "$stage/tier.txt"
: > "$stage/task.env"
while IFS= read -r key; do
    printf 'export %s=%q\n' "$key" "${!key}" >> "$stage/task.env"
done < <(compgen -A variable SCI_MATRIX_)
tar -C "$stage" -cf - src mirrors revision.txt tier.txt task.env | zstd -q -T2 -o "$out"
sha256sum "$out"
