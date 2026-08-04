#!/usr/bin/env bash
# Regenerate fortad derivative kernels for the fortnum testbed.
#
# fortad is a build-time source generator, not a runtime dependency: the
# generated .f90 files are committed and compile with any conforming Fortran
# compiler, so fortnum gains no new link-time or plugin dependency.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root=$(cd "$here/../.." && pwd)
fortad_repo=${FORTAD_REPO:-"$root/../fortad"}

fortad_bin=$(find "$fortad_repo/build" -name fortad -type f -perm -u+x 2>/dev/null | head -1)
if [ -z "$fortad_bin" ]; then
    ( cd "$fortad_repo" && fpm build >/dev/null )
    fortad_bin=$(find "$fortad_repo/build" -name fortad -type f -perm -u+x | head -1)
fi

out="$root/src/generated"
mkdir -p "$out"

"$fortad_bin" --indep a,b --name fortnum_dot_sin_jvp \
    -o "$out/fortnum_dot_sin_jvp.f90" "$here/dot_sin_kernel.f90"
"$fortad_bin" --indep a,b --name fortnum_dot_sin_jvp_v -d n_dir \
    -o "$out/fortnum_dot_sin_jvp_v.f90" "$here/dot_sin_kernel.f90"
"$fortad_bin" --mode reverse --indep a,b --name fortnum_dot_sin_vjp \
    -o "$out/fortnum_dot_sin_vjp.f90" "$here/dot_sin_kernel.f90"

echo "regenerated fortad kernels in src/generated"
