#!/usr/bin/env bash
# Regenerate fortad derivative kernels for the fortnum testbed.
#
# fortad is a build-time source generator, not a runtime dependency: the
# generated .f90 files are committed and compile with any conforming Fortran
# compiler, so fortnum gains no new link-time or plugin dependency.
#
# Each kernel here mirrors one fortsym generator, argument for argument, so the
# two derivative kernels can be compared entry by entry. Where fortsym emits a
# tangent-only product, fortad is asked for one too: a signature difference
# would make the comparison a port of the caller rather than of the derivative.
set -euo pipefail

here=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
root=$(cd "$here/../.." && pwd)
fortad_repo=${FORTAD_REPO:-"$root/../fortad"}

fortad_bin=$(find "$fortad_repo/build" -name fortad -type f -perm -u+x 2>/dev/null | head -1)
if [ -z "$fortad_bin" ]; then
    ( cd "$fortad_repo" && fpm build >/dev/null )
    fortad_bin=$(find "$fortad_repo/build" -name fortad -type f -perm -u+x | head -1)
fi

out="$root/src/generated/fortad"
mkdir -p "$out"
kernels="$here/kernels"

# operator : independents : does the fortsym kernel return the value too?
#
# The answer differs per product and per operator, and it is not a detail: a
# fortad kernel that returns the value where fortsym does not is a different
# contract, and the comparison would then be of the caller rather than of the
# derivative.
specs=(
    "det2:a,b,c,d:no:no"
    "det3:a,b,c,d,f,g,h,j,k:no:no"
    "multi_input_p2:x1,x2:yes:yes"
    "multi_input_p4:x1,x2,x3,x4:yes:yes"
    "multi_input_p8:x1,x2,x3,x4,x5,x6,x7,x8:yes:yes"
    "multi_input_p16:x1,x2,x3,x4,x5,x6,x7,x8,x9,x10,x11,x12,x13,x14,x15,x16:yes:yes"
    "lagrange4:x,y1,y2,y3,y4:yes:yes"
    "erf:x:no:no"
    "erfc:x:no:no"
    "sqrt1pm1_stable:x:no:no"
    "sqrt1pm1_raw:x:no:no"
    "dawson_outer_value:f:no:no"
    "toroidal_order:x,current,next_order:no:no"
    "scaled_jacobi_recurrence:x,scale,previous,current:no:no"
    "inv2:a,b,c,d:no:skip"
)

# Operators with more than one intent(out): the dependent has to be named.
# operator : independents : dependent
multi_output=(
    "legendre_recurrence:x,previous,current:next"
    "legendre_recurrence_derivative:x,previous,current:derivative"
    "scalar_root_residual:x,p1,p2:residual"
    "implicit_root_residual:x,p:residual"
    "hypergeom_2f1_term:z,term:next_term"
)

for spec in "${specs[@]}"; do
    IFS=: read -r op indep jvp_value vjp_value <<< "$spec"
    src="$kernels/$op.f90"
    jvp_flags=(--indep "$indep" --name "fortnum_${op}_jvp_fortad")
    [ "$jvp_value" = "no" ] && jvp_flags+=(--no-primal)
    "$fortad_bin" "${jvp_flags[@]}" \
        --module "fortnum_fortad_${op}_jvp" \
        -o "$out/fortnum_${op}_jvp_fortad.f90" "$src"
    # A reverse product needs one scalar dependent. An operator with several
    # outputs and no natural single one - the matrix inverse - is forward only.
    if [ "$vjp_value" != "skip" ]; then
        vjp_flags=(--mode reverse --indep "$indep" --name "fortnum_${op}_vjp_fortad")
        [ "$vjp_value" = "no" ] && vjp_flags+=(--no-primal)
        "$fortad_bin" "${vjp_flags[@]}" \
            --module "fortnum_fortad_${op}_vjp" \
            -o "$out/fortnum_${op}_vjp_fortad.f90" "$src"
    fi
done

for spec in "${multi_output[@]}"; do
    IFS=: read -r op indep dep <<< "$spec"
    src="$kernels/${op%%_derivative}.f90"
    "$fortad_bin" --indep "$indep" --dep "$dep" --no-primal \
        --name "fortnum_${op}_jvp_fortad" \
        --module "fortnum_fortad_${op}_jvp" \
        -o "$out/fortnum_${op}_jvp_fortad.f90" "$src"
    "$fortad_bin" --mode reverse --indep "$indep" --dep "$dep" --no-primal \
        --name "fortnum_${op}_vjp_fortad" \
        --module "fortnum_fortad_${op}_vjp" \
        -o "$out/fortnum_${op}_vjp_fortad.f90" "$src"
done

# The primals themselves are compiled into the library, so the tests have
# something to difference the generated products against. The generator reads
# the copies under tools/; these are what link.
for src in "$kernels"/*.f90; do
    name=$(basename "$src")
    {
        echo "! Copied from tools/fortad/kernels/$name by tools/fortad/generate.sh."
        echo "! It is compiled into the library so the derivative products have a"
        echo "! primal to be differenced against. Edit the original, not this."
        echo
        cat "$src"
    } > "$out/$name"
done

# The original dot_sin testbed kernel, kept as the vector-mode example. Its
# primal is copied into the library so the tests have something to link
# against; the generator still reads the copy that lives beside it here.
{
    echo "! Copied from tools/fortad/dot_sin_kernel.f90 by tools/fortad/generate.sh."
    echo "! It is compiled into the library so the testbed's derivative kernels have a"
    echo "! primal to be checked against. Edit the original, not this."
    echo
    cat "$here/dot_sin_kernel.f90"
} > "$out/fortnum_dot_sin_primal.f90"
# Every kernel is wrapped in a module. fortnum compiles what is in src/, and a
# bare subroutine there gives the caller an unchecked external declaration
# instead of an interface the compiler verifies.
"$fortad_bin" --indep a,b --name fortnum_dot_sin_jvp \
    --module fortnum_fortad_dot_sin_jvp \
    -o "$root/src/generated/fortnum_dot_sin_jvp.f90" "$here/dot_sin_kernel.f90"
"$fortad_bin" --indep a,b --name fortnum_dot_sin_jvp_v -d n_dir \
    --module fortnum_fortad_dot_sin_jvp_v \
    -o "$root/src/generated/fortnum_dot_sin_jvp_v.f90" "$here/dot_sin_kernel.f90"
"$fortad_bin" --mode reverse --indep a,b --name fortnum_dot_sin_vjp \
    --module fortnum_fortad_dot_sin_vjp \
    -o "$root/src/generated/fortnum_dot_sin_vjp.f90" "$here/dot_sin_kernel.f90"

echo "regenerated fortad kernels in src/generated and src/generated/fortad"
