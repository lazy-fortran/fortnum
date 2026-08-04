! Copied from tools/fortad/kernels/sqrt1pm1_raw.f90 by tools/fortad/generate.sh.
! It is compiled into the library so the derivative products have a
! primal to be differenced against. Edit the original, not this.

subroutine fortnum_sqrt1pm1_raw(x, value)
    !! sqrt(1+x) - 1 written directly. Kept beside the stable form because the
    !! pair is what the accuracy comparison is about.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x
    real(dp), intent(out) :: value

    value = sqrt(x + 1) - 1
end subroutine fortnum_sqrt1pm1_raw
