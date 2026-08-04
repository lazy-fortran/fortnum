! Copied from tools/fortad/kernels/erfc.f90 by tools/fortad/generate.sh.
! It is compiled into the library so the derivative products have a
! primal to be differenced against. Edit the original, not this.

subroutine fortnum_erfc(x, value)
    !! Elementwise erfc over an assumed-shape array, matching the shape fortnum's
    !! own erfc kernel takes.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value(:)

    value = erfc(x)
end subroutine fortnum_erfc
