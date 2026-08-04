! Copied from tools/fortad/kernels/det2.f90 by tools/fortad/generate.sh.
! It is compiled into the library so the derivative products have a
! primal to be differenced against. Edit the original, not this.

subroutine fortnum_det2(a, b, c, d, value)
    !! Determinant of a 2x2 matrix, row major.
    !!
    !! The primal for the fortad path. It states the same expression the
    !! fortsym generator states symbolically, so the two derivative kernels are
    !! comparable entry by entry rather than only in aggregate.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: a, b, c, d
    real(dp), intent(out) :: value

    value = a*d - c*b
end subroutine fortnum_det2
