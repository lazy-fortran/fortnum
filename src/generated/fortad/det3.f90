! Copied from tools/fortad/kernels/det3.f90 by tools/fortad/generate.sh.
! It is compiled into the library so the derivative products have a
! primal to be differenced against. Edit the original, not this.

subroutine fortnum_det3(a, b, c, d, f, g, h, j, k, value)
    !! Determinant of a 3x3 matrix, in the generator's cofactor arrangement.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: a, b, c, d, f, g, h, j, k
    real(dp), intent(out) :: value

    value = a*(f*k - j*g) - d*(b*k - j*c) + h*(b*g - f*c)
end subroutine fortnum_det3
