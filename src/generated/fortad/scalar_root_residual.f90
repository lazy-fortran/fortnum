! Copied from tools/fortad/kernels/scalar_root_residual.f90 by tools/fortad/generate.sh.
! It is compiled into the library so the derivative products have a
! primal to be differenced against. Edit the original, not this.

subroutine fortnum_scalar_root_residual(x, p1, p2, residual)
    !! The residual whose root the scalar solver chases.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x, p1, p2
    real(dp), intent(out) :: residual

    residual = x**3 + p1*x - p2
end subroutine fortnum_scalar_root_residual
