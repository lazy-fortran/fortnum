! Copied from tools/fortad/kernels/implicit_root_residual.f90 by tools/fortad/generate.sh.
! It is compiled into the library so the derivative products have a
! primal to be differenced against. Edit the original, not this.

subroutine fortnum_implicit_root_residual(x, p, residual)
    !! The residual behind the implicit-function-theorem tangent.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x, p
    real(dp), intent(out) :: residual

    residual = x*x - p
end subroutine fortnum_implicit_root_residual
