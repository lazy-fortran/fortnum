subroutine fortnum_implicit_root_residual(x, p, residual)
    !! The residual behind the implicit-function-theorem tangent.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x, p
    real(dp), intent(out) :: residual

    residual = x*x - p
end subroutine fortnum_implicit_root_residual
