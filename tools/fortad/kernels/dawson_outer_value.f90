subroutine fortnum_dawson_outer_value(f, value)
    !! The outer function applied to a Dawson evaluation.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: f
    real(dp), intent(out) :: value

    value = sin(f) + f**2
end subroutine fortnum_dawson_outer_value
