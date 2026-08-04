subroutine fortnum_erf(x, value)
    !! Elementwise erf over an assumed-shape array, matching the shape fortnum's
    !! own erf kernel takes.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value(:)

    value = erf(x)
end subroutine fortnum_erf
