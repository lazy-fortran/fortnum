subroutine fortnum_erfc(x, value)
    !! Elementwise erfc over an assumed-shape array, matching the shape fortnum's
    !! own erfc kernel takes.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: value(:)

    value = erfc(x)
end subroutine fortnum_erfc
