subroutine fortnum_toroidal_order(degree, order, x, current, next_order, &
                                  following_order)
    !! One step of the toroidal-harmonic order recurrence.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: degree, order, x, current, next_order
    real(dp), intent(out) :: following_order

    following_order = current*(degree + order + 1)*(degree - order) &
                      - next_order*x*(order + 1)*2/sqrt(x*x - 1)
end subroutine fortnum_toroidal_order
