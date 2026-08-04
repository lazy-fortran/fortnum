subroutine fortnum_legendre_recurrence(degree, order, x, previous, current, &
                                       next, derivative)
    !! One step of the associated Legendre recurrence, with the derivative
    !! relation alongside it. The degree and order select the polynomial and
    !! are inactive; the argument and the two previous values are active.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: degree, order, x, previous, current
    real(dp), intent(out) :: next, derivative

    next = (current*x*(degree*2 + 1) - previous*(degree + order))/ &
           (degree - order + 1)
    derivative = (current*degree*x - previous*(degree + order))/(x*x - 1)
end subroutine fortnum_legendre_recurrence
