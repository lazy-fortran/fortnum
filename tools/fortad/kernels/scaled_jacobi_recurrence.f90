subroutine fortnum_scaled_jacobi_recurrence(degree, alpha, beta, x, scale, &
                                            previous, current, next)
    !! One step of the scaled Jacobi recurrence. The degree and the two Jacobi
    !! parameters select the polynomial and are inactive.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: degree, alpha, beta, x, scale, previous, current
    real(dp), intent(out) :: next
    real(dp) :: two_degree, lower, upper

    two_degree = degree*2
    lower = alpha + beta + two_degree - 2
    upper = alpha + beta + two_degree
    next = (current*(alpha + beta + two_degree - 1)* &
            (scale*(alpha*alpha - beta*beta) + x*upper*lower) &
            - previous*scale*scale*upper*(alpha + degree - 1)* &
            (beta + degree - 1)*2)/ &
           (degree*(alpha + beta + degree)*lower*2)
end subroutine fortnum_scaled_jacobi_recurrence
