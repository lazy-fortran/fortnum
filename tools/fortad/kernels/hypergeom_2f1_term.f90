subroutine fortnum_hypergeom_2f1_term(a, b, c, k, z, term, next_term)
    !! One step of the 2F1 series recurrence. The parameters and the term index
    !! select the series; the argument and the running term are active.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: a, b, c, k, z, term
    real(dp), intent(out) :: next_term

    next_term = term*z*(a + k)*(b + k)/((c + k)*(k + 1))
end subroutine fortnum_hypergeom_2f1_term
