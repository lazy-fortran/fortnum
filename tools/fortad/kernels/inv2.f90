subroutine fortnum_inv2(a, b, c, d, r1, r2, r3, r4)
    !! Entries of the inverse of a 2x2 matrix, row major.
    !!
    !! fortsym generates the inverse's derivative as a rule taking the inverse
    !! entries as input - the identity d(A^-1) = -A^-1 dA A^-1. fortad works
    !! the other way: differentiate the closed-form inverse of the matrix. The
    !! two products are the same map, so the test checks fortad's against the
    !! rule as well as against differences of this primal.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: a, b, c, d
    real(dp), intent(out) :: r1, r2, r3, r4
    real(dp) :: determinant

    determinant = a*d - b*c
    r1 = d/determinant
    r2 = -b/determinant
    r3 = -c/determinant
    r4 = a/determinant
end subroutine fortnum_inv2
