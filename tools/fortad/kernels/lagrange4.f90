subroutine fortnum_lagrange4(x, y1, y2, y3, y4, value)
    !! Lagrange interpolation through four samples at the fixed nodes
    !! -1, 0, 1, 2. The nodes are constants, so only the sample values and the
    !! evaluation point are active.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x, y1, y2, y3, y4
    real(dp), intent(out) :: value
    real(dp) :: b1, b2, b3, b4

    b1 = (x - 0.0_dp)*(x - 1.0_dp)*(x - 2.0_dp)/((-1.0_dp)*(-2.0_dp)*(-3.0_dp))
    b2 = (x + 1.0_dp)*(x - 1.0_dp)*(x - 2.0_dp)/((1.0_dp)*(-1.0_dp)*(-2.0_dp))
    b3 = (x + 1.0_dp)*(x - 0.0_dp)*(x - 2.0_dp)/((2.0_dp)*(1.0_dp)*(-1.0_dp))
    b4 = (x + 1.0_dp)*(x - 0.0_dp)*(x - 1.0_dp)/((3.0_dp)*(2.0_dp)*(1.0_dp))
    value = y1*b1 + y2*b2 + y3*b3 + y4*b4
end subroutine fortnum_lagrange4
