subroutine fortnum_multi_input_p2(x1, x2, value)
    !! The generator's multi-input scalar: sum of sines plus half the square of
    !! the sum. It exists to scale the input count while holding the shape of
    !! the expression fixed, which is what makes the 2-input timing comparable
    !! with the others.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x1, x2
    real(dp), intent(out) :: value
    real(dp) :: total

    value = sin(x1) + sin(x2)
    total = x1 + x2
    value = value + total*total/2.0_dp

end subroutine fortnum_multi_input_p2
