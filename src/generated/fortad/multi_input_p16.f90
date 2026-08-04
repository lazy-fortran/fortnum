! Copied from tools/fortad/kernels/multi_input_p16.f90 by tools/fortad/generate.sh.
! It is compiled into the library so the derivative products have a
! primal to be differenced against. Edit the original, not this.

subroutine fortnum_multi_input_p16(x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16, value)
    !! The generator's multi-input scalar: sum of sines plus half the square of
    !! the sum. It exists to scale the input count while holding the shape of
    !! the expression fixed, which is what makes the 16-input timing comparable
    !! with the others.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11, x12, x13, x14, x15, x16
    real(dp), intent(out) :: value
    real(dp) :: total

    value = sin(x1) + sin(x2) + sin(x3) + sin(x4) + sin(x5) + sin(x6) + sin(x7) + &
        sin(x8) + sin(x9) + sin(x10) + sin(x11) + sin(x12) + sin(x13) + sin(x14) + &
        sin(x15) + sin(x16)
    total = x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11 + x12 + x13 + &
        x14 + x15 + x16
    value = value + total*total/2.0_dp

end subroutine fortnum_multi_input_p16
