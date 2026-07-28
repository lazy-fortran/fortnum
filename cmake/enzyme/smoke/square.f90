module square_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    private
    public :: square
contains
    pure function square(x) result(y) bind(c, name="fortnum_smoke_square")
        real(c_double), intent(in), value :: x
        real(c_double) :: y

        y = x*x
    end function square
end module square_kernel

program test_square
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_generated_enzyme_square, only: &
        fortnum_enzyme_square_vjp_scalar
    implicit none

    real(c_double), parameter :: tolerance = 1.0e-12_c_double
    real(c_double) :: cotangent, got, reference, x
    integer :: i

    do i = 1, 257
        x = -4.0_c_double + 8.0_c_double*real(i - 1, c_double)/256.0_c_double
        cotangent = -0.7_c_double + 0.003_c_double*real(i, c_double)
        got = fortnum_enzyme_square_vjp_scalar(x, cotangent)
        reference = 2.0_c_double*x*cotangent
        if (abs(got - reference) > tolerance*max(1.0_c_double, &
            abs(reference))) then
            error stop "square autodiff VJP mismatch"
        end if
    end do
    write (*, "(a)") "PASS"

end program test_square
