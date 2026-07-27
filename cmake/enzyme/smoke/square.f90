! Enzyme smoke test: reverse-mode derivative of f(x) = x*x.
!
! The differentiable kernel is an ordinary bind(c) Fortran function. The raw
! __enzyme_autodiff entry point is declared only here, inside the test
! wrapper, and never leaks into the public fortnum API (see ad.md sec. 2).
! Conservative subset: a single scalar real(real64) active argument passed by
! value (ad.md sec. 3). The kernel asserts d/dx (x*x) = 2*x against the
! analytic value and exits nonzero on mismatch.
module square_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    private
    public :: square
contains
    pure function square(x) result(y) bind(c, name="fortnum_smoke_square")
        real(c_double), intent(in), value :: x
        real(c_double) :: y
        y = x * x
    end function square
end module square_kernel

program test_square
    use, intrinsic :: iso_c_binding, only: c_double, c_funptr, c_funloc
    use, intrinsic :: iso_fortran_env, only: error_unit
    use square_kernel, only: square
    implicit none

    interface
        function enzyme_autodiff(f, x) result(dx) &
                bind(c, name="__enzyme_autodiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), value :: x
            real(c_double) :: dx
        end function enzyme_autodiff
    end interface

    real(c_double), parameter :: tol = 1.0e-12_c_double
    real(c_double) :: x, dx, analytic

    x = 3.0_c_double
    dx = enzyme_autodiff(c_funloc(square), x)
    analytic = 2.0_c_double * x

    if (abs(dx - analytic) > tol) then
        write (error_unit, "(a,2es24.16)") "FAIL square vjp ", dx, analytic
        stop 1
    end if
    write (*, "(a)") "PASS"
end program test_square
