! Enzyme smoke test: reverse-mode VJP of a weighted array sum/dot.
!
! The kernel computes s = sum(w*x) over a contiguous explicit-shape array
! (ad.md sec. 3). Reverse mode accumulates the gradient of s into the shadow
! array dx; with a unit output seed the VJP is d s / d x_i = w_i. The weight
! array w is inactive (its shadow dw is supplied but the kernel does not write
! it). The raw __enzyme_autodiff symbol stays inside this wrapper.
module asum_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    implicit none
    private
    public :: wsum
contains
    pure function wsum(x, w, n) result(s) bind(c, name="fortnum_smoke_wsum")
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n), w(n)
        real(c_double) :: s
        integer :: i
        s = 0.0_c_double
        do i = 1, n
            s = s + w(i) * x(i)
        end do
    end function wsum
end module asum_kernel

program test_array_sum
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_funptr, c_funloc
    use, intrinsic :: iso_fortran_env, only: error_unit
    use asum_kernel, only: wsum
    implicit none

    interface
        function enzyme_autodiff(f, x, dx, w, dw, n) result(s) &
                bind(c, name="__enzyme_autodiff")
            import :: c_double, c_int, c_funptr
            type(c_funptr), value :: f
            real(c_double), intent(in)    :: x(*), w(*), dw(*)
            real(c_double), intent(inout) :: dx(*)
            integer(c_int), value :: n
            real(c_double) :: s
        end function enzyme_autodiff
    end interface

    integer(c_int), parameter :: n = 5
    real(c_double), parameter :: tol = 1.0e-12_c_double
    real(c_double) :: x(n), w(n), dx(n), dw(n), s
    integer :: i

    do i = 1, n
        x(i) = real(i, c_double)
        w(i) = real(2 * i, c_double)
    end do
    dx = 0.0_c_double
    dw = 0.0_c_double

    s = enzyme_autodiff(c_funloc(wsum), x, dx, w, dw, n)

    do i = 1, n
        if (abs(dx(i) - w(i)) > tol) then
            write (error_unit, "(a,i0,2es24.16)") "FAIL sum vjp ", i, dx(i), w(i)
            stop 1
        end if
    end do
    write (*, "(a)") "PASS"
end program test_array_sum
