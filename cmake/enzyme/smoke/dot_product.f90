! Enzyme smoke test: forward-mode JVP of a dot product.
!
! The kernel computes s = sum(x*y) over a contiguous explicit-shape array
! (ad.md sec. 3). Forward mode seeds the array tangents dx, dy and the
! returned directional derivative must match the analytic JVP
! sum(dx*y + x*dy). The raw __enzyme_fwddiff symbol stays inside this wrapper.
module dot_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_int
    implicit none
    private
    public :: dotp
contains
    pure function dotp(x, y, n) result(s) bind(c, name="fortnum_smoke_dotp")
        integer(c_int), intent(in), value :: n
        real(c_double), intent(in) :: x(n), y(n)
        real(c_double) :: s
        integer :: i
        s = 0.0_c_double
        do i = 1, n
            s = s + x(i) * y(i)
        end do
    end function dotp
end module dot_kernel

program test_dot
    use, intrinsic :: iso_c_binding, only: c_double, c_int, c_funptr, c_funloc
    use, intrinsic :: iso_fortran_env, only: error_unit
    use dot_kernel, only: dotp
    implicit none

    interface
        function enzyme_fwddiff(f, x, dx, y, dy, n) result(ds) &
                bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_int, c_funptr
            type(c_funptr), value :: f
            real(c_double), intent(in) :: x(*), dx(*), y(*), dy(*)
            integer(c_int), value :: n
            real(c_double) :: ds
        end function enzyme_fwddiff
    end interface

    integer(c_int), parameter :: n = 4
    real(c_double), parameter :: tol = 1.0e-12_c_double
    real(c_double) :: x(n), y(n), dx(n), dy(n)
    real(c_double) :: ds, analytic
    integer :: i

    x = [1.0_c_double, 2.0_c_double, 3.0_c_double, 4.0_c_double]
    y = [5.0_c_double, 6.0_c_double, 7.0_c_double, 8.0_c_double]
    dx = [1.0_c_double, 0.0_c_double, 0.5_c_double, 0.0_c_double]
    dy = [0.0_c_double, 2.0_c_double, 0.0_c_double, 0.0_c_double]

    ds = enzyme_fwddiff(c_funloc(dotp), x, dx, y, dy, n)

    analytic = 0.0_c_double
    do i = 1, n
        analytic = analytic + dx(i) * y(i) + x(i) * dy(i)
    end do

    if (abs(ds - analytic) > tol) then
        write (error_unit, "(a,2es24.16)") "FAIL dot jvp ", ds, analytic
        stop 1
    end if
    write (*, "(a)") "PASS"
end program test_dot
