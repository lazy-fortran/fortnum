program test_fortnum_cubic_spline
    ! Behavioural tests for fortnum_cubic_spline:
    !   1. A cubic polynomial on nonuniform nodes is reproduced exactly
    !      (value and derivative) between and at the nodes, for every column.
    !   2. A smooth non-polynomial converges at fourth order in the value.
    !   3. Invalid input is refused with its status.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_cubic_spline, only: cubic_spline_t, cubic_spline_init, &
        cubic_spline_eval
    implicit none
    integer :: nfail
    nfail = 0
    call test_cubic_exact(nfail)
    call test_convergence(nfail)
    call test_refusal(nfail)
    if (nfail > 0) error stop 'test_fortnum_cubic_spline: FAIL'
    print '(a)', 'test_fortnum_cubic_spline: all passed'
contains

    subroutine test_cubic_exact(nfail)
        integer, intent(inout) :: nfail
        type(cubic_spline_t) :: sp
        real(dp) :: x(7), y(7, 2), v(2), d(2), xq, err
        integer :: i, status
        x = [0.0_dp, 0.3_dp, 0.45_dp, 1.0_dp, 1.7_dp, 1.75_dp, 2.6_dp]
        y(:, 1) = poly(x)
        y(:, 2) = 2 - x
        call cubic_spline_init(x, y, sp, status)
        if (status /= 0) then
            nfail = nfail + 1
            return
        end if
        err = 0
        do i = 0, 40
            xq = 2.6_dp*i/40
            call cubic_spline_eval(sp, xq, v, d)
            err = max(err, abs(v(1) - poly(xq)), abs(d(1) - dpoly(xq)), &
                abs(v(2) - (2 - xq)), abs(d(2) + 1))
        end do
        if (err > 1e-13_dp) then
            print '(a,es10.2)', 'FAIL cubic not reproduced: ', err
            nfail = nfail + 1
        end if
    end subroutine test_cubic_exact

    subroutine test_convergence(nfail)
        integer, intent(inout) :: nfail
        type(cubic_spline_t) :: sp
        real(dp), allocatable :: x(:), y(:, :)
        real(dp) :: v(1), d(1), err(2), xq
        integer :: k, n, i, status
        do k = 1, 2
            n = 10*2**k
            x = [(real(i, dp)/(n - 1), i=0, n - 1)]
            allocate (y(n, 1))
            y(:, 1) = sin(3*x)
            call cubic_spline_init(x, y, sp, status)
            err(k) = 0
            do i = 0, 200
                xq = real(i, dp)/200
                call cubic_spline_eval(sp, xq, v, d)
                err(k) = max(err(k), abs(v(1) - sin(3*xq)))
            end do
            deallocate (y)
        end do
        if (log(err(1)/err(2))/log(2.0_dp) < 3.5_dp) then
            print '(a,2es10.2)', 'FAIL value convergence below fourth order: ', err
            nfail = nfail + 1
        end if
    end subroutine test_convergence

    subroutine test_refusal(nfail)
        integer, intent(inout) :: nfail
        type(cubic_spline_t) :: sp
        real(dp) :: y3(3, 1), y4(4, 1), y5(5, 1)
        integer :: status
        y3 = 0; y4 = 0; y5 = 0
        call cubic_spline_init([0.0_dp, 1.0_dp, 2.0_dp], y3, sp, status)
        if (status /= 1) nfail = nfail + 1
        call cubic_spline_init([0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], y4, sp, status)
        if (status /= 2) nfail = nfail + 1
        call cubic_spline_init([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], y5, sp, status)
        if (status /= 3) nfail = nfail + 1
    end subroutine test_refusal

    elemental real(dp) function poly(x)
        real(dp), intent(in) :: x
        poly = 1.5_dp - 2*x + 0.7_dp*x**2 - 0.4_dp*x**3
    end function poly

    elemental real(dp) function dpoly(x)
        real(dp), intent(in) :: x
        dpoly = -2 + 1.4_dp*x - 1.2_dp*x**2
    end function dpoly
end program test_fortnum_cubic_spline
