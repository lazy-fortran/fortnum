program test_fortnum_lapack_gtsv
    ! Tridiagonal solves against the dense LU solve of the same matrix, for
    ! diagonally dominant systems and for systems that need row interchanges,
    ! several right-hand sides, n = 1 and 2, a singular matrix and argument
    ! errors.
    use fortnum_kinds, only: dp
    use fortnum_lapack_gtsv, only: fortnum_dgtsv, fortnum_zgtsv
    use fortnum_lapack_lu, only: fortnum_dgesv, fortnum_zgesv
    implicit none
    real(dp), parameter :: eps = epsilon(1.0_dp)
    integer :: n

    do n = 1, 40, 13
        call real_case(n, 4.0_dp)
        call real_case(n, 0.1_dp)
        call complex_case(n, 4.0_dp)
        call complex_case(n, 0.1_dp)
    end do
    call real_case(2, 0.1_dp)
    call complex_case(2, 0.1_dp)
    call singular()
    call argument_errors()
    print '(a)', 'tridiagonal solves match dense LU'

contains

    ! diag scales the main diagonal: 4 is dominant, 0.1 forces interchanges
    subroutine real_case(n, diag)
        integer, intent(in) :: n
        real(dp), intent(in) :: diag
        real(dp) :: dl(max(1, n - 1)), d(n), du(max(1, n - 1)), a(n, n), b(n, 3), x(n, 3)
        integer :: i, info, ipiv(n)
        a = 0
        do i = 1, n
            d(i) = diag*(1 + 0.3_dp*sin(1.7_dp*i))
            a(i, i) = d(i)
            if (i < n) then
                dl(i) = cos(0.9_dp*i) + 0.5_dp
                du(i) = sin(1.1_dp*i) - 0.7_dp
                a(i + 1, i) = dl(i)
                a(i, i + 1) = du(i)
            end if
            b(i, :) = [cos(0.3_dp*i), real(i, dp), 1.0_dp]
        end do
        x = b
        call fortnum_dgtsv(n, 3, dl, d, du, x, n, info)
        call require(info == 0, 'dgtsv info')
        call fortnum_dgesv(n, 3, a, n, ipiv, b, n, info)
        call require(maxval(abs(x - b)) <= 100*n*eps*maxval(abs(b)), 'dgtsv = dgesv')
    end subroutine real_case

    subroutine complex_case(n, diag)
        integer, intent(in) :: n
        real(dp), intent(in) :: diag
        complex(dp) :: dl(max(1, n - 1)), d(n), du(max(1, n - 1)), a(n, n), b(n, 2), x(n, 2)
        integer :: i, info, ipiv(n)
        a = 0
        do i = 1, n
            d(i) = diag*cmplx(1 + 0.3_dp*sin(1.7_dp*i), 0.2_dp*cos(i*1.0_dp), dp)
            a(i, i) = d(i)
            if (i < n) then
                dl(i) = cmplx(cos(0.9_dp*i) + 0.5_dp, 0.4_dp*sin(0.2_dp*i), dp)
                du(i) = cmplx(sin(1.1_dp*i) - 0.7_dp, -0.3_dp, dp)
                a(i + 1, i) = dl(i)
                a(i, i + 1) = du(i)
            end if
            b(i, :) = [cmplx(cos(0.3_dp*i), 1.0_dp, dp), cmplx(i, -i, dp)]
        end do
        x = b
        call fortnum_zgtsv(n, 2, dl, d, du, x, n, info)
        call require(info == 0, 'zgtsv info')
        call fortnum_zgesv(n, 2, a, n, ipiv, b, n, info)
        call require(maxval(abs(x - b)) <= 100*n*eps*maxval(abs(b)), 'zgtsv = zgesv')
    end subroutine complex_case

    ! rows 2 and 3 equal: the elimination meets a zero pivot
    subroutine singular()
        real(dp) :: dl(3), d(4), du(3), b(4, 1)
        complex(dp) :: zdl(3), zd(4), zdu(3), zb(4, 1)
        integer :: info
        d = [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp]
        dl = [0.0_dp, 0.0_dp, 0.0_dp]
        du = [0.0_dp, 0.0_dp, 0.0_dp]
        d(3) = 0
        b = 1
        call fortnum_dgtsv(4, 1, dl, d, du, b, 4, info)
        call require(info == 3, 'dgtsv singular info')
        zd = 1
        zd(2) = 0
        zdl = 0
        zdu = 0
        zb = 1
        call fortnum_zgtsv(4, 1, zdl, zd, zdu, zb, 4, info)
        call require(info == 2, 'zgtsv singular info')
    end subroutine singular

    subroutine argument_errors()
        real(dp) :: dl(2), d(3), du(2), b(3, 1)
        integer :: info
        d = 1
        dl = 0
        du = 0
        b = 0
        call fortnum_dgtsv(-1, 1, dl, d, du, b, 3, info)
        call require(info == -1, 'bad n')
        call fortnum_dgtsv(3, -1, dl, d, du, b, 3, info)
        call require(info == -2, 'bad nrhs')
        call fortnum_dgtsv(3, 1, dl, d, du, b, 2, info)
        call require(info == -7, 'bad ldb')
    end subroutine argument_errors

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAIL: '//what
            error stop 1
        end if
    end subroutine require

end program test_fortnum_lapack_gtsv
