program test_fortnum_blas12
    ! Level-1 and level-2 BLAS against direct definitions in quad precision,
    ! with unit, strided and negative increments (a negative increment walks
    ! the vector from its far end, as in reference BLAS).
    use fortnum_kinds, only: dp, sp
    use fortnum_blas12
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(qp), parameter :: tol = 64*real(epsilon(1.0_dp), qp)
    integer, parameter :: incs(4) = [1, 2, -1, -3]
    integer :: ix, iy, n

    do ix = 1, 4
        do iy = 1, 4
            do n = 0, 9, 3
                call vectors(n, incs(ix), incs(iy))
            end do
        end do
    end do
    call norms()
    do ix = 1, 4
        do iy = 1, 4
            call matrix_vector(7, 5, incs(ix), incs(iy))
            call matrix_vector(1, 9, incs(ix), incs(iy))
            call matrix_vector(13, 11, incs(ix), incs(iy))
        end do
    end do
    print '(a)', 'level-1 and level-2 BLAS match quad precision'

contains

    pure integer function pos(i, n, inc)
        integer, intent(in) :: i, n, inc
        pos = merge(1 + (i - 1)*inc, 1 + (n - i)*abs(inc), inc > 0)
    end function pos

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    subroutine vectors(n, incx, incy)
        integer, intent(in) :: n, incx, incy
        real(dp) :: x(40), y(40), y0(40), x0(40)
        complex(dp) :: zx(40), zy(40), zy0(40)
        real(sp) :: sx(40), sy(40)
        real(qp) :: s
        complex(qp) :: zs, zt
        integer :: i
        do i = 1, 40
            x(i) = sin(1.3_dp*i)
            y(i) = cos(0.7_dp*i)
            zx(i) = cmplx(x(i), y(i)**2, dp)
            zy(i) = cmplx(y(i), -x(i), dp)
        end do
        sx = real(x, sp)
        sy = real(y, sp)
        s = 0
        zs = 0
        zt = 0
        do i = 1, n
            s = s + real(x(pos(i, n, incx)), qp)*y(pos(i, n, incy))
            zs = zs + conjg(cmplx(zx(pos(i, n, incx)), kind=qp))*zy(pos(i, n, incy))
            zt = zt + cmplx(zx(pos(i, n, incx)), kind=qp)*zy(pos(i, n, incy))
        end do
        call require(abs(fortnum_ddot(n, x, incx, y, incy) - s) <= tol*(1 + n), 'ddot')
        call require(abs(fortnum_sdot(n, sx, incx, sy, incy) - s) <= 1e-5_qp*(1 + n), 'sdot')
        call require(abs(fortnum_zdotc(n, zx, incx, zy, incy) - zs) <= tol*(1 + n), 'zdotc')
        call require(abs(fortnum_zdotu(n, zx, incx, zy, incy) - zt) <= tol*(1 + n), 'zdotu')
        y0 = y
        call fortnum_daxpy(n, 0.5_dp, x, incx, y, incy)
        do i = 1, n
            call require(abs(y(pos(i, n, incy)) - (y0(pos(i, n, incy)) + 0.5_dp*x(pos(i, n, incx)))) <= tol, 'daxpy')
        end do
        zy0 = zy
        call fortnum_zaxpy(n, (0.5_dp, 1.0_dp), zx, incx, zy, incy)
        do i = 1, n
            call require(abs(zy(pos(i, n, incy)) - (zy0(pos(i, n, incy)) + (0.5_dp, 1.0_dp)*zx(pos(i, n, incx)))) &
                <= tol*4, 'zaxpy')
        end do
        y = y0
        call fortnum_dcopy(n, x, incx, y, incy)
        do i = 1, n
            call require(y(pos(i, n, incy)) == x(pos(i, n, incx)), 'dcopy')
        end do
        x0 = x
        y = y0
        call fortnum_dswap(n, x, incx, y, incy)
        do i = 1, n
            call require(y(pos(i, n, incy)) == x0(pos(i, n, incx)) .and. x(pos(i, n, incx)) == y0(pos(i, n, incy)), &
                'dswap')
        end do
        x = x0
        y = y0
        call fortnum_drot(n, x, incx, y, incy, 0.6_dp, 0.8_dp)
        do i = 1, n
            call require(abs(x(pos(i, n, incx)) - (0.6_dp*x0(pos(i, n, incx)) + 0.8_dp*y0(pos(i, n, incy)))) <= tol &
                .and. abs(y(pos(i, n, incy)) - (0.6_dp*y0(pos(i, n, incy)) - 0.8_dp*x0(pos(i, n, incx)))) <= tol, &
                'drot')
        end do
        if (incx > 0) then
            x = x0
            call fortnum_dscal(n, -2.0_dp, x, incx)
            do i = 1, n
                call require(x(1 + (i - 1)*incx) == -2*x0(1 + (i - 1)*incx), 'dscal')
            end do
            call require(fortnum_idamax(n, x0, incx) == maxloc_first(n, x0, incx), 'idamax')
        end if
    end subroutine vectors

    integer function maxloc_first(n, x, inc)
        integer, intent(in) :: n, inc
        real(dp), intent(in) :: x(*)
        integer :: i
        maxloc_first = 0
        if (n < 1) return
        maxloc_first = 1
        do i = 2, n
            if (abs(x(1 + (i - 1)*inc)) > abs(x(1 + (maxloc_first - 1)*inc))) maxloc_first = i
        end do
    end function maxloc_first

    subroutine norms()
        real(dp) :: x(5)
        complex(dp) :: z(3)
        real(qp) :: want
        x = [3.0_dp, -4.0_dp, 0.0_dp, 12.0_dp, 0.0_dp]
        call require(abs(fortnum_dnrm2(5, x, 1) - 13) <= tol*13, 'dnrm2 plain')
        x = x*1.0e300_dp
        call require(abs(fortnum_dnrm2(5, x, 1)/1.0e300_dp - 13) <= tol*13, 'dnrm2 without overflow')
        x = (x*1.0e-300_dp)*1.0e-300_dp
        want = 13
        call require(abs(fortnum_dnrm2(5, x, 1)/1.0e-300_dp - want) <= 1e-12_qp*want, 'dnrm2 without underflow')
        call require(fortnum_dnrm2(3, x, 2) == fortnum_dnrm2(3, [x(1), x(3), x(5)], 1), 'dnrm2 stride')
        z = [(3.0_dp, 4.0_dp), (0.0_dp, 0.0_dp), (12.0e200_dp, 0.0_dp)]
        call require(abs(fortnum_dznrm2(2, z, 1) - 5) <= tol*5, 'dznrm2')
        call require(abs(fortnum_dznrm2(3, z, 1) - 12.0e200_dp) <= 1e-12_dp*12e200_dp, 'dznrm2 large')
        call require(fortnum_izamax(3, [(1.0_dp, -3.0_dp), (3.0_dp, 1.0_dp), (0.0_dp, 2.0_dp)], 1) == 1, &
            'izamax keeps the first of equal |re|+|im|')
        x = [3.0_dp, -4.0_dp, 0.5_dp, -12.0_dp, 2.0_dp]
        call require(fortnum_dasum(5, x, 1) == 21.5_dp, 'dasum')
        call require(fortnum_dasum(3, x, 2) == 5.5_dp, 'dasum stride')
        call require(fortnum_dasum(2, x, -1) == 0 .and. fortnum_dasum(0, x, 1) == 0, 'dasum empty')
        z = [(3.0_dp, -4.0_dp), (-0.5_dp, 0.0_dp), (1.0_dp, 2.0_dp)]
        call require(fortnum_dzasum(3, z, 1) == 10.5_dp .and. fortnum_dzasum(2, z, 2) == 10, 'dzasum')
    end subroutine norms

    subroutine matrix_vector(m, n, incx, incy)
        integer, intent(in) :: m, n, incx, incy
        real(dp) :: a(m + 2, n), x(60), y(60), y0(60), s(n + 2, n + 2), t(n + 2, n + 2)
        complex(dp) :: za(m + 2, n), zx(60), zy(60), zy0(60), h(n + 2, n + 2)
        real(qp) :: want, ref(n + 2, n + 2)
        complex(qp) :: zwant
        integer :: i, j, k, l, nn
        character :: tr
        do j = 1, n
            do i = 1, m + 2
                a(i, j) = sin(0.3_dp*i + 1.7_dp*j)
                za(i, j) = cmplx(a(i, j), cos(0.9_dp*i*j), dp)
            end do
        end do
        do i = 1, 60
            x(i) = cos(0.41_dp*i)
            y(i) = sin(0.23_dp*i)
            zx(i) = cmplx(x(i), -y(i), dp)
            zy(i) = cmplx(y(i), x(i), dp)
        end do
        y0 = y
        zy0 = zy
        ! dgemv / zgemv, every transpose
        do l = 1, 3
            tr = 'NTC'(l:l)
            if (l < 3) then
                y = y0
                call fortnum_dgemv(tr, m, n, 1.5_dp, a, m + 2, x, incx, -0.5_dp, y, incy)
                do i = 1, merge(m, n, tr == 'N')
                    want = -0.5_qp*y0(pos(i, merge(m, n, tr == 'N'), incy))
                    do j = 1, merge(n, m, tr == 'N')
                        if (tr == 'N') then
                            want = want + 1.5_qp*a(i, j)*x(pos(j, n, incx))
                        else
                            want = want + 1.5_qp*a(j, i)*x(pos(j, m, incx))
                        end if
                    end do
                    call require(abs(y(pos(i, merge(m, n, tr == 'N'), incy)) - want) <= tol*(m + n), 'dgemv '//tr)
                end do
            end if
            zy = zy0
            call fortnum_zgemv(tr, m, n, (0.5_dp, 1.0_dp), za, m + 2, zx, incx, (0.0_dp, 0.0_dp), zy, incy)
            do i = 1, merge(m, n, tr == 'N')
                zwant = 0
                do j = 1, merge(n, m, tr == 'N')
                    select case (tr)
                    case ('N')
                        zwant = zwant + za(i, j)*cmplx(zx(pos(j, n, incx)), kind=qp)
                    case ('T')
                        zwant = zwant + za(j, i)*cmplx(zx(pos(j, m, incx)), kind=qp)
                    case ('C')
                        zwant = zwant + conjg(cmplx(za(j, i), kind=qp))*zx(pos(j, m, incx))
                    end select
                end do
                zwant = (0.5_qp, 1.0_qp)*zwant
                call require(abs(zy(pos(i, merge(m, n, tr == 'N'), incy)) - zwant) <= 4*tol*(m + n), 'zgemv '//tr)
            end do
        end do
        ! dsymv, zhemv, dtrsv, ztrsv on an n x n matrix built from a
        nn = n
        do j = 1, nn
            do i = 1, nn
                s(i, j) = sin(0.5_dp*i*j) + merge(4.0_dp, 0.0_dp, i == j)
                h(i, j) = cmplx(s(i, j), 0.3_dp*(i - j), dp)
            end do
        end do
        do l = 1, 2
            y = y0
            call fortnum_dsymv('UL'(l:l), nn, 2.0_dp, s, nn + 2, x, incx, 1.0_dp, y, incy)
            zy = zy0
            call fortnum_zhemv('UL'(l:l), nn, (1.0_dp, -1.0_dp), h, nn + 2, zx, incx, (0.0_dp, 0.0_dp), zy, incy)
            do i = 1, nn
                want = y0(pos(i, nn, incy))
                zwant = 0
                do j = 1, nn
                    want = want + 2*sym(s, i, j, l == 1)*x(pos(j, nn, incx))
                    zwant = zwant + her(h, i, j, l == 1)*cmplx(zx(pos(j, nn, incx)), kind=qp)
                end do
                call require(abs(y(pos(i, nn, incy)) - want) <= tol*nn*8, 'dsymv')
                call require(abs(zy(pos(i, nn, incy)) - (1.0_qp, -1.0_qp)*zwant) <= tol*nn*8, 'zhemv')
            end do
            ! triangular solves: check op(T) x_new = x_old
            do j = 1, 2
                tr = 'NT'(j:j)
                t = s
                y = x
                call fortnum_dtrsv('UL'(l:l), tr, 'N', nn, t, nn + 2, y, incx)
                ref = 0
                do i = 1, nn
                    want = 0
                    do k = 1, n
                        want = want + tri(t, i, k, l == 1, tr == 'T')*y(pos(k, n, incx))
                    end do
                    call require(abs(want - x(pos(i, n, incx))) <= tol*n*16, 'dtrsv')
                end do
            end do
        end do
    end subroutine matrix_vector

    real(qp) function sym(s, i, j, upper)
        real(dp), intent(in) :: s(:, :)
        integer, intent(in) :: i, j
        logical, intent(in) :: upper
        if ((i <= j) .eqv. upper .or. i == j) then
            sym = s(i, j)
        else
            sym = s(j, i)
        end if
    end function sym

    complex(qp) function her(h, i, j, upper)
        complex(dp), intent(in) :: h(:, :)
        integer, intent(in) :: i, j
        logical, intent(in) :: upper
        if (i == j) then
            her = real(h(i, i), qp)
        else if ((i < j) .eqv. upper) then
            her = h(i, j)
        else
            her = conjg(cmplx(h(j, i), kind=qp))
        end if
    end function her

    real(qp) function tri(t, i, j, upper, trans)
        real(dp), intent(in) :: t(:, :)
        integer, intent(in) :: i, j
        logical, intent(in) :: upper, trans
        integer :: r, c
        r = merge(j, i, trans)
        c = merge(i, j, trans)
        tri = 0
        if ((upper .and. r <= c) .or. (.not. upper .and. r >= c)) tri = t(r, c)
    end function tri

end program test_fortnum_blas12
