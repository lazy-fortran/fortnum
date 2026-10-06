program test_fortnum_special_riccati_bessel
    ! Riccati-Bessel tables against independent quad-precision oracles:
    !   * Maclaurin series DLMF 10.53.1 (jhat) and 10.53.2 (yhat), used only
    !     where sum|t_k|/|f| stays below 1e17;
    !   * the exact finite Hankel sum DLMF 10.49.6,
    !       x h_l^(1)(x) = (-i)^(l+1) e^(ix) sum_k (i/(2x))^k (l+k)!/(k!(l-k)!),
    !     whose real and imaginary parts are jhat_l and yhat_l, under the
    !     same conditioning guard;
    !   * for l <= x, quad-precision upward recurrence from sin x and
    !     sin x/x - cos x, the stable direction below the turning point
    !     (the implementation uses Miller's backward recurrence whenever
    !     lmax > x, so this is a different algorithm on that branch);
    !   * the Wronskian jhat_{l+1} yhat_l - jhat_l yhat_{l+1} = 1 at every l
    !     (DLMF 10.50.3), which the implementation imposes only at l = 0.
    ! Every jhat outside the turning band |l - x| <= 3 x^(1/3) + 2 must meet
    ! at least one direct oracle; inside the band only the Wronskian applies.
    ! The grid covers the Miller region l > x (including l >> x where jhat
    ! is far below yhat) and the upward region l <= x.
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use fortnum_special, only: riccati_bessel
    implicit none

    integer, parameter :: lmax = 300
    real(dp), parameter :: xs(10) = [1.0e-2_dp, 0.37_dp, 1.0_dp, 3.0_dp, &
        3.14159265358979_dp, 9.5_dp, 25.0_dp, 60.0_dp, 299.5_dp, 1000.25_dp]
    integer, parameter :: lmaxes(4) = [0, 1, 40, lmax]
    real(dp), parameter :: tol = 1.0e-13_dp, tol_wronskian = 1.0e-13_dp
    real(dp) :: jhat(0:lmax), yhat(0:lmax)
    integer :: nfail, i, m, l, top
    real(dp) :: x, ref, w, scale
    logical :: ok, any_ok

    nfail = 0
    do i = 1, size(xs)
        x = xs(i)
        do m = 1, size(lmaxes)
            top = lmaxes(m)
            call riccati_bessel(top, x, jhat(0:top), yhat(0:top))
            do l = 0, top
                if (abs(jhat(l)) < 1.0e-280_dp .or. &
                        abs(yhat(l)) > 1.0e280_dp) cycle
                any_ok = .false.
                call jhat_series(l, x, ref, ok)
                if (ok) call check("jhat series", l, jhat(l), ref)
                any_ok = any_ok .or. ok
                call yhat_series(l, x, ref, ok)
                if (ok) call check("yhat series", l, yhat(l), ref)
                call hankel(l, x, .true., ref, ok)
                if (ok) call check("jhat Hankel", l, jhat(l), ref)
                any_ok = any_ok .or. ok
                if (real(l, dp) <= x) then
                    ref = upward_quad(l, x)
                    call check("jhat upward quad", l, jhat(l), ref)
                    any_ok = .true.
                end if
                call hankel(l, x, .false., ref, ok)
                if (ok) call check("yhat Hankel", l, yhat(l), ref)
                if (l < top .and. abs(yhat(l + 1)) < 1.0e280_dp) then
                    w = jhat(l + 1)*yhat(l) - jhat(l)*yhat(l + 1)
                    scale = abs(jhat(l + 1)*yhat(l)) + abs(jhat(l)*yhat(l + 1))
                    if (.not. abs(w - 1.0_dp) <= tol_wronskian*scale) &
                        call check("Wronskian", l, w, 1.0_dp)
                end if
                if (top == lmax .and. .not. any_ok .and. l <= 2*x + 40 .and. &
                        abs(l - x) > 3.0_dp*x**(1.0_dp/3.0_dp) + 2.0_dp .and. &
                        abs(jhat(l)) > 1.0e-3_dp*hypot(jhat(l), yhat(l))) &
                    call check("no jhat oracle", l, 0.0_dp, 1.0_dp)
            end do
        end do
    end do

    call riccati_bessel(-1, 1.0_dp, jhat(0:0), yhat(0:0))
    if (.not. ieee_is_nan(jhat(0))) call check("invalid lmax", 0, 0.0_dp, 1.0_dp)
    call riccati_bessel(2, 0.0_dp, jhat(0:2), yhat(0:2))
    if (.not. ieee_is_nan(yhat(2))) call check("invalid x", 0, 0.0_dp, 1.0_dp)

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) FAILED"
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine jhat_series(l, xd, value, ok)
        ! x^(l+1)/(2l+1)!! sum_k (-x^2/2)^k/(k! (2l+3)(2l+5)...(2l+2k+1))
        integer, intent(in) :: l
        real(dp), intent(in) :: xd
        real(dp), intent(out) :: value
        logical, intent(out) :: ok
        real(qp) :: xq, term, total, magnitude, prefactor
        integer :: k

        xq = real(xd, qp)
        term = 1.0_qp
        total = 1.0_qp
        magnitude = 1.0_qp
        do k = 1, 2000
            term = -term*xq*xq/(2.0_qp*k*real(2*l + 2*k + 1, qp))
            total = total + term
            magnitude = magnitude + abs(term)
            if (abs(term) < 1.0e-36_qp*magnitude .and. k > xq) exit
        end do
        prefactor = (l + 1)*log(xq) - log_double_factorial(2*l + 1)
        ok = magnitude < 1.0e17_qp*abs(total) .and. prefactor > -640.0_qp
        value = 0.0_dp
        if (ok) value = real(exp(prefactor)*total, dp)
    end subroutine jhat_series

    subroutine yhat_series(l, xd, value, ok)
        ! -(2l-1)!!/x^l sum_k (-x^2/2)^k/(k! (1-2l)(3-2l)...(2k-1-2l))
        integer, intent(in) :: l
        real(dp), intent(in) :: xd
        real(dp), intent(out) :: value
        logical, intent(out) :: ok
        real(qp) :: xq, term, total, magnitude, prefactor
        integer :: k

        xq = real(xd, qp)
        term = 1.0_qp
        total = 1.0_qp
        magnitude = 1.0_qp
        do k = 1, 2000
            term = -term*xq*xq/(2.0_qp*k*real(2*k - 1 - 2*l, qp))
            total = total + term
            magnitude = magnitude + abs(term)
            if (abs(term) < 1.0e-36_qp*magnitude .and. k > xq + l) exit
        end do
        prefactor = log_double_factorial(2*l - 1) - l*log(xq)
        ok = magnitude < 1.0e17_qp*abs(total) .and. prefactor < 640.0_qp
        value = 0.0_dp
        if (ok) value = real(-exp(prefactor)*total, dp)
    end subroutine yhat_series

    subroutine hankel(l, xd, real_part, value, ok)
        integer, intent(in) :: l
        real(dp), intent(in) :: xd
        logical, intent(in) :: real_part
        real(dp), intent(out) :: value
        logical, intent(out) :: ok
        complex(qp) :: term, total, phase
        real(qp) :: xq, magnitude
        integer :: k

        xq = real(xd, qp)
        term = (1.0_qp, 0.0_qp)
        total = term
        magnitude = 1.0_qp
        do k = 1, l
            ! (l+k)!/(k!(l-k)!) ratio (l+k)(l-k+1)/k
            term = term*(0.0_qp, 1.0_qp)/(2.0_qp*xq)* &
                real(l + k, qp)*real(l - k + 1, qp)/real(k, qp)
            total = total + term
            magnitude = magnitude + abs(term)
        end do
        phase = (0.0_qp, -1.0_qp)**(l + 1)*exp(cmplx(0.0_qp, xq, qp))
        total = phase*total
        if (real_part) then
            value = real(real(total), dp)
            ok = magnitude < 1.0e17_qp*abs(real(total))
        else
            value = real(aimag(total), dp)
            ok = magnitude < 1.0e17_qp*abs(aimag(total))
        end if
        ok = ok .and. magnitude < 1.0e280_qp
    end subroutine hankel

    function upward_quad(l, xd) result(value)
        integer, intent(in) :: l
        real(dp), intent(in) :: xd
        real(dp) :: value
        real(qp) :: xq, previous, current, next
        integer :: k

        xq = real(xd, qp)
        previous = sin(xq)
        current = sin(xq)/xq - cos(xq)
        if (l == 0) current = previous
        do k = 1, l - 1
            next = real(2*k + 1, qp)/xq*current - previous
            previous = current
            current = next
        end do
        value = real(current, dp)
    end function upward_quad

    pure real(qp) function log_double_factorial(n) result(s)
        integer, intent(in) :: n
        integer :: k
        s = 0.0_qp
        do k = n, 2, -2
            s = s + log(real(k, qp))
        end do
    end function log_double_factorial

    subroutine check(label, l, got, expected)
        character(*), intent(in) :: label
        integer, intent(in) :: l
        real(dp), intent(in) :: got, expected
        real(dp) :: s
        ! Below the turning point both functions oscillate; relative error
        ! is measured against the Hankel modulus there.
        s = abs(expected)
        if (real(l, dp) + 0.5_dp <= x) s = max(s, hypot(jhat(l), yhat(l)))
        if (.not. abs(got - expected) <= tol*s) then
            nfail = nfail + 1
            write (error_unit, "(a,' l=',i0,' x=',es12.5,2es24.16)") label, l, &
                x, got, expected
        end if
    end subroutine check

end program test_fortnum_special_riccati_bessel
