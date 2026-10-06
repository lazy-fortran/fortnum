program test_fortnum_special_legendre_q_large_degree
    ! Large-degree reproducer for legendre_q near the cut x = 1.
    !
    ! Upward recurrence lost accuracy for Q_l because Q_l is the minimal
    ! solution of the degree recurrence on x > 1 (relative error 5e-4 at
    ! l = 400, x = 1 + 1/1800).  The oracles below share no recurrence with
    ! the implementation and run in quad precision:
    !   * the Gauss series of DLMF 14.3.7 with 14.3.10,
    !       Q_l(x) = sqrt(pi) l!/(2^(l+1) Gamma(l+3/2)) x^(-l-1)
    !                2F1((l+2)/2, (l+1)/2; l+3/2; 1/x^2),
    !     summed term by term (all terms positive, no cancellation);
    !   * Christoffel formula (Abramowitz-Stegun 8.6.19)
    !       Q_l = P_l Q_0 - sum_{k=1}^{l} P_{k-1} P_{l-k}/k
    !     for x so close to 1 that the series converges too slowly; it is
    !     only used while its cancellation P_l Q_0/|Q_l| stays below 1e12.
    ! Near the cut both recurrence directions are almost neutral, so double
    ! rounding accumulates like l^(3/2) eps (about 3e-12 at l = 2000); the
    ! tolerance 1e-11 sits well below the 6e-7 to 5e-4 upward-recurrence loss.
    ! Derivatives are compared with (l/(x^2-1)) (x Q_l - Q_{l-1}) from the
    ! quad-precision oracle values.
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_special, only: legendre_q, legendre_q_derivative
    implicit none

    integer, parameter :: degrees(14) = [0, 1, 2, 5, 10, 29, 30, 100, 300, &
        400, 700, 1000, 1500, 2000]
    real(qp), parameter :: offsets(6) = [1.0_qp/1800.0_qp, 1.0e-4_qp, &
        1.0e-2_qp, 0.1_qp, 0.5_qp, 1.0_qp]
    real(qp), parameter :: near_offsets(2) = [1.0e-8_qp, 1.0e-6_qp]
    real(dp), parameter :: tol = 1.0e-11_dp, tol_derivative = 2.0e-11_dp
    integer :: nfail, i, j
    real(qp) :: x

    nfail = 0
    do j = 1, size(offsets)
        x = real(1.0_dp + real(offsets(j), dp), qp)
        do i = 1, size(degrees)
            call compare(degrees(i), x, .true.)
        end do
    end do
    do j = 1, size(near_offsets)
        x = real(1.0_dp + real(near_offsets(j), dp), qp)
        do i = 1, size(degrees)
            call compare(degrees(i), x, .false.)
        end do
    end do

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) FAILED"
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    subroutine compare(l, xq, use_series)
        integer, intent(in) :: l
        real(qp), intent(in) :: xq
        logical, intent(in) :: use_series
        real(qp) :: q, qm, dq
        real(dp) :: xd, got, got_derivative
        logical :: usable

        xd = real(xq, dp)
        if (use_series) then
            q = q_series(l, xq)
            qm = 0.0_qp
            if (l > 0) qm = q_series(l - 1, xq)
            usable = .true.
        else
            call q_christoffel(l, xq, q, qm, usable)
        end if
        if (.not. usable) return
        if (l == 0) then
            dq = -1.0_qp/(xq*xq - 1.0_qp)
        else
            dq = real(l, qp)*(xq*q - qm)/(xq*xq - 1.0_qp)
        end if

        got = legendre_q(l, xd)
        got_derivative = legendre_q_derivative(l, xd)
        if (abs(q) < 1.0e-290_qp) then
            if (.not. abs(got) <= 1.0e-280_dp) call fail("underflow", l, xq, &
                got, q)
            return
        end if
        if (.not. abs(got - real(q, dp)) <= tol*abs(real(q, dp))) &
            call fail("Q", l, xq, got, q)
        if (.not. abs(got_derivative - real(dq, dp)) <= &
                tol_derivative*abs(real(dq, dp))) &
            call fail("dQ", l, xq, got_derivative, dq)
    end subroutine compare

    function q_series(l, xq) result(q)
        integer, intent(in) :: l
        real(qp), intent(in) :: xq
        real(qp) :: q
        real(qp) :: a, b, c, z, term, total, prefactor
        integer :: k

        a = real(l + 2, qp)/2.0_qp
        b = real(l + 1, qp)/2.0_qp
        c = real(l, qp) + 1.5_qp
        z = 1.0_qp/(xq*xq)
        term = 1.0_qp
        total = 1.0_qp
        k = 0
        do while (term > 1.0e-36_qp*total)
            term = term*(a + k)*(b + k)/((c + k)*(k + 1))*z
            total = total + term
            k = k + 1
        end do
        prefactor = exp(0.5_qp*log(acos(-1.0_qp)) + log_gamma(real(l + 1, qp)) &
            - log_gamma(c) - real(l + 1, qp)*log(2.0_qp*xq))
        q = prefactor*total
    end function q_series

    subroutine q_christoffel(l, xq, q, qm, usable)
        integer, intent(in) :: l
        real(qp), intent(in) :: xq
        real(qp), intent(out) :: q, qm
        logical, intent(out) :: usable
        real(qp) :: p(0:2000), q0
        integer :: n

        p(0) = 1.0_qp
        p(1) = xq
        do n = 1, l - 1
            p(n + 1) = (real(2*n + 1, qp)*xq*p(n) - real(n, qp)*p(n - 1)) &
                /real(n + 1, qp)
        end do
        q0 = 0.5_qp*log((xq + 1.0_qp)/(xq - 1.0_qp))
        q = christoffel(l, p, q0)
        qm = 0.0_qp
        if (l > 0) qm = christoffel(l - 1, p, q0)
        usable = p(l)*q0 < 1.0e12_qp*abs(q)
    end subroutine q_christoffel

    function christoffel(m, p, q0) result(value)
        integer, intent(in) :: m
        real(qp), intent(in) :: p(0:), q0
        real(qp) :: value
        integer :: k

        value = 0.0_qp
        do k = 1, m
            value = value + p(k - 1)*p(m - k)/real(k, qp)
        end do
        value = p(m)*q0 - value
    end function christoffel

    subroutine fail(label, l, xq, got, expected)
        character(*), intent(in) :: label
        integer, intent(in) :: l
        real(qp), intent(in) :: xq, expected
        real(dp), intent(in) :: got
        nfail = nfail + 1
        write (error_unit, "(a,' l=',i0,' x-1=',es10.3,' got=',es24.16, &
            &' expected=',es24.16,' rel=',es9.2)") label, l, &
            real(xq - 1.0_qp, dp), got, real(expected, dp), &
            real(abs(got - expected)/max(abs(expected), tiny(1.0_qp)), dp)
    end subroutine fail

end program test_fortnum_special_legendre_q_large_degree
