module ode_residual_certificate_fixture
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use fortnum_interval, only: interval_t, interval, operator(+), &
        operator(-), operator(*)
    use fortnum_ode_cash_karp, only: cash_karp_step
    implicit none

    !> Linear test system y' = a y with exact binary64 matrix entries.
    type :: linear_system_t
        real(dp), allocatable :: a(:, :)
        logical :: fail = .false.
    end type linear_system_t

contains

    subroutine linear_rhs(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx

        dydt = 0.0_dp
        if (.not. present(ctx)) return
        select type (ctx)
        type is (linear_system_t)
            dydt = matmul(ctx%a, y)
        end select
        if (t < -1.0_dp) dydt = 0.0_dp
    end subroutine linear_rhs

    !> Interval action of the point matrix, independent of the candidate.
    subroutine linear_action(x, ax, ok, ctx)
        type(interval_t), intent(in) :: x(:)
        type(interval_t), intent(out) :: ax(:)
        logical, intent(out) :: ok
        class(*), intent(in), optional :: ctx
        integer :: i, j

        ok = .false.
        if (.not. present(ctx)) return
        select type (ctx)
        type is (linear_system_t)
            if (ctx%fail) return
            do i = 1, size(x)
                ax(i) = interval(0.0_dp)
            end do
            do j = 1, size(x)
                do i = 1, size(x)
                    ax(i) = ax(i) + ctx%a(i, j)*x(j)
                end do
            end do
            ok = .true.
        end select
    end subroutine linear_action

    subroutine logistic_rhs(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx

        dydt = y*(1.0_dp - y)
        if (t < -1.0_dp .or. present(ctx)) dydt = dydt
    end subroutine logistic_rhs

    subroutine logistic_box(t, x, fx, ok, ctx)
        type(interval_t), intent(in) :: t, x(:)
        type(interval_t), intent(out) :: fx(:)
        logical, intent(out) :: ok
        class(*), intent(in), optional :: ctx

        fx = x*(1.0_dp - x)
        ok = t%hi >= t%lo .or. present(ctx)
    end subroutine logistic_box

    !> Uniform-step Cash-Karp candidate: nodes t, values y, slopes f.
    subroutine make_candidate(rhs, y0, duration, steps, t, y, f, ctx)
        interface
            subroutine rhs(t, y, dydt, ctx)
                import :: dp
                real(dp), intent(in) :: t, y(:)
                real(dp), intent(out) :: dydt(:)
                class(*), intent(in), optional :: ctx
            end subroutine rhs
        end interface
        real(dp), intent(in) :: y0(:), duration
        integer, intent(in) :: steps
        real(dp), allocatable, intent(out) :: t(:), y(:, :), f(:, :)
        class(*), intent(in), optional :: ctx
        real(dp), allocatable :: k(:, :)
        real(dp) :: h
        integer :: j, nfev

        allocate (t(0:steps), y(size(y0), 0:steps), f(size(y0), 0:steps))
        allocate (k(size(y0), 9))
        h = duration/real(steps, dp)
        y(:, 0) = y0
        nfev = 0
        do j = 0, steps
            t(j) = real(j, dp)*h
        end do
        call rhs(t(0), y(:, 0), f(:, 0), ctx)
        do j = 1, steps
            call cash_karp_step(rhs, t(j - 1), y(:, j - 1), t(j) - t(j - 1), &
                .false., k(:, 1), k(:, 2), k(:, 3), k(:, 4), k(:, 5), k(:, 6), &
                k(:, 7), k(:, 8), k(:, 9), nfev, ctx)
            y(:, j) = k(:, 8)
            call rhs(t(j), y(:, j), f(:, j), ctx)
        end do
    end subroutine make_candidate

    !> exp(t a) y0 in binary128 by scaling and squaring of a Taylor series.
    function expm_apply_qp(a, t, y0) result(y)
        real(dp), intent(in) :: a(:, :), y0(:)
        real(qp), intent(in) :: t
        real(qp) :: y(size(y0))
        real(qp) :: e(size(y0), size(y0)), term(size(y0), size(y0))
        real(qp) :: b(size(y0), size(y0))
        integer :: i, k, s

        b = t*real(a, qp)
        s = 0
        do while (maxval(sum(abs(b), dim=1)) > 0.25_qp)
            b = b/2.0_qp
            s = s + 1
        end do
        e = 0.0_qp
        term = 0.0_qp
        do i = 1, size(y0)
            e(i, i) = 1.0_qp
            term(i, i) = 1.0_qp
        end do
        do k = 1, 40
            term = matmul(term, b)/real(k, qp)
            e = e + term
        end do
        do k = 1, s
            e = matmul(e, e)
        end do
        y = matmul(e, real(y0, qp))
    end function expm_apply_qp

    !> Hermite basis evaluation in binary128 (independent of the
    !> production power-coefficient form).
    function hermite_qp(y0, y1, f0, f1, ta, tb, t) result(p)
        real(dp), intent(in) :: y0(:), y1(:), f0(:), f1(:), ta, tb
        real(qp), intent(in) :: t
        real(qp) :: p(size(y0)), h, s, h00, h10, h01, h11

        h = real(tb, qp) - real(ta, qp)
        s = (t - real(ta, qp))/h
        h00 = (1 + 2*s)*(1 - s)**2
        h10 = s*(1 - s)**2
        h01 = s**2*(3 - 2*s)
        h11 = s**2*(s - 1)
        p = h00*real(y0, qp) + h10*h*real(f0, qp) + h01*real(y1, qp) &
            + h11*h*real(f1, qp)
    end function hermite_qp

end module ode_residual_certificate_fixture

program test_fortnum_ode_residual_certificate
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf, ieee_is_finite
    use fortnum_interval, only: interval_t, interval, sqrt, operator(-), &
        operator(/)
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_ode_residual_certificate, only: certify_linear_trajectory, &
        certify_nonlinear_trajectory
    use ode_residual_certificate_fixture
    implicit none

    integer, parameter :: oracle_expm = 1, oracle_oscillator = 2, &
        oracle_logistic = 3
    character(len=32) :: arg

    call test_unitary_complex()
    call test_dissipative()
    call test_oscillator()
    call test_logistic()
    call test_corruption()
    call test_failures()
    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        if (arg == '--benchmark') call benchmark_unitary(125, 20)
    end if
    print '(a)', 'test_fortnum_ode_residual_certificate: all checks passed'

contains

    !> Real/imaginary split of -i H for Hermitian H = hr + i hi.
    subroutine hermitian_generator(n, sys)
        integer, intent(in) :: n
        type(linear_system_t), intent(out) :: sys
        real(dp) :: hr(n, n), hi(n, n)
        integer :: i, j

        do j = 1, n
            do i = 1, n
                hr(i, j) = sin(1.3_dp*real(i + j, dp)) + merge(1.0_dp, 0.0_dp, &
                                                               i == j)*real(i, dp)/n
                hi(i, j) = cos(0.7_dp*real(i*j, dp))*real(i - j, dp)/n
            end do
        end do
        hi = 0.5_dp*(hi - transpose(hi))
        hr = 0.5_dp*(hr + transpose(hr))
        allocate (sys%a(2*n, 2*n))
        sys%a(1:n, 1:n) = hi
        sys%a(1:n, n + 1:2*n) = hr
        sys%a(n + 1:2*n, 1:n) = -hr
        sys%a(n + 1:2*n, n + 1:2*n) = hi
    end subroutine hermitian_generator

    function exact_qp(kind, sys, y0, t) result(y)
        integer, intent(in) :: kind
        type(linear_system_t), intent(in) :: sys
        real(dp), intent(in) :: y0(:)
        real(qp), intent(in) :: t
        real(qp) :: y(size(y0))

        select case (kind)
        case (oracle_expm)
            y = expm_apply_qp(sys%a, t, y0)
        case (oracle_oscillator)
            y(1) = y0(1)*cos(t) + y0(2)*sin(t)
            y(2) = -y0(1)*sin(t) + y0(2)*cos(t)
        case default
            y(1) = 1.0_qp/(1.0_qp + (1.0_qp/real(y0(1), qp) - 1.0_qp)*exp(-t))
        end select
    end function exact_qp

    !> ||y(t) - p(t)|| <= radius(k) at nodes and interior samples.
    subroutine check_containment(kind, sys, y0, t, y, f, radius, label)
        integer, intent(in) :: kind
        type(linear_system_t), intent(in) :: sys
        real(dp), intent(in) :: y0(:), t(0:), y(:, 0:), f(:, 0:), radius(0:)
        character(*), intent(in) :: label
        real(qp), parameter :: frac(5) = [0.0_qp, 0.137_qp, 0.5_qp, &
                                          0.81_qp, 1.0_qp]
        real(qp) :: tt, err
        integer :: k, j

        if (radius(0) /= 0.0_dp) error stop 'initial radius must be zero'
        do k = 1, size(t) - 1
            do j = 1, size(frac)
                tt = real(t(k - 1), qp) + frac(j)*(real(t(k), qp) - real(t(k - 1), qp))
                err = norm2(exact_qp(kind, sys, y0, tt) - hermite_qp(y(:, k - 1), &
                    y(:, k), f(:, k - 1), f(:, k), t(k - 1), t(k), tt))
                if (err > real(radius(k), qp)) then
                    print *, label, k, j, real(err), radius(k)
                    error stop 'exact solution escapes certified radius'
                end if
            end do
        end do
    end subroutine check_containment

    subroutine run_linear(sys, y0, duration, steps, mu, kind, label, &
                          final_radius, max_residual)
        type(linear_system_t), intent(in) :: sys
        real(dp), intent(in) :: y0(:), duration
        integer, intent(in) :: steps, kind
        type(interval_t), intent(in) :: mu
        character(*), intent(in) :: label
        real(dp), intent(out) :: final_radius, max_residual
        real(dp), allocatable :: t(:), y(:, :), f(:, :), residual(:), radius(:)
        type(interval_t), allocatable :: c(:, :), w(:, :)
        type(fortnum_status_t) :: status

        call make_candidate(linear_rhs, y0, duration, steps, t, y, f, sys)
        allocate (c(size(y0), 0:3), w(size(y0), 0:3), residual(steps), &
                  radius(0:steps))
        call certify_linear_trajectory(linear_action, t, y, f, mu, 0.0_dp, &
                                       c, w, residual, radius, status, sys)
        if (status%code /= FORTNUM_OK) then
            print *, label, trim(status%msg)
            error stop 'linear certificate failed'
        end if
        call check_containment(kind, sys, y0, t, y, f, radius, label)
        final_radius = radius(steps)
        max_residual = maxval(residual)
    end subroutine run_linear

    subroutine test_unitary_complex()
        type(linear_system_t) :: sys
        real(dp) :: y0(6), r1, r2, m1, m2

        call hermitian_generator(3, sys)
        y0 = [0.6_dp, -0.2_dp, 0.1_dp, 0.3_dp, 0.0_dp, -0.7_dp]
        call run_linear(sys, y0, 3.0_dp, 40, interval(0.0_dp), oracle_expm, &
                        'unitary40', r1, m1)
        call run_linear(sys, y0, 3.0_dp, 80, interval(0.0_dp), oracle_expm, &
                        'unitary80', r2, m2)
        ! r = dp/dt - A p is O(h^3) per cell (h R is O(h^4)); halving h
        ! divides the residual and the accumulated radius by about 8.
        if (m1/m2 < 6.0_dp .or. m1/m2 > 10.0_dp) error stop 'residual order lost'
        if (r1/r2 < 6.0_dp .or. r1/r2 > 10.0_dp) error stop 'radius order lost'
        if (r2 > 5.0e-4_dp) error stop 'unitary radius needlessly broad'
    end subroutine test_unitary_complex

    subroutine test_dissipative()
        type(linear_system_t) :: sys
        type(interval_t) :: mu
        real(dp) :: r1, m1, r0, m0

        ! Non-normal decaying A: mu_2(A) = (sqrt(17) - 3)/2 > 0 exercises the
        ! Gronwall branch although every eigenvalue is negative.
        allocate (sys%a(2, 2))
        sys%a = reshape([-1.0_dp, 0.0_dp, 4.0_dp, -2.0_dp], [2, 2])
        mu = (sqrt(interval(17.0_dp)) - 3.0_dp)/2.0_dp
        call run_linear(sys, [1.0_dp, 0.5_dp], 2.0_dp, 50, mu, oracle_expm, &
                        'dissipative', r1, m1)
        call run_linear(sys, [1.0_dp, 0.5_dp], 2.0_dp, 50, interval(0.0_dp), &
                        oracle_expm, 'dissipative-mu0', r0, m0)
        if (.not. (r1 > r0)) error stop 'Gronwall factor not applied'
        if (r1 > 1.0e-3_dp) error stop 'dissipative radius needlessly broad'
    end subroutine test_dissipative

    subroutine test_oscillator()
        type(linear_system_t) :: sys
        real(dp) :: r1, m1

        allocate (sys%a(2, 2))
        sys%a = reshape([0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
        call run_linear(sys, [1.0_dp, 0.25_dp], 20.0_dp, 200, interval(0.0_dp), &
                        oracle_oscillator, 'oscillator', r1, m1)
        if (r1 > 1.0e-3_dp) error stop 'oscillator radius needlessly broad'
    end subroutine test_oscillator

    subroutine test_logistic()
        type(linear_system_t) :: none
        real(dp), allocatable :: t(:), y(:, :), f(:, :), residual(:), radius(:)
        type(interval_t) :: c(1, 0:3), w(1, 0:2)
        type(fortnum_status_t) :: status
        real(dp) :: coarse
        integer :: steps, pass

        coarse = 0.0_dp
        do pass = 1, 2
            steps = 60*pass
            call make_candidate(logistic_rhs, [0.1_dp], 3.0_dp, steps, t, y, f)
            if (allocated(residual)) deallocate (residual, radius)
            allocate (residual(steps), radius(0:steps))
            ! mu = 1 bounds f'(y) = 1 - 2y for y >= 0; checked a posteriori:
            ! the solution stays >= 0.1 and every radius is below 0.1.
            call certify_nonlinear_trajectory(logistic_box, t, y, f, 32, &
                interval(1.0_dp), 0.0_dp, c, w, residual, radius, status)
            if (status%code /= FORTNUM_OK) error stop 'logistic certificate failed'
            if (maxval(radius) >= 0.1_dp) error stop 'logistic mu premise fails'
            call check_containment(oracle_logistic, none, [0.1_dp], t, y, f, &
                                   radius, 'logistic')
            if (pass == 2 .and. radius(steps) > 0.6_dp*coarse) &
                error stop 'logistic first-order refinement lost'
            coarse = radius(steps)
        end do
    end subroutine test_logistic

    subroutine test_corruption()
        type(linear_system_t) :: sys
        real(dp), allocatable :: t(:), y(:, :), f(:, :), residual(:), radius(:)
        type(interval_t) :: c(6, 0:3), w(6, 0:3)
        type(fortnum_status_t) :: status
        real(dp) :: y0(6), clean
        integer :: pass

        call hermitian_generator(3, sys)
        y0 = [0.6_dp, -0.2_dp, 0.1_dp, 0.3_dp, 0.0_dp, -0.7_dp]
        allocate (residual(80), radius(0:80))
        do pass = 0, 2
            call make_candidate(linear_rhs, y0, 3.0_dp, 80, t, y, f, sys)
            if (pass == 1) y(2, 40) = y(2, 40) + 1.0e-3_dp
            if (pass == 2) f(5, 10) = f(5, 10) - 1.0e-1_dp
            call certify_linear_trajectory(linear_action, t, y, f, &
                interval(0.0_dp), 0.0_dp, c, w, residual, radius, status, sys)
            if (status%code /= FORTNUM_OK) error stop 'finite corruption rejected'
            call check_containment(oracle_expm, sys, y0, t, y, f, radius, 'corrupt')
            if (pass == 0) clean = radius(80)
            ! A target ten times the clean radius must reject corrupted traces.
            if (pass > 0 .and. radius(80) < 10.0_dp*clean) &
                error stop 'corrupted trace certified as accurate'
        end do
    end subroutine test_corruption

    subroutine expect_failure(status, radius, label)
        type(fortnum_status_t), intent(in) :: status
        real(dp), intent(in) :: radius(0:)
        character(*), intent(in) :: label

        if (status%code == FORTNUM_OK) then
            print *, label
            error stop 'invalid input accepted'
        end if
        if (any(ieee_is_finite(radius(1:)))) then
            print *, label
            error stop 'failed certificate left finite radii'
        end if
    end subroutine expect_failure

    subroutine test_failures()
        type(linear_system_t) :: sys
        real(dp), allocatable :: t(:), y(:, :), f(:, :), ts(:)
        real(dp) :: residual(10), radius(0:10), short(9)
        type(interval_t) :: c(2, 0:3), w(2, 0:3)
        type(fortnum_status_t) :: status

        allocate (sys%a(2, 2))
        sys%a = reshape([0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
        call make_candidate(linear_rhs, [1.0_dp, 0.0_dp], 1.0_dp, 10, t, y, f, sys)
        y(1, 3) = ieee_value(0.0_dp, ieee_quiet_nan)
        call certify_linear_trajectory(linear_action, t, y, f, interval(0.0_dp), &
            0.0_dp, c, w, residual, radius, status, sys)
        call expect_failure(status, radius, 'nan')
        call make_candidate(linear_rhs, [1.0_dp, 0.0_dp], 1.0_dp, 10, t, y, f, sys)
        ts = t
        ts(5) = ts(4)
        call certify_linear_trajectory(linear_action, ts, y, f, interval(0.0_dp), &
            0.0_dp, c, w, residual, radius, status, sys)
        call expect_failure(status, radius, 'mesh')
        call certify_linear_trajectory(linear_action, t, y, f, &
            interval(ieee_value(0.0_dp, ieee_positive_inf)), 0.0_dp, c, w, &
            residual, radius, status, sys)
        call expect_failure(status, radius, 'mu')
        call certify_linear_trajectory(linear_action, t, y, f, interval(0.0_dp), &
            -1.0_dp, c, w, residual, radius, status, sys)
        call expect_failure(status, radius, 'initial')
        call certify_linear_trajectory(linear_action, t, y, f, interval(0.0_dp), &
            0.0_dp, c, w, short, radius, status, sys)
        if (status%code == FORTNUM_OK) error stop 'shape mismatch accepted'
        call certify_linear_trajectory(linear_action, t, y, f, interval(0.0_dp), &
            0.0_dp, c, w, residual, radius, status)
        call expect_failure(status, radius, 'missing context')
        sys%fail = .true.
        call certify_linear_trajectory(linear_action, t, y, f, interval(0.0_dp), &
            0.0_dp, c, w, residual, radius, status, sys)
        call expect_failure(status, radius, 'action failure')
        call certify_nonlinear_trajectory(logistic_box, t, y(1:1, :), &
            f(1:1, :), 0, interval(1.0_dp), 0.0_dp, c(1:1, :), w(1:1, 0:2), &
            residual, radius, status)
        call expect_failure(status, radius, 'nsub')
    end subroutine test_failures

    !> Warm per-cell cost of the strict certificate for a 125-dimensional
    !> complex unitary system (250 real unknowns, interval dense action).
    subroutine benchmark_unitary(n, steps)
        integer, intent(in) :: n, steps
        type(linear_system_t) :: sys
        real(dp), allocatable :: t(:), y(:, :), f(:, :), residual(:), radius(:)
        real(dp), allocatable :: y0(:)
        type(interval_t), allocatable :: c(:, :), w(:, :)
        type(fortnum_status_t) :: status
        integer(int64) :: c0, c1, rate
        real(dp) :: best, seconds
        integer :: trial, j

        call hermitian_generator(n, sys)
        allocate (y0(2*n))
        do j = 1, 2*n
            y0(j) = cos(0.37_dp*real(j, dp))/sqrt(real(n, dp))
        end do
        call make_candidate(linear_rhs, y0, 0.02_dp, steps, t, y, f, sys)
        allocate (c(2*n, 0:3), w(2*n, 0:3), residual(steps), radius(0:steps))
        best = huge(1.0_dp)
        do trial = 1, 5
            call system_clock(c0, rate)
            call certify_linear_trajectory(linear_action, t, y, f, &
                interval(0.0_dp), 0.0_dp, c, w, residual, radius, status, sys)
            call system_clock(c1)
            if (status%code /= FORTNUM_OK) error stop 'benchmark certificate failed'
            seconds = real(c1 - c0, dp)/real(rate, dp)
            best = min(best, seconds)
        end do
        print '(a,i0,a,i0,a,es10.3,a,es10.3)', 'benchmark n=', n, ' steps=', &
            steps, ' best seconds per cell=', best/steps, ' final radius=', &
            radius(steps)
    end subroutine benchmark_unitary

end program test_fortnum_ode_residual_certificate
