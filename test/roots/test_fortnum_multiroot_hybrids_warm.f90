program test_fortnum_multiroot_hybrids_warm
    ! Warm start of multiroot_hybrids from the previous solve's Jacobian.
    !
    ! A sequence of Broyden tridiagonal systems (More, Garbow & Hillstrom,
    ! ACM TOMS 7 (1981) 17, problem 30) with a slowly drifting constant term,
    ! each started from the previous root as an implicit time stepper would.
    ! Oracles: the residual |F|_inf at every returned point, agreement of the
    ! warm and cold roots (the root is isolated, so both must find it), and
    ! the evaluation count, which must drop clearly below the cold count, also
    ! with a residual noise floor above ftol.  A garbage warm estimate must
    ! still converge (refresh after poor steps), a mismatched size is
    ! ignored, and only non-finite factors are dropped after a failure.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_multiroot, only: multiroot_hybrids, multiroot_hybrids_warm_t
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    integer, parameter :: N = 8, NSTEP = 40
    real(dp), parameter :: FTOL = 1.0e-12_dp, DRIFT = 2.0e-3_dp
    integer, parameter :: NO_ROOT = 2, NAN_RESIDUAL = 3
    integer(int64), parameter :: MASK52 = 2_int64**52 - 1
    integer :: nfev, nfail
    real(dp) :: noise

    nfail = 0
    call check_sequence(0.0_dp)
    ! A residual noise floor above ftol, as when a caller's tolerance is below
    ! the rounding level of its residual: every solve ends stagnating there.
    call check_sequence(1.0e-9_dp)
    noise = 0.0_dp
    call check_garbage()
    call check_size_mismatch()
    call check_failure_paths()

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " check(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine check_sequence(noise_level)
        real(dp), intent(in) :: noise_level
        real(dp) :: xc(N), xw(N), xc_next(N), xw_next(N), f(N), shift
        type(fortnum_status_t) :: sc, sw
        type(multiroot_hybrids_warm_t) :: warm
        integer :: k, cold_total, warm_total, max_ratio
        real(dp) :: max_dx

        noise = noise_level
        max_ratio = merge(6, 9, noise == 0.0_dp)
        xc = -1.0_dp
        xw = -1.0_dp
        cold_total = 0
        warm_total = 0
        max_dx = 0.0_dp
        do k = 1, NSTEP
            shift = DRIFT * k
            nfev = 0
            call multiroot_hybrids(system_fn, N, xc, xc_next, sc, ftol=FTOL, &
                ctx=shift)
            cold_total = cold_total + nfev
            nfev = 0
            call multiroot_hybrids(system_fn, N, xw, xw_next, sw, ftol=FTOL, &
                ctx=shift, warm=warm)
            warm_total = warm_total + nfev
            if (noise == 0.0_dp .and. &
                (sc%code /= FORTNUM_OK .or. sw%code /= FORTNUM_OK)) then
                write (error_unit, "(a,i0,1x,a,1x,a)") "FAIL status step ", k, &
                    trim(sc%msg), trim(sw%msg)
                nfail = nfail + 1
                return
            end if
            call system_fn(xw_next, f, shift)
            if (maxval(abs(f)) > max(FTOL, 10.0_dp * noise)) then
                write (error_unit, "(a,i0,es12.4)") "FAIL warm residual step ", &
                    k, maxval(abs(f))
                nfail = nfail + 1
            end if
            if (.not. warm%valid) then
                write (error_unit, "(a,i0)") "FAIL warm not valid after step ", k
                nfail = nfail + 1
            end if
            max_dx = max(max_dx, maxval(abs(xw_next - xc_next)))
            xc = xc_next
            xw = xw_next
        end do
        write (*, "(a,es8.1,a,i0,a,i0,a,es10.3)") "sequence noise=", noise, &
            ": cold evaluations=", cold_total, " warm=", warm_total, " max |x_warm - x_cold|=", max_dx
        ! Residuals within FTOL (or the noise) and a well-conditioned Jacobian
        ! bound the root difference by a small multiple of that.
        if (max_dx > max(1.0e-10_dp, 100.0_dp * noise)) then
            write (error_unit, "(a,es12.4)") "FAIL warm and cold roots differ ", &
                max_dx
            nfail = nfail + 1
        end if
        ! A warm solve skips the N evaluations of the difference Jacobian and
        ! pays a few extra secant iterations; cold costs 1 + N + ~3 here.  At
        ! the noise floor, secant updates from noise must not reach the next
        ! solve, or the saving shrinks.
        if (10 * warm_total > max_ratio * cold_total) then
            write (error_unit, "(a,i0,a,i0)") "FAIL warm evaluations ", &
                warm_total, " vs cold ", cold_total
            nfail = nfail + 1
        end if
    end subroutine check_sequence

    subroutine check_garbage()
        real(dp) :: x(N), f(N), shift
        type(fortnum_status_t) :: status
        type(multiroot_hybrids_warm_t) :: warm
        integer :: i, j

        ! Orthogonal q and an ill-scaled r of the wrong sign: a useless model.
        allocate (warm%q(N, N), warm%r(N, N))
        warm%q = 0.0_dp
        warm%r = 0.0_dp
        do i = 1, N
            warm%q(N + 1 - i, i) = 1.0_dp
            do j = i, N
                warm%r(i, j) = -1.0e3_dp / (i + j)
            end do
        end do
        warm%valid = .true.
        shift = 0.05_dp
        nfev = 0
        call multiroot_hybrids(system_fn, N, spread(-1.0_dp, 1, N), x, status, &
            ftol=FTOL, ctx=shift, warm=warm)
        call system_fn(x, f, shift)
        write (*, "(a,i0,1x,a)") "garbage warm: evaluations=", nfev, trim(status%msg)
        if (status%code /= FORTNUM_OK .or. maxval(abs(f)) > FTOL) then
            write (error_unit, "(a,es12.4)") "FAIL garbage warm residual ", &
                maxval(abs(f))
            nfail = nfail + 1
        end if
    end subroutine check_garbage

    subroutine check_size_mismatch()
        real(dp) :: x(N), xref(N), shift
        type(fortnum_status_t) :: status, sref
        type(multiroot_hybrids_warm_t) :: warm
        integer :: nfev_ref

        allocate (warm%q(N - 1, N - 1), warm%r(N - 1, N - 1))
        warm%q = 0.0_dp
        warm%r = 0.0_dp
        warm%valid = .true.
        shift = 0.05_dp
        nfev = 0
        call multiroot_hybrids(system_fn, N, spread(-1.0_dp, 1, N), xref, sref, &
            ftol=FTOL, ctx=shift)
        nfev_ref = nfev
        nfev = 0
        call multiroot_hybrids(system_fn, N, spread(-1.0_dp, 1, N), x, status, &
            ftol=FTOL, ctx=shift, warm=warm)
        ! Ignoring it means the cold path, evaluation for evaluation.
        if (status%code /= FORTNUM_OK .or. nfev /= nfev_ref .or. &
            any(x /= xref) .or. size(warm%q, 1) /= N) then
            write (error_unit, "(a,i0,a,i0)") "FAIL size mismatch: evaluations ", &
                nfev, " vs ", nfev_ref
            nfail = nfail + 1
        end if
    end subroutine check_size_mismatch

    ! F(x) = x^2 + 1 stagnates at x = 0 with a finite estimate, which is
    ! kept; a NaN residual yields non-finite factors, which are dropped.
    subroutine check_failure_paths()
        real(dp) :: x(1)
        type(fortnum_status_t) :: status
        type(multiroot_hybrids_warm_t) :: warm

        call multiroot_hybrids(system_fn, 1, [1.0_dp], x, status, ctx=NO_ROOT, &
            warm=warm)
        if (status%code == FORTNUM_OK .or. .not. warm%valid) then
            write (error_unit, "(a)") "FAIL no root: finite estimate dropped"
            nfail = nfail + 1
        end if
        call multiroot_hybrids(system_fn, 1, [1.0_dp], x, status, &
            ctx=NAN_RESIDUAL, warm=warm)
        if (status%code == FORTNUM_OK .or. warm%valid) then
            write (error_unit, "(a)") "FAIL NaN residual: estimate kept"
            nfail = nfail + 1
        end if
    end subroutine check_failure_paths

    ! ctx: real shift of the constant term (Broyden tridiagonal), or the
    ! integer NO_ROOT for F(x) = x^2 + 1 or NAN_RESIDUAL for F = NaN.
    subroutine system_fn(x, f, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:)
        class(*), intent(in), optional :: ctx
        integer :: i, m

        nfev = nfev + 1
        m = size(x)
        select type (ctx)
            type is (integer)
            f(1) = x(1)**2 + 1.0_dp
            if (ctx == NAN_RESIDUAL) f(1) = ieee_value(1.0_dp, ieee_quiet_nan)
            type is (real(dp))
            do i = 1, m
                f(i) = (3.0_dp - 2.0_dp*x(i))*x(i) + 1.0_dp + ctx
                if (i > 1) f(i) = f(i) - x(i - 1)
                if (i < m) f(i) = f(i) - 2.0_dp*x(i + 1)
                f(i) = f(i) + noise * bit_noise(x(i))
            end do
        end select
    end subroutine system_fn

    ! Deterministic value in [-1/2, 1/2) from the bits of x: unlike a smooth
    ! function it changes arbitrarily under the smallest step, like rounding.
    pure real(dp) function bit_noise(x)
        real(dp), intent(in) :: x
        integer(int64) :: h

        h = transfer(x, h)
        h = ieor(h, ishft(h, 13))
        h = ieor(h, ishft(h, -7))
        h = ieor(h, ishft(h, 17))
        bit_noise = real(iand(h, MASK52), dp) / 2.0_dp**52 - 0.5_dp
    end function bit_noise

end program test_fortnum_multiroot_hybrids_warm
