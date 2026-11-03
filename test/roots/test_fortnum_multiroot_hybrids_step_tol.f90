program test_fortnum_multiroot_hybrids_step_tol
    ! Step-size convergence of multiroot_hybrids (optional step_tol).
    !
    ! Broyden tridiagonal system (More, Garbow & Hillstrom, ACM TOMS 7 (1981)
    ! 17, problem 30) with a deterministic residual noise of amplitude 1e-9,
    ! far above ftol = 1e-12, as for a caller whose residual has a rounding
    ! floor above any useful ftol.  Oracle: the root of the noiseless system
    ! solved to |F|_inf <= 1e-14.  Its Jacobian is well conditioned near the
    ! root, so the noise moves the root by about 1e-9, below the step
    ! tolerance.  Without step_tol every noisy solve ends in a failure status;
    ! with it the solve must succeed within step_tol of the oracle root and
    ! need fewer evaluations.  Without noise it must agree with the ftol
    ! solve to step_tol, and a step_tol of zero must change nothing.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use fortnum_multiroot, only: multiroot_hybrids, multiroot_hybrids_warm_t
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    integer, parameter :: N = 8
    real(dp), parameter :: FTOL = 1.0e-12_dp, STEP = 1.0e-8_dp
    integer(int64), parameter :: MASK52 = 2_int64**52 - 1
    integer :: nfev, nfail
    real(dp) :: noise, xref(N), x(N), xftol(N), x0(N)
    type(fortnum_status_t) :: status
    integer :: nfev_ftol, nfev_step

    nfail = 0
    x0 = -1.0_dp

    noise = 0.0_dp
    call multiroot_hybrids(system_fn, N, x0, xref, status, ftol=1.0e-14_dp)
    call expect(status%code == FORTNUM_OK, "noiseless reference solve")

    noise = 1.0e-9_dp
    nfev = 0
    call multiroot_hybrids(system_fn, N, x0, xftol, status, ftol=FTOL)
    nfev_ftol = nfev
    call expect(status%code /= FORTNUM_OK, &
        "noisy solve without step_tol fails at the noise floor")

    nfev = 0
    call multiroot_hybrids(system_fn, N, x0, x, status, ftol=FTOL, &
        step_tol=spread(STEP, 1, N))
    nfev_step = nfev
    write (*, "(a,i0,a,i0,a,es10.3)") "noisy: evaluations ftol-only=", &
        nfev_ftol, " step_tol=", nfev_step, " |x - xref|=", &
        maxval(abs(x - xref))
    call expect(status%code == FORTNUM_OK, "noisy solve with step_tol")
    call expect(maxval(abs(x - xref)) <= STEP, "noisy root within step_tol")
    call expect(nfev_step < nfev_ftol, "step_tol stops before stagnation")

    call check_warm_sequence()

    noise = 0.0_dp
    nfev = 0
    call multiroot_hybrids(system_fn, N, x0, xftol, status, ftol=FTOL)
    nfev_ftol = nfev
    nfev = 0
    call multiroot_hybrids(system_fn, N, x0, x, status, ftol=FTOL, &
        step_tol=spread(STEP, 1, N))
    nfev_step = nfev
    call expect(status%code == FORTNUM_OK, "noiseless solve with step_tol")
    call expect(maxval(abs(x - xref)) <= STEP, "noiseless root within step_tol")
    call expect(nfev_step <= nfev_ftol, "noiseless step_tol costs no more")

    nfev = 0
    call multiroot_hybrids(system_fn, N, x0, x, status, ftol=FTOL, &
        step_tol=spread(0.0_dp, 1, N))
    call expect(status%code == FORTNUM_OK .and. nfev == nfev_ftol .and. &
        all(x == xftol), "zero step_tol leaves the ftol solve unchanged")

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " check(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    ! A drifting sequence of noisy solves with a warm start, as an implicit
    ! time stepper runs it: each must succeed, and the warm estimate saved on
    ! a step_tol exit must stay usable for the next solve.
    subroutine check_warm_sequence()
        type(multiroot_hybrids_warm_t) :: warm
        real(dp) :: xa(N), xb(N), shift_ref(N)
        integer :: k, total

        noise = 1.0e-9_dp
        xa = x0
        total = 0
        do k = 1, 20
            nfev = 0
            call multiroot_hybrids(system_fn, N, xa, xb, status, ftol=FTOL, &
                ctx=2.0e-3_dp*k, warm=warm, step_tol=spread(STEP, 1, N))
            total = total + nfev
            call expect(status%code == FORTNUM_OK, "warm noisy sequence")
            call expect(warm%valid, "warm estimate kept on step_tol exit")
            noise = 0.0_dp
            call multiroot_hybrids(system_fn, N, xb, shift_ref, status, &
                ftol=1.0e-14_dp, ctx=2.0e-3_dp*k)
            noise = 1.0e-9_dp
            call expect(maxval(abs(xb - shift_ref)) <= STEP, &
                "warm noisy root within step_tol")
            xa = xb
        end do
        write (*, "(a,i0)") "warm noisy sequence: evaluations=", total
    end subroutine check_warm_sequence

    subroutine expect(ok, what)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: what
        if (.not. ok) then
            write (error_unit, "(a,a)") "FAIL ", what
            nfail = nfail + 1
        end if
    end subroutine expect

    subroutine system_fn(x, f, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: shift
        integer :: i, m

        nfev = nfev + 1
        shift = 0.0_dp
        if (present(ctx)) then
            select type (ctx)
                type is (real(dp))
                shift = ctx
            end select
        end if
        m = size(x)
        do i = 1, m
            f(i) = (3.0_dp - 2.0_dp*x(i))*x(i) + 1.0_dp + shift
            if (i > 1) f(i) = f(i) - x(i - 1)
            if (i < m) f(i) = f(i) - 2.0_dp*x(i + 1)
            f(i) = f(i) + noise * bit_noise(x(i))
        end do
    end subroutine system_fn

    pure real(dp) function bit_noise(x)
        real(dp), intent(in) :: x
        integer(int64) :: h

        h = transfer(x, h)
        h = ieor(h, ishft(h, 13))
        h = ieor(h, ishft(h, -7))
        h = ieor(h, ishft(h, 17))
        bit_noise = real(iand(h, MASK52), dp) / 2.0_dp**52 - 0.5_dp
    end function bit_noise

end program test_fortnum_multiroot_hybrids_step_tol
