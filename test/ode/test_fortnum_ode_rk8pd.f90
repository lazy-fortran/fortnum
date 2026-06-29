program test_fortnum_ode_rk8pd
    ! Behavioural tests for the clean-room RK8(7)13M integrator
    ! (fortnum_ode_rk8pd).
    !   - scalar decay y' = -y reaches exp(-t) to near machine precision under
    !     the adaptive error-per-step control
    !   - one step matches the analytic eighth-order propagation on y' = -y
    !   - the adaptive control accepts a tight tolerance and grows the step on an
    !     easy stretch (the carried step ends larger than it started)
    !   - continuous re-entrant evolution across several output points lands on
    !     the same endpoint, to machine precision, as a single evolve over the
    !     whole interval (the schedule is carried, not reset)
    !   - backward integration lands on t1
    !   - bad input (uninitialised state, non-positive tolerance) is rejected

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_ode_rk8pd, only: rk8pd_state_t, rk8pd_step, &
        rk8pd_evolve_init, rk8pd_evolve_apply
    implicit none

    integer :: nfail

    nfail = 0
    call check_decay(nfail)
    call check_single_step_order(nfail)
    call check_control_grows(nfail)
    call check_reentrant_matches_single(nfail)
    call check_single_step_mesh(nfail)
    call check_backward(nfail)
    call check_bad_input(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "fortnum_ode_rk8pd: all tests passed"
    stop 0

contains

    subroutine rhs_decay(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
            dydt(1) = -y(1)
        end subroutine rhs_decay

        subroutine rhs_osc(t, y, dydt, ctx)
            real(dp), intent(in)  :: t
            real(dp), intent(in)  :: y(:)
            real(dp), intent(out) :: dydt(:)
            class(*), intent(in), optional :: ctx
            associate (unused_t => t); end associate
                dydt(1) =  y(2)
                dydt(2) = -y(1)
            end subroutine rhs_osc

            subroutine check_decay(nfail)
                integer, intent(inout) :: nfail
                type(rk8pd_state_t)    :: st
                type(fortnum_status_t) :: status
                real(dp) :: y(1), t, exact, errabs
                integer  :: nfev
                call rk8pd_evolve_init(st, 1, 1.0e-3_dp, status)
                y = [1.0_dp]
                t = 0.0_dp
                nfev = 0
                call rk8pd_evolve_apply(rhs_decay, st, t, 4.0_dp, y, 1.0e-13_dp, &
                    1.0e-13_dp, 100000, nfev, status)
                exact = exp(-4.0_dp)
                errabs = abs(y(1) - exact)
                if (status%code /= FORTNUM_OK .or. errabs > 1.0e-11_dp .or. &
                    abs(t - 4.0_dp) > 1.0e-13_dp) then
                    write (error_unit, "(a,es14.6,a,es14.6,a,i0)") &
                        "FAIL check_decay: y=", y(1), " err=", errabs, &
                        " code=", status%code
                    nfail = nfail + 1
                end if
            end subroutine check_decay

            ! One step of size h on y' = -y from y0 returns the eighth-order value; the
            ! local error vs exp(-h) y0 must be at the eighth-order level (~h^9).
            subroutine check_single_step_order(nfail)
                integer, intent(inout) :: nfail
                real(dp) :: y0(1), h
                real(dp) :: k1(1), k2(1), k3(1), k4(1), k5(1), k6(1)
                real(dp) :: k7(1), k8(1), k9(1), k10(1), k11(1), k12(1), k13(1)
                real(dp) :: ytmp(1), y8(1), yerr(1)
                real(dp) :: e1, e2, ratio
                integer  :: nfev
                ! Steps large enough that the eighth-order leading term dominates and
                ! both local errors stay above the rounding floor on y' = -y.
                y0 = [1.0_dp]
                h = 1.0_dp
                nfev = 0
                call rk8pd_step(rhs_decay, 0.0_dp, y0, h, .false., k1, k2, k3, k4, k5, &
                    k6, k7, k8, k9, k10, k11, k12, k13, ytmp, y8, yerr, nfev)
                e1 = abs(y8(1) - y0(1) * exp(-h))
                h = 0.5_dp
                nfev = 0
                call rk8pd_step(rhs_decay, 0.0_dp, y0, h, .false., k1, k2, k3, k4, k5, &
                    k6, k7, k8, k9, k10, k11, k12, k13, ytmp, y8, yerr, nfev)
                e2 = abs(y8(1) - y0(1) * exp(-h))
                ratio = e1 / e2
                write (*, "(a,es12.4,a,f10.2)") "rk8pd local err h=1.0:", e1, &
                    "  halving ratio=", ratio
                ! Eighth order: halving h shrinks the local error by ~2^9 = 512. Require
                ! at least 2^7 = 128 (one order of safety below theory).
                if (ratio < 128.0_dp) then
                    write (error_unit, "(a,f12.2)") &
                        "FAIL single_step_order: halving ratio too low ", ratio
                    nfail = nfail + 1
                end if
            end subroutine check_single_step_order

            ! On the easy decay stretch the adaptive control should accept and grow the
            ! carried step: starting from a tiny h0 the state step ends much larger.
            subroutine check_control_grows(nfail)
                integer, intent(inout) :: nfail
                type(rk8pd_state_t)    :: st
                type(fortnum_status_t) :: status
                real(dp) :: y(1), t, h_start
                integer  :: nfev
                h_start = 1.0e-4_dp
                call rk8pd_evolve_init(st, 1, h_start, status)
                y = [1.0_dp]
                t = 0.0_dp
                nfev = 0
                call rk8pd_evolve_apply(rhs_decay, st, t, 5.0_dp, y, 1.0e-9_dp, &
                    1.0e-9_dp, 100000, nfev, status)
                ! After integrating an easy problem with a loose tolerance the carried
                ! step must have grown well past its tiny start.
                if (status%code /= FORTNUM_OK .or. st%h <= 10.0_dp * h_start) then
                    write (error_unit, "(a,es14.6,a,i0)") &
                        "FAIL control_grows: final st%h=", st%h, " code=", status%code
                    nfail = nfail + 1
                end if
            end subroutine check_control_grows

            ! Re-entrant continuity: evolving across many intermediate output points,
            ! carrying state%h, must reproduce a single evolve over the whole interval
            ! to machine precision. Uses the harmonic oscillator so the schedule is
            ! non-trivial.
            subroutine check_reentrant_matches_single(nfail)
                integer, intent(inout) :: nfail
                type(rk8pd_state_t)    :: st_seg, st_one
                type(fortnum_status_t) :: status
                real(dp) :: y_seg(2), y_one(2), t, t1, h0, gap
                integer  :: nfev, i, nseg
                real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
                h0 = 1.0e-3_dp
                t1 = 6.0_dp * PI
                nseg = 24

                ! Single evolve over [0, t1].
                call rk8pd_evolve_init(st_one, 2, h0, status)
                y_one = [1.0_dp, 0.0_dp]
                t = 0.0_dp
                nfev = 0
                call rk8pd_evolve_apply(rhs_osc, st_one, t, t1, y_one, 1.0e-12_dp, &
                    1.0e-12_dp, 1000000, nfev, status)

                ! Segmented evolve over nseg equal sub-intervals, carrying the state.
                call rk8pd_evolve_init(st_seg, 2, h0, status)
                y_seg = [1.0_dp, 0.0_dp]
                t = 0.0_dp
                nfev = 0
                do i = 1, nseg
                    call rk8pd_evolve_apply(rhs_osc, st_seg, t, &
                        t1 * real(i, dp) / real(nseg, dp), y_seg, 1.0e-12_dp, &
                        1.0e-12_dp, 1000000, nfev, status)
                end do

                gap = max(abs(y_seg(1) - y_one(1)), abs(y_seg(2) - y_one(2)))
                write (*, "(a,es12.4)") "rk8pd re-entrant vs single endpoint gap=", gap
                ! Carrying the schedule should give a near-identical trajectory: the clip
                ! onto each intermediate point is local and the step continues unchanged.
                ! Allow only a little drift from the extra clip arithmetic.
                if (status%code /= FORTNUM_OK .or. gap > 1.0e-12_dp) then
                    write (error_unit, "(a,es14.6,a,i0)") &
                        "FAIL reentrant_matches_single: gap=", gap, " code=", status%code
                    nfail = nfail + 1
                end if
            end subroutine check_reentrant_matches_single

            ! single_step mode returns after one accepted step, mirroring a single
            ! adaptive evolve advance. Looping until t == t1 and recording each
            ! step must (a) advance t strictly toward t1 each call, (b) land exactly on
            ! t1, and (c) reach the same endpoint as one full continuous evolve. This is
            ! the contract KiLCA's imhd zone solver relies on to rebuild the evolve mesh.
            subroutine check_single_step_mesh(nfail)
                integer, intent(inout) :: nfail
                type(rk8pd_state_t)    :: st_step, st_one
                type(fortnum_status_t) :: status
                real(dp) :: y_step(2), y_one(2), t, t_prev, t1, h0, gap
                integer  :: nfev, nsteps
                real(dp), parameter :: PI = 3.14159265358979323846264338327950288_dp
                h0 = 1.0e-3_dp
                t1 = 6.0_dp * PI

                ! One full continuous evolve over [0, t1].
                call rk8pd_evolve_init(st_one, 2, h0, status)
                y_one = [1.0_dp, 0.0_dp]
                t = 0.0_dp
                nfev = 0
                call rk8pd_evolve_apply(rhs_osc, st_one, t, t1, y_one, 1.0e-12_dp, &
                    1.0e-12_dp, 1000000, nfev, status)

                ! Step-by-step evolve recording each accepted step.
                call rk8pd_evolve_init(st_step, 2, h0, status)
                y_step = [1.0_dp, 0.0_dp]
                t = 0.0_dp
                nfev = 0
                nsteps = 0
                do
                    if (t >= t1) exit
                    t_prev = t
                    call rk8pd_evolve_apply(rhs_osc, st_step, t, t1, y_step, &
                        1.0e-12_dp, 1.0e-12_dp, 1000000, nfev, status, &
                        single_step=.true.)
                    if (status%code /= FORTNUM_OK) exit
                    if (.not. (t > t_prev)) then
                        write (error_unit, "(a,es14.6)") &
                            "FAIL single_step_mesh: t did not advance, t=", t
                        nfail = nfail + 1
                        return
                    end if
                    nsteps = nsteps + 1
                    if (nsteps > 1000000) exit
                end do

                gap = max(abs(y_step(1) - y_one(1)), abs(y_step(2) - y_one(2)))
                write (*, "(a,i0,a,es12.4)") &
                    "rk8pd single-step mesh: nsteps=", nsteps, "  endpoint gap=", gap
                if (status%code /= FORTNUM_OK .or. abs(t - t1) > 1.0e-13_dp .or. &
                    gap > 1.0e-12_dp .or. nsteps < 2) then
                    write (error_unit, "(a,es14.6,a,es14.6,a,i0)") &
                        "FAIL single_step_mesh: t=", t, " gap=", gap, &
                        " code=", status%code
                    nfail = nfail + 1
                end if
            end subroutine check_single_step_mesh

            subroutine check_backward(nfail)
                integer, intent(inout) :: nfail
                type(rk8pd_state_t)    :: st
                type(fortnum_status_t) :: status
                real(dp) :: y(1), t, exact
                integer  :: nfev
                call rk8pd_evolve_init(st, 1, 1.0e-3_dp, status)
                y = [exp(-2.0_dp)]
                t = 2.0_dp
                nfev = 0
                call rk8pd_evolve_apply(rhs_decay, st, t, 0.0_dp, y, 1.0e-12_dp, &
                    1.0e-12_dp, 100000, nfev, status)
                exact = 1.0_dp
                if (status%code /= FORTNUM_OK .or. abs(y(1) - exact) > 1.0e-10_dp .or. &
                    abs(t) > 1.0e-13_dp) then
                    write (error_unit, "(a,es14.6,a,i0)") &
                        "FAIL check_backward: y=", y(1), " code=", status%code
                    nfail = nfail + 1
                end if
            end subroutine check_backward

            subroutine check_bad_input(nfail)
                integer, intent(inout) :: nfail
                type(rk8pd_state_t)    :: st
                type(fortnum_status_t) :: status
                real(dp) :: y(1), t
                integer  :: nfev
                ! evolve_init rejects a non-positive starting step.
                call rk8pd_evolve_init(st, 1, 0.0_dp, status)
                if (status%code /= FORTNUM_DOMAIN_ERROR) then
                    write (error_unit, "(a,i0)") &
                        "FAIL bad_input (h0=0): code=", status%code
                    nfail = nfail + 1
                end if
                ! evolve_apply on an uninitialised state for this neq is rejected.
                y = [1.0_dp]
                t = 0.0_dp
                nfev = 0
                call rk8pd_evolve_apply(rhs_decay, st, t, 1.0_dp, y, 1.0e-12_dp, &
                    1.0e-12_dp, 100000, nfev, status)
                if (status%code /= FORTNUM_DOMAIN_ERROR) then
                    write (error_unit, "(a,i0)") &
                        "FAIL bad_input (uninit state): code=", status%code
                    nfail = nfail + 1
                end if
                ! Initialised, but a non-positive eps_abs is rejected.
                call rk8pd_evolve_init(st, 1, 1.0e-3_dp, status)
                y = [1.0_dp]
                t = 0.0_dp
                nfev = 0
                call rk8pd_evolve_apply(rhs_decay, st, t, 1.0_dp, y, 0.0_dp, &
                    1.0e-12_dp, 100000, nfev, status)
                if (status%code /= FORTNUM_DOMAIN_ERROR) then
                    write (error_unit, "(a,i0)") &
                        "FAIL bad_input (eps_abs=0): code=", status%code
                    nfail = nfail + 1
                end if
            end subroutine check_bad_input

        end program test_fortnum_ode_rk8pd
