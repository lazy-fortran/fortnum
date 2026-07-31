program test_fortnum_ode_extrapolation
    ! Behavioural tests for the Gragg-Bulirsch-Stoer extrapolation integrator.
    !
    ! The claim that matters is the order claim: column j of the extrapolation
    ! table has order 2j. That is checked directly, by running the modified
    ! midpoint rule at a fixed column count and halving the step, on a NONLINEAR
    ! problem -- a linear one can mask order defects because the extrapolation
    ! becomes exact for the wrong reason.
    !
    ! Oracles: closed-form solutions, an exact invariant (Kepler energy), and
    ! comparison against DOP853, an unrelated method. No recorded output.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_ode, only: ode_problem_t, ode_solution_t
    use fortnum_ode_extrapolation, only: ode_integrate_gbs, ode_solve_gbs, &
                                         gbs_modified_midpoint, gbs_step_sequence
    use fortnum_ode_dop853, only: ode_solve_dop
    implicit none

    integer :: nfail

    nfail = 0
    call check_midpoint_order(nfail)
    call check_extrapolated_order(nfail)
    call check_nonlinear_scalar(nfail)
    call check_kepler_energy(nfail)
    call check_vs_dop853(nfail)
    call check_backward(nfail)
    call check_bad_input(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "fortnum_ode_extrapolation: all tests passed"

contains

    ! y' = -y^2, y(0) = 1 => y(t) = 1/(1+t)
    subroutine rhs_nl(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) = -y(1)*y(1)
    end subroutine rhs_nl

    subroutine rhs_kepler(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: r3
        associate (unused_t => t); end associate
        r3 = (sqrt(y(1)**2 + y(2)**2))**3
        dydt(1) = y(3)
        dydt(2) = y(4)
        dydt(3) = -y(1)/r3
        dydt(4) = -y(2)/r3
    end subroutine rhs_kepler

    real(dp) function kepler_energy(y) result(e)
        real(dp), intent(in) :: y(:)
        e = 0.5_dp*(y(3)**2 + y(4)**2) - 1.0_dp/sqrt(y(1)**2 + y(2)**2)
    end function kepler_energy

    real(dp) function nl_exact(t) result(v)
        real(dp), intent(in) :: t
        v = 1.0_dp/(1.0_dp + t)
    end function nl_exact

    ! One un-extrapolated modified midpoint sweep is second order.
    subroutine check_midpoint_order(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: y(1), zprev(1), zcur(1), ftmp(1), yout(1)
        real(dp) :: err(3), ordobs, tend
        integer  :: nfev, i, nsub(3)

        tend = 1.0_dp
        nsub = [8, 16, 32]
        do i = 1, 3
            y(1) = 1.0_dp
            nfev = 0
            call gbs_modified_midpoint(rhs_nl, 0.0_dp, y, tend, nsub(i), zprev, &
                                       zcur, ftmp, yout, nfev)
            err(i) = abs(yout(1) - nl_exact(tend))
        end do
        ordobs = log(err(2)/err(3))/log(2.0_dp)
        write (*, "(a,f6.2,a,es10.3)") "  modified midpoint observed order ", &
            ordobs, "   err ", err(3)
        if (abs(ordobs - 2.0_dp) > 0.3_dp) then
            write (error_unit, "(a,f6.2)") "  midpoint order not 2: ", ordobs
            nfail = nfail + 1
        end if
    end subroutine check_midpoint_order

    ! Extrapolating j columns must give order 2j. Verified for j = 2 and j = 3
    ! by building the Neville table by hand from midpoint sweeps, which keeps
    ! the check independent of the integrator's own step control.
    subroutine check_extrapolated_order(nfail)
        integer, intent(inout) :: nfail

        integer, parameter :: NCOL = 3
        real(dp) :: y(1), zprev(1), zcur(1), ftmp(1), yout(1)
        real(dp) :: tab(NCOL), aux, cnew, ratio, denom
        real(dp) :: e2(2), e3(2), ord2, ord3, tend, base
        integer  :: nseq(NCOL), nfev, i, j, k, ih
        real(dp) :: hscale(2)

        call gbs_step_sequence(nseq)
        tend = 1.0_dp
        hscale = [1.0_dp, 0.5_dp]

        do ih = 1, 2
            ! Whole interval covered in 1/hscale steps of the extrapolation, but
            ! for an order check a single step of size h suffices.
            base = tend*hscale(ih)
            tab = 0.0_dp
            do k = 1, NCOL
                y(1) = 1.0_dp
                nfev = 0
                call gbs_modified_midpoint(rhs_nl, 0.0_dp, y, base, nseq(k), &
                                           zprev, zcur, ftmp, yout, nfev)
                aux = tab(1)
                tab(1) = yout(1)
                do j = 2, k
                    ratio = real(nseq(k), dp)/real(nseq(k - j + 1), dp)
                    denom = ratio*ratio - 1.0_dp
                    cnew = tab(j - 1) + (tab(j - 1) - aux)/denom
                    aux = tab(j)
                    tab(j) = cnew
                end do
                if (k == 2) e2(ih) = abs(tab(2) - nl_exact(base))
                if (k == 3) e3(ih) = abs(tab(3) - nl_exact(base))
            end do
        end do

        ord2 = log(e2(1)/e2(2))/log(2.0_dp)
        ord3 = log(e3(1)/e3(2))/log(2.0_dp)
        write (*, "(a,f6.2,a,f6.2)") "  extrapolated order: 2 columns ", ord2, &
            "   3 columns ", ord3

        ! Column j gains two orders over the previous one. Measured this way --
        ! one step of size H, with the substep count held fixed per column, so
        ! the substep scales with H -- the error of column j goes as H^(2j).
        ! The un-extrapolated sweep above is the j = 1 case and gives 2, so the
        ! whole ladder is consistent: 2, 4, 6.
        if (abs(ord2 - 4.0_dp) > 0.6_dp) then
            write (error_unit, "(a,f6.2,a)") "  2-column order ", ord2, &
                ", expected about 4"
            nfail = nfail + 1
        end if
        if (abs(ord3 - 6.0_dp) > 0.8_dp) then
            write (error_unit, "(a,f6.2,a)") "  3-column order ", ord3, &
                ", expected about 6"
            nfail = nfail + 1
        end if
    end subroutine check_extrapolated_order

    subroutine check_nonlinear_scalar(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: status
        real(dp) :: y0(1), err
        real(dp), allocatable :: t_out(:), y_out(:, :)

        y0(1) = 1.0_dp
        call ode_solve_gbs(rhs_nl, 0.0_dp, 4.0_dp, y0, t_out, y_out, status, &
                           rtol=1.0e-12_dp)
        err = abs(y_out(1, size(t_out)) - nl_exact(4.0_dp))
        write (*, "(a,es10.3)") "  nonlinear scalar endpoint error ", err
        if (status%code /= FORTNUM_OK .or. .not. (err < 1.0e-11_dp)) then
            write (error_unit, "(a)") "  nonlinear scalar: endpoint inaccurate"
            nfail = nfail + 1
        end if
    end subroutine check_nonlinear_scalar

    subroutine check_kepler_energy(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp) :: y0(4), e0, e1, drift
        real(dp), parameter :: ecc = 0.5_dp
        real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)
        integer :: n

        y0(1) = 1.0_dp - ecc
        y0(2) = 0.0_dp
        y0(3) = 0.0_dp
        y0(4) = sqrt((1.0_dp + ecc)/(1.0_dp - ecc))
        e0 = kepler_energy(y0)

        problem%rhs => rhs_kepler
        problem%t0 = 0.0_dp
        problem%t1 = 50.0_dp*twopi
        problem%y0 = y0
        problem%rtol = 1.0e-12_dp
        problem%max_steps = 200000

        call ode_integrate_gbs(problem, solution, status)
        n = size(solution%t)
        e1 = kepler_energy(solution%y(:, n))
        drift = abs((e1 - e0)/e0)

        write (*, "(a,es10.3,a,i0,a,i0)") "  kepler 50 orbits: rel energy drift ", &
            drift, "   steps ", solution%nsteps, "  nfev ", solution%nfev
        ! GBS is not symplectic and carries no special structure, so its energy
        ! does drift -- around 1e-8 over 50 orbits here, against 1e-14 over 200
        ! orbits for the Gauss-Radau integrator. That gap is a benchmark result,
        ! not a defect; the assertion only catches catastrophic drift.
        if (status%code /= FORTNUM_OK .or. .not. (drift < 1.0e-6_dp)) then
            write (error_unit, "(a)") "  kepler: energy drift is catastrophic"
            nfail = nfail + 1
        end if
    end subroutine check_kepler_energy

    subroutine check_vs_dop853(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s1, s2
        real(dp) :: y0(1), eg, ed
        real(dp), allocatable :: tg(:), yg(:, :), td(:), yd(:, :)

        y0(1) = 1.0_dp
        call ode_solve_gbs(rhs_nl, 0.0_dp, 4.0_dp, y0, tg, yg, s1, rtol=1.0e-10_dp)
        call ode_solve_dop(rhs_nl, 0.0_dp, 4.0_dp, y0, td, yd, s2, rtol=1.0e-10_dp)
        eg = abs(yg(1, size(tg)) - 0.2_dp)
        ed = abs(yd(1, size(td)) - 0.2_dp)
        write (*, "(a,es10.3,a,es10.3)") "  vs dop853 at rtol 1e-10: gbs ", eg, &
            "   dop853 ", ed
        if (.not. (eg <= max(ed*50.0_dp, 1.0e-11_dp))) then
            write (error_unit, "(a)") "  gbs markedly less accurate than dop853"
            nfail = nfail + 1
        end if
    end subroutine check_vs_dop853

    subroutine check_backward(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: status
        real(dp) :: y0(1), yend, tend
        real(dp), allocatable :: t_out(:), y_out(:, :)

        y0(1) = nl_exact(2.0_dp)
        call ode_solve_gbs(rhs_nl, 2.0_dp, 0.0_dp, y0, t_out, y_out, status, &
                           rtol=1.0e-12_dp)
        tend = t_out(size(t_out))
        yend = y_out(1, size(t_out))
        write (*, "(a,es10.3)") "  backward endpoint error ", abs(yend - 1.0_dp)
        if (status%code /= FORTNUM_OK .or. abs(tend) > 1.0e-13_dp .or. &
            .not. (abs(yend - 1.0_dp) < 1.0e-9_dp)) then
            write (error_unit, "(a)") "  backward integration failed"
            nfail = nfail + 1
        end if
    end subroutine check_backward

    subroutine check_bad_input(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp) :: y0(1)

        y0(1) = 1.0_dp
        problem%rhs => rhs_nl
        problem%t0 = 0.0_dp
        problem%t1 = 0.0_dp
        problem%y0 = y0
        call ode_integrate_gbs(problem, solution, status)
        if (status%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, "(a)") "  zero-span problem was not rejected"
            nfail = nfail + 1
        end if

        problem%t1 = 1.0_dp
        problem%rtol = -1.0_dp
        call ode_integrate_gbs(problem, solution, status)
        if (status%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, "(a)") "  negative rtol was not rejected"
            nfail = nfail + 1
        end if
    end subroutine check_bad_input

end program test_fortnum_ode_extrapolation
