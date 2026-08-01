program test_fortnum_ode_tdrk
    ! Tests the claim the TDRK module rests on: an RKN (Nystrom) tableau keeps
    ! its order when applied to a FIRST-order system via the reformulation
    ! zddot = F'(z) F(z), with the velocity re-synchronised each step.
    !
    ! The claim is provable -- an order-p RKN method on zddot = G(z) with exact
    ! initial velocity has O(h^(p+1)) local error by definition, and the re-sync
    ! supplies exactly that -- so what is really under test is the
    ! implementation and the tableau construction.
    !
    ! Oracles, none of them a recording of this code:
    !   - the Nystrom position order conditions, checked arithmetically, and
    !     checked to be SHARP (the next order must fail, or the comparison
    !     would be vacuous)
    !   - closed-form solutions for the scalar problem
    !   - a high-resolution classical RK4 run for Kepler and Lorentz, which is a
    !     different method
    !   - exact conservation of |v|^2 in a static magnetic field, since a purely
    !     magnetic force does no work
    !
    ! Order is measured on NONLINEAR problems only. Linear problems can mask
    ! order-condition defects.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_ode_tdrk, only: tdrk_tableau, tdrk_integrate_fixed, &
                                rkng_integrate_fixed, tdrk_integrate_adaptive, &
                                rkng_integrate_adaptive, TDRK_MAX_STAGES
    implicit none

    integer :: nfail

    nfail = 0
    call check_order_conditions(nfail)
    call check_transplanted_order(nfail)
    call check_lorentz_energy(nfail)
    call check_rkng_lorentz(nfail)
    call check_tdrk_adaptive(nfail)
    call check_rkng_adaptive(nfail)
    call report_cost_model()

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "fortnum_ode_tdrk: all tests passed"

contains

    ! ------------------------------------------------------------- problems
    ! prob 1: nonlinear scalar  zdot = -z^2,  z(0)=1,  z(t)=1/(1+t)
    ! prob 2: Kepler as a first-order 4D system, eccentricity 0.5
    ! prob 3: Lorentz 6D in a static divergence-free non-uniform B

    subroutine f_scalar(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (u => t); end associate
        dydt(1) = -y(1)*y(1)
    end subroutine f_scalar

    subroutine g_scalar(t, y, g, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: g(:)
        class(*), intent(in), optional :: ctx
        associate (u => t); end associate
        ! F = -z^2, F' = -2z, so G = F'F = 2 z^3
        g(1) = 2.0_dp*y(1)**3
    end subroutine g_scalar

    subroutine f_kepler(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: r3
        associate (u => t); end associate
        r3 = (sqrt(y(1)**2 + y(2)**2))**3
        dydt(1) = y(3)
        dydt(2) = y(4)
        dydt(3) = -y(1)/r3
        dydt(4) = -y(2)/r3
    end subroutine f_kepler

    subroutine g_kepler(t, y, g, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: g(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: r, r3, r5, rv
        associate (u => t); end associate
        r = sqrt(y(1)**2 + y(2)**2)
        r3 = r**3
        r5 = r**5
        rv = y(1)*y(3) + y(2)*y(4)
        g(1) = -y(1)/r3
        g(2) = -y(2)/r3
        g(3) = -y(3)/r3 + 3.0_dp*y(1)*rv/r5
        g(4) = -y(4)/r3 + 3.0_dp*y(2)*rv/r5
    end subroutine g_kepler

    ! B(x) = B0 * (-x/(2L), -y/(2L), 1 + z/L): divergence free and non-uniform.
    subroutine lorentz_b(x, b)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: b(3)
        real(dp), parameter :: b0 = 1.7_dp, l = 2.5_dp
        b(1) = b0*(-x(1)/(2.0_dp*l))
        b(2) = b0*(-x(2)/(2.0_dp*l))
        b(3) = b0*(1.0_dp + x(3)/l)
    end subroutine lorentz_b

    subroutine lorentz_vgradb(v, db)
        real(dp), intent(in) :: v(3)
        real(dp), intent(out) :: db(3)
        real(dp), parameter :: b0 = 1.7_dp, l = 2.5_dp
        db(1) = b0*(-v(1)/(2.0_dp*l))
        db(2) = b0*(-v(2)/(2.0_dp*l))
        db(3) = b0*(v(3)/l)
    end subroutine lorentz_vgradb

    function cross(a, b) result(c)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: c(3)
        c(1) = a(2)*b(3) - a(3)*b(2)
        c(2) = a(3)*b(1) - a(1)*b(3)
        c(3) = a(1)*b(2) - a(2)*b(1)
    end function cross

    subroutine f_lorentz(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: b(3)
        associate (u => t); end associate
        call lorentz_b(y(1:3), b)
        dydt(1:3) = y(4:6)
        dydt(4:6) = cross(y(4:6), b)
    end subroutine f_lorentz

    subroutine g_lorentz(t, y, g, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: g(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: b(3), a(3), db(3)
        associate (u => t); end associate
        call lorentz_b(y(1:3), b)
        a = cross(y(4:6), b)
        call lorentz_vgradb(y(4:6), db)
        g(1:3) = a
        g(4:6) = cross(a, b) + cross(y(4:6), db)
    end subroutine g_lorentz

    ! Second-order form of the same Lorentz problem: y'' = f(t, y, y').
    subroutine f2_lorentz(t, y, yp, ypp, ctx)
        real(dp), intent(in) :: t, y(:), yp(:)
        real(dp), intent(out) :: ypp(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: b(3)
        associate (u => t); end associate
        call lorentz_b(y(1:3), b)
        ypp(1:3) = cross(yp(1:3), b)
    end subroutine f2_lorentz

    ! ----------------------------------------------------------- RK4 oracle
    subroutine rk4(prob, neq, y0, t0, t1, nsteps, yend)
        integer, intent(in) :: prob, neq, nsteps
        real(dp), intent(in) :: y0(:), t0, t1
        real(dp), intent(out) :: yend(:)

        real(dp) :: y(neq), k1(neq), k2(neq), k3(neq), k4(neq), h, t
        integer :: k

        h = (t1 - t0)/real(nsteps, dp)
        y = y0(1:neq)
        t = t0
        do k = 1, nsteps
            call f_any(prob, t, y, k1)
            call f_any(prob, t + 0.5_dp*h, y + 0.5_dp*h*k1, k2)
            call f_any(prob, t + 0.5_dp*h, y + 0.5_dp*h*k2, k3)
            call f_any(prob, t + h, y + h*k3, k4)
            y = y + (h/6.0_dp)*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)
            t = t0 + real(k, dp)*h
        end do
        yend(1:neq) = y
    end subroutine rk4

    subroutine f_any(prob, t, y, dydt)
        integer, intent(in) :: prob
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        select case (prob)
        case (1); call f_scalar(t, y, dydt)
        case (2); call f_kepler(t, y, dydt)
        case (3); call f_lorentz(t, y, dydt)
        end select
    end subroutine f_any

    ! ------------------------------------------------------------- checks

    ! The order claim must hold AND be sharp: if the next order's condition also
    ! held, matching the measured order would prove nothing.
    subroutine check_order_conditions(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: c(TDRK_MAX_STAGES), abar(TDRK_MAX_STAGES, TDRK_MAX_STAGES)
        real(dp) :: bbar(TDRK_MAX_STAGES)
        real(dp) :: s0, s1, s2, s3, stage
        integer :: s, i, icase, order
        logical :: ok, three
        real(dp), parameter :: tol = 1.0e-14_dp

        do icase = 1, 3
            select case (icase)
            case (1); order = 3; three = .false.
            case (2); order = 4; three = .false.
            case (3); order = 4; three = .true.
            end select

            call tdrk_tableau(order, three, s, c, abar, bbar, ok)
            if (.not. ok) then
                write (error_unit, "(a)") "  tableau construction failed"
                nfail = nfail + 1
                cycle
            end if

            do i = 1, s
                stage = sum(abar(i, 1:s))
                if (abs(stage - 0.5_dp*c(i)**2) > tol) then
                    write (error_unit, "(a,i0)") "  stage consistency fails at ", i
                    nfail = nfail + 1
                end if
            end do

            s0 = sum(bbar(1:s))
            s1 = sum(bbar(1:s)*c(1:s))
            s2 = sum(bbar(1:s)*c(1:s)**2)
            s3 = sum(bbar(1:s)*c(1:s)**3)

            if (abs(s0 - 0.5_dp) > tol .or. abs(s1 - 1.0_dp/6.0_dp) > tol) then
                write (error_unit, "(a,i0)") "  order-3 conditions fail, case ", icase
                nfail = nfail + 1
            end if

            if (order >= 4) then
                if (abs(s2 - 1.0_dp/12.0_dp) > tol) then
                    write (error_unit, "(a,i0)") "  order-4 condition fails, case ", &
                        icase
                    nfail = nfail + 1
                end if
                if (abs(s3 - 0.05_dp) <= tol) then
                    write (error_unit, "(a,i0)") "  order-4 claim not sharp, case ", &
                        icase
                    nfail = nfail + 1
                end if
            else
                if (abs(s2 - 1.0_dp/12.0_dp) <= tol) then
                    write (error_unit, "(a,i0)") "  order-3 claim not sharp, case ", &
                        icase
                    nfail = nfail + 1
                end if
            end if
        end do
        write (*, "(a)") "  tableau order conditions: hold and are sharp"
    end subroutine check_order_conditions

    subroutine check_transplanted_order(nfail)
        integer, intent(inout) :: nfail

        integer :: nsteps(4), prob, icase, neq, order, i, nff, nfg
        real(dp) :: y0(6), yref(6), yend(6), t1, err(4), ordobs, meanord
        type(fortnum_status_t) :: status
        logical :: three
        character(len=28) :: label
        character(len=30) :: pname

        nsteps = [200, 400, 800, 1600]

        do prob = 1, 3
            call problem_setup(prob, neq, y0, t1, pname)
            if (prob == 1) then
                yref = 0.0_dp
                yref(1) = 1.0_dp/(1.0_dp + t1)
            else
                call rk4(prob, neq, y0, 0.0_dp, t1, 2000000, yref)
            end if
            write (*, "(a,a)") "  problem: ", trim(pname)

            do icase = 1, 3
                select case (icase)
                case (1); order = 3; three = .false.; label = "RKN3 (2-stage)"
                case (2); order = 4; three = .false.; label = "RKN4 (2-stage)"
                case (3); order = 4; three = .true.; label = "RKN4 (3-stage coupled)"
                end select

                do i = 1, 4
                    call integrate_case(prob, y0, neq, t1, nsteps(i), order, three, &
                                        yend, nff, nfg, status)
                    err(i) = maxval(abs(yend(1:neq) - yref(1:neq)))
                end do

                meanord = 0.5_dp*(log(err(2)/err(3)) + log(err(3)/err(4)))/log(2.0_dp)
                write (*, "(a,a24,a,i0,a,f6.2,a,es9.2)") "    ", label, &
                    "  claimed ", order, "  observed ", meanord, "   err ", err(4)

                if (abs(meanord - real(order, dp)) > 0.3_dp) then
                    write (error_unit, "(a,f6.2,a,i0)") "    observed order ", &
                        meanord, " does not match ", order
                    nfail = nfail + 1
                end if
            end do
        end do
    end subroutine check_transplanted_order

    subroutine problem_setup(prob, neq, y0, t1, pname)
        integer, intent(in) :: prob
        integer, intent(out) :: neq
        real(dp), intent(out) :: y0(:), t1
        character(len=*), intent(out) :: pname

        real(dp), parameter :: ecc = 0.5_dp

        y0 = 0.0_dp
        select case (prob)
        case (1)
            neq = 1; t1 = 2.0_dp; y0(1) = 1.0_dp
            pname = "nonlinear scalar"
        case (2)
            neq = 4; t1 = 3.0_dp
            y0(1) = 1.0_dp - ecc
            y0(4) = sqrt((1.0_dp + ecc)/(1.0_dp - ecc))
            pname = "Kepler 4D (e=0.5)"
        case (3)
            neq = 6; t1 = 2.0_dp
            y0(1:6) = [0.3_dp, -0.2_dp, 0.1_dp, 0.4_dp, 0.9_dp, 0.35_dp]
            pname = "Lorentz 6D (non-uniform B)"
        end select
    end subroutine problem_setup

    subroutine integrate_case(prob, y0, neq, t1, nsteps, order, three, yend, &
                              nff, nfg, status)
        integer, intent(in) :: prob, neq, nsteps, order
        real(dp), intent(in) :: y0(:), t1
        logical, intent(in) :: three
        real(dp), intent(out) :: yend(:)
        integer, intent(out) :: nff, nfg
        type(fortnum_status_t), intent(out) :: status

        select case (prob)
        case (1)
            call tdrk_integrate_fixed(f_scalar, g_scalar, 0.0_dp, t1, y0(1:neq), &
                                      nsteps, order, three, yend, nff, nfg, status)
        case (2)
            call tdrk_integrate_fixed(f_kepler, g_kepler, 0.0_dp, t1, y0(1:neq), &
                                      nsteps, order, three, yend, nff, nfg, status)
        case (3)
            call tdrk_integrate_fixed(f_lorentz, g_lorentz, 0.0_dp, t1, y0(1:neq), &
                                      nsteps, order, three, yend, nff, nfg, status)
        end select
    end subroutine integrate_case

    ! |v|^2 is exactly conserved by the true flow, so the drift is pure
    ! truncation error and must decay at least as fast as the method order.
    ! A bare threshold would be arbitrary; the decay rate is the real claim.
    subroutine check_lorentz_energy(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: y0(6), yend(6), e0, drift(2), rate
        integer :: nff, nfg, i, ns(2), order, icase
        logical :: three
        type(fortnum_status_t) :: status
        character(len=28) :: label

        y0 = [0.3_dp, -0.2_dp, 0.1_dp, 0.4_dp, 0.9_dp, 0.35_dp]
        e0 = sum(y0(4:6)**2)
        ns = [2500, 5000]

        do icase = 1, 2
            select case (icase)
            case (1); order = 3; three = .false.; label = "RKN3 (2-stage)"
            case (2); order = 4; three = .false.; label = "RKN4 (2-stage)"
            end select
            do i = 1, 2
                call tdrk_integrate_fixed(f_lorentz, g_lorentz, 0.0_dp, 40.0_dp, y0, &
                                          ns(i), order, three, yend, nff, nfg, status)
                drift(i) = abs(sum(yend(4:6)**2) - e0)/e0
            end do
            rate = log(drift(1)/drift(2))/log(2.0_dp)
            write (*, "(a,a24,a,es10.3,a,f6.2)") "  ", label, "  |v|^2 drift ", &
                drift(2), "   decay rate ", rate
            if (rate < real(order, dp) - 0.4_dp) then
                write (error_unit, "(a,f6.2)") "  drift decays slower than order: ", &
                    rate
                nfail = nfail + 1
            end if
        end do
    end subroutine check_lorentz_energy

    ! RKNG on the same problem in genuine second-order form. Full orbit is the
    ! one place a Nystrom method applies with no reformulation. The order here
    ! is measured, not claimed: the tableau satisfies the special (velocity
    ! independent) conditions, and the general case has strictly more.
    subroutine check_rkng_lorentz(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: y0(3), yp0(3), yend(3), ypend(3), yref(6), z0(6)
        real(dp) :: err(3), ordobs, e0, drift
        integer :: nfev, i, ns(3)
        type(fortnum_status_t) :: status

        z0 = [0.3_dp, -0.2_dp, 0.1_dp, 0.4_dp, 0.9_dp, 0.35_dp]
        y0 = z0(1:3)
        yp0 = z0(4:6)
        e0 = sum(yp0**2)
        call rk4(3, 6, z0, 0.0_dp, 2.0_dp, 2000000, yref)

        ns = [400, 800, 1600]
        do i = 1, 3
            call rkng_integrate_fixed(f2_lorentz, 0.0_dp, 2.0_dp, y0, yp0, ns(i), &
                                      yend, ypend, nfev, status)
            err(i) = maxval(abs(yend - yref(1:3)))
        end do
        ordobs = log(err(2)/err(3))/log(2.0_dp)
        drift = abs(sum(ypend**2) - e0)/e0

        write (*, "(a,f6.2,a,es10.3,a,es10.3)") "  RKNG (2nd-order form) order ", &
            ordobs, "   err ", err(3), "   |v|^2 drift ", drift
        if (status%code /= FORTNUM_OK .or. .not. (ordobs > 1.7_dp)) then
            write (error_unit, "(a,f6.2)") "  RKNG order below 2: ", ordobs
            nfail = nfail + 1
        end if
        if (.not. (err(3) < 1.0e-3_dp)) then
            write (error_unit, "(a)") "  RKNG endpoint error too large"
            nfail = nfail + 1
        end if
    end subroutine check_rkng_lorentz

    ! ------------------------------------------------------------- error control
    !
    ! What error control has to deliver is not "the error is small" -- a fixed
    ! step can do that -- but that the achieved error TRACKS the requested
    ! tolerance. So the test measures the log-log slope of global error against
    ! rtol over four decades and requires it near 1. A stepper whose embedded
    ! estimate is broken (say, identically zero, or one order too high) still
    ! produces small errors at small rtol; it fails this.
    !
    ! Oracles: the closed form 1/(1+t) for the scalar problem, and a
    ! high-resolution classical RK4 -- a different method -- for Kepler.
    subroutine check_tdrk_adaptive(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: y0(4), yend(4), yref(4), hlast
        real(dp) :: tols(4), errs(4), slope, ecc
        integer  :: nff, nfg, nacc(4), nrej, i
        type(fortnum_status_t) :: status

        tols = [1.0e-4_dp, 1.0e-6_dp, 1.0e-8_dp, 1.0e-10_dp]

        ! --- nonlinear scalar against its closed form
        do i = 1, 4
            y0 = 0.0_dp
            y0(1) = 1.0_dp
            call tdrk_integrate_adaptive(f_scalar, g_scalar, 0.0_dp, 2.0_dp, &
                                         y0(1:1), tols(i), tols(i)*1.0e-3_dp, &
                                         0.0_dp, yend(1:1), hlast, nff, nfg, &
                                         nacc(i), nrej, status)
            if (status%code /= FORTNUM_OK) then
                write (error_unit, "(a)") "  adaptive TDRK failed on scalar problem"
                nfail = nfail + 1
                return
            end if
            errs(i) = abs(yend(1) - 1.0_dp/3.0_dp)
        end do
        slope = log10(errs(1)/errs(4))/log10(tols(1)/tols(4))
        write (*, "(a,f6.2,a,es10.3,a,i0)") &
            "  adaptive TDRK  scalar: err/tol slope ", slope, "   err ", errs(4), &
            "   steps ", nacc(4)
        if (.not. (slope > 0.6_dp .and. slope < 1.4_dp)) then
            write (error_unit, "(a,f6.2)") &
                "  achieved error does not track tolerance, slope ", slope
            nfail = nfail + 1
        end if
        if (.not. (nacc(4) > nacc(1))) then
            write (error_unit, "(a)") &
                "  tightening the tolerance did not cost more steps"
            nfail = nfail + 1
        end if

        ! --- Kepler against a high-resolution RK4 run
        ecc = 0.5_dp
        y0 = 0.0_dp
        y0(1) = 1.0_dp - ecc
        y0(4) = sqrt((1.0_dp + ecc)/(1.0_dp - ecc))
        call rk4(2, 4, y0, 0.0_dp, 3.0_dp, 2000000, yref)
        do i = 1, 4
            call tdrk_integrate_adaptive(f_kepler, g_kepler, 0.0_dp, 3.0_dp, y0, &
                                         tols(i), tols(i)*1.0e-3_dp, 0.0_dp, &
                                         yend, hlast, nff, nfg, nacc(i), nrej, &
                                         status)
            if (status%code /= FORTNUM_OK) then
                write (error_unit, "(a)") "  adaptive TDRK failed on Kepler"
                nfail = nfail + 1
                return
            end if
            errs(i) = maxval(abs(yend - yref))
        end do
        slope = log10(errs(1)/errs(4))/log10(tols(1)/tols(4))
        write (*, "(a,f6.2,a,es10.3,a,i0)") &
            "  adaptive TDRK  Kepler: err/tol slope ", slope, "   err ", errs(4), &
            "   steps ", nacc(4)
        if (.not. (slope > 0.6_dp .and. slope < 1.4_dp)) then
            write (error_unit, "(a,f6.2)") &
                "  Kepler error does not track tolerance, slope ", slope
            nfail = nfail + 1
        end if
    end subroutine check_tdrk_adaptive

    ! Same claim for RKNG, with |v|^2 conservation as the second, exact oracle:
    ! a static magnetic force does no work, so the drift is pure truncation
    ! error and must shrink with the requested tolerance.
    subroutine check_rkng_adaptive(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: y0(3), yp0(3), yend(3), ypend(3), yref(6), z0(6), hlast
        real(dp) :: tols(4), errs(4), drift(4), slope, e0
        integer  :: nfev, nacc(4), nrej, i
        type(fortnum_status_t) :: status

        z0 = [0.3_dp, -0.2_dp, 0.1_dp, 0.4_dp, 0.9_dp, 0.35_dp]
        y0 = z0(1:3)
        yp0 = z0(4:6)
        e0 = sum(yp0**2)
        call rk4(3, 6, z0, 0.0_dp, 2.0_dp, 2000000, yref)

        tols = [1.0e-4_dp, 1.0e-6_dp, 1.0e-8_dp, 1.0e-10_dp]
        do i = 1, 4
            call rkng_integrate_adaptive(f2_lorentz, 0.0_dp, 2.0_dp, y0, yp0, &
                                         tols(i), tols(i)*1.0e-3_dp, 0.0_dp, &
                                         yend, ypend, hlast, nfev, nacc(i), nrej, &
                                         status)
            if (status%code /= FORTNUM_OK) then
                write (error_unit, "(a)") "  adaptive RKNG failed"
                nfail = nfail + 1
                return
            end if
            errs(i) = maxval(abs(yend - yref(1:3)))
            drift(i) = abs(sum(ypend**2) - e0)/e0
        end do
        slope = log10(errs(1)/errs(4))/log10(tols(1)/tols(4))
        write (*, "(a,f6.2,a,es10.3,a,es10.3,a,i0)") &
            "  adaptive RKNG  Lorentz: err/tol slope ", slope, "   err ", errs(4), &
            "   |v|^2 drift ", drift(4), "   steps ", nacc(4)
        if (.not. (slope > 0.6_dp .and. slope < 1.4_dp)) then
            write (error_unit, "(a,f6.2)") &
                "  RKNG error does not track tolerance, slope ", slope
            nfail = nfail + 1
        end if
        if (.not. (drift(4) < drift(1))) then
            write (error_unit, "(a)") &
                "  |v|^2 drift did not improve with a tighter tolerance"
            nfail = nfail + 1
        end if
        if (.not. (nacc(4) > nacc(1))) then
            write (error_unit, "(a)") &
                "  RKNG: tightening the tolerance did not cost more steps"
            nfail = nfail + 1
        end if
    end subroutine check_rkng_adaptive

    ! The number the SIMPLE side needs: TDRK beats classical RK of equal order
    ! iff cost(G)/cost(F) is below the break-even printed here.
    subroutine report_cost_model()
        write (*, "(a)") "  cost model (TDRK 1F + sG vs classical RK s_rk F):"
        write (*, "(a)") "    order 3, 2 stages: 1F+2G vs 3F  =>  break-even rho 1.00"
        write (*, "(a)") "    order 4, 2 stages: 1F+2G vs 4F  =>  break-even rho 1.50"
        write (*, "(a)") "    order 4, 3 stages: 1F+3G vs 4F  =>  break-even rho 1.00"
    end subroutine report_cost_model

end program test_fortnum_ode_tdrk
