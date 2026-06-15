module fortnum_ode_ddeabm
    ! Variable-order (1..12) Adams-Bashforth-Moulton PECE integrator,
    ! functionally equivalent to SLATEC DDEABM. A modified divided-difference
    ! Adams predictor-corrector with local extrapolation, local error control
    ! per unit step, and order/step adaptation over orders 1 through 12. The
    ! integrator is re-entrant: the state type (ddeabm_state_t) carries the
    ! phi divided-difference array, the psi/alpha/beta/sig/g/v/w coefficient
    ! work, the order, and the step history across calls so a consumer can
    ! continue from grid point to grid point (DDEABM INFO(1)=1 continuation).
    ! An optional tstop bound forbids stepping past a set endpoint (DDEABM
    ! INFO(4)=1, RWORK(1)=tstop); output is delivered at tout (INFO(3)=0) by
    ! interpolating the Adams polynomial back onto the requested point.
    !
    ! Provenance: clean-room implementation of the published Shampine-Gordon
    !   variable-order Adams PECE algorithm.
    !   L. F. Shampine and M. K. Gordon, "Computer Solution of Ordinary
    !   Differential Equations: The Initial Value Problem", W. H. Freeman,
    !   1975 (the DE/STEP/INTRP structure); and L. F. Shampine and
    !   M. K. Gordon, "Solving Ordinary Differential Equations with ODE, STEP,
    !   and INTRP", Report SLA-73-1060, Sandia Laboratories, 1973. The smooth
    !   interpolant follows H. A. Watts, "A smoother interpolant for DE/STEP,
    !   INTRP II", Report SAND84-0293, Sandia Laboratories, 1984. Written
    !   natively from the algorithm description, MIT licensed,
    !   differentiation-ready style. No SLATEC ddeabm/ddes/dsteps/dintp source
    !   was copied; the control constants reproduce the published method.
    !
    ! Derivative policy: this module records no trace and exposes no AD path
    !   (the AD surface lives in fortnum_ode via the Cash-Karp trace). The
    !   public dispatch in fortnum_ode treats ddeabm as a primal-only method.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_ode, only: ode_rhs_t
    implicit none
    private

    public :: ddeabm_state_t
    public :: ddeabm_init, ddeabm_integrate_to
    public :: ode_solve_ddeabm
    public :: DDEABM_MAX_ORDER

    ! Maximum Adams order. The phi divided-difference array carries 16 columns
    ! (orders 1..12 plus the prediction/round-off bookkeeping columns 13..16).
    integer, parameter :: DDEABM_MAX_ORDER = 12

    ! gstr(k) error coefficients for the order-k Adams formula and the
    ! step-doubling caps two(k) = 2**k, both as published in Shampine-Gordon.
    real(dp), parameter :: GSTR(13) = [ &
        0.5_dp, 0.0833_dp, 0.0417_dp, 0.0264_dp, 0.0188_dp, 0.0143_dp, &
        0.0114_dp, 0.00936_dp, 0.00789_dp, 0.00679_dp, 0.00592_dp, &
        0.00524_dp, 0.00468_dp]
    real(dp), parameter :: TWO_POW(13) = [ &
        2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp, 32.0_dp, 64.0_dp, 128.0_dp, &
        256.0_dp, 512.0_dp, 1024.0_dp, 2048.0_dp, 4096.0_dp, 8192.0_dp]

    ! Re-entrant integrator state. Mirrors the SLATEC DDES/DSTEPS persistent
    ! work split across one type. phi(:,1..16) is the modified divided
    ! difference array; columns 13..16 hold prediction and round-off control.
    type :: ddeabm_state_t
        integer  :: neq    = 0
        integer  :: kord   = 1   ! order for the next step (K in DSTEPS)
        integer  :: kold   = 0   ! order used for the last successful step
        integer  :: kprev  = 0
        integer  :: ns     = 0   ! steps taken at the current step size
        integer  :: ksteps = 0   ! attempted steps since last reset
        integer  :: ivc    = 0
        integer  :: kgi    = 0
        integer  :: init   = 0   ! 0 fresh, 1 yp set, 2 fully started
        integer  :: nsteps = 0   ! accepted steps (diagnostic)
        integer  :: nfev   = 0
        integer  :: nrejected = 0
        integer  :: max_steps = 500000
        logical  :: start  = .true.
        logical  :: phase1 = .true.
        logical  :: nornd  = .true.
        logical  :: intout = .false.
        real(dp) :: x      = 0.0_dp ! internal mesh top
        real(dp) :: xold   = 0.0_dp
        real(dp) :: t      = 0.0_dp ! reported current point
        real(dp) :: told   = 0.0_dp
        real(dp) :: h      = 0.0_dp
        real(dp) :: hold   = 0.0_dp
        real(dp) :: eps    = 1.0_dp
        real(dp) :: delsgn = 0.0_dp
        real(dp) :: twou   = 0.0_dp
        real(dp) :: fouru  = 0.0_dp
        real(dp) :: tstop  = 0.0_dp
        logical  :: has_tstop = .false.
        real(dp), allocatable :: yy(:), yp(:), ypout(:), wt(:), p(:)
        real(dp), allocatable :: phi(:,:) ! neq x 16
        integer  :: iv(10)    = 0
        real(dp) :: psi(12)   = 0.0_dp
        real(dp) :: alpha(12) = 0.0_dp
        real(dp) :: beta(12)  = 0.0_dp
        real(dp) :: sig(13)   = 0.0_dp
        real(dp) :: v(12)     = 0.0_dp
        real(dp) :: w(12)     = 0.0_dp
        real(dp) :: g(13)     = 0.0_dp
        real(dp) :: gi(11)    = 0.0_dp
    end type ddeabm_state_t

contains

    ! Initialise the integrator state at (t0, y0). Allocates the divided
    ! difference array and work vectors. Tolerances and bounds are passed per
    ! integration call. Call once per new problem (DDEABM INFO(1)=0); a later
    ! ddeabm_integrate_to continues from the carried state (INFO(1)=1).
    subroutine ddeabm_init(state, neq, t0, y0)
        type(ddeabm_state_t), intent(inout) :: state
        integer,  intent(in) :: neq
        real(dp), intent(in) :: t0
        real(dp), intent(in) :: y0(:)

        real(dp) :: u

        u = epsilon(1.0_dp)
        state%neq    = neq
        state%kord   = 1
        state%kold   = 0
        state%kprev  = 0
        state%ns     = 0
        state%ksteps = 0
        state%ivc    = 0
        state%kgi    = 0
        state%init   = 0
        state%nsteps = 0
        state%nfev   = 0
        state%nrejected = 0
        state%start  = .true.
        state%phase1 = .true.
        state%nornd  = .true.
        state%intout = .false.
        state%x      = t0
        state%xold   = t0
        state%t      = t0
        state%told   = t0
        state%h      = 0.0_dp
        state%hold   = 0.0_dp
        state%eps    = 1.0_dp
        state%delsgn = 0.0_dp
        state%twou   = 2.0_dp * u
        state%fouru  = 4.0_dp * u
        state%tstop     = 0.0_dp
        state%has_tstop = .false.
        state%iv    = 0
        state%psi   = 0.0_dp
        state%alpha = 0.0_dp
        state%beta  = 0.0_dp
        state%sig   = 0.0_dp
        state%v     = 0.0_dp
        state%w     = 0.0_dp
        state%g     = 0.0_dp
        state%gi    = 0.0_dp

        call ensure_work(state, neq)
        state%yy(:) = y0(1:neq)
        state%phi(:, :) = 0.0_dp
    end subroutine ddeabm_init

    subroutine ensure_work(state, neq)
        type(ddeabm_state_t), intent(inout) :: state
        integer, intent(in) :: neq
        if (allocated(state%yy)) then
            if (size(state%yy) == neq) return
            deallocate(state%yy, state%yp, state%ypout, state%wt, state%p, &
                state%phi)
        end if
        allocate(state%yy(neq), state%yp(neq), state%ypout(neq), &
            state%wt(neq), state%p(neq), state%phi(neq, 16))
    end subroutine ensure_work

    ! Optional tstop bound (DDEABM INFO(4)=1, RWORK(1)=tstop): forbid the
    ! integrator from stepping past tstop. Applied from the ddeabm_integrate_to
    ! optional argument.
    subroutine ddeabm_set_tstop(state, tstop)
        type(ddeabm_state_t), intent(inout) :: state
        real(dp), intent(in) :: tstop
        state%tstop = tstop
        state%has_tstop = .true.
    end subroutine ddeabm_set_tstop

    ! Integrate the current state to tout (DDEABM ITASK interval mode,
    ! INFO(3)=0): take internal Adams steps until tout is reached or passed,
    ! then interpolate the Adams polynomial back to tout exactly. relerr is
    ! scalar; abserr is a vector of length neq or length 1 (scalar broadcast).
    ! On return y_out holds y(tout) and state continues from the internal mesh.
    !
    ! If tstop is present (or set on an earlier call), no step lands past it
    ! (DDEABM INFO(4)=1): the final approach to tstop is by extrapolation and
    ! the return is exactly at tout (which must satisfy |tout-t| <= |tstop-t|).
    subroutine ddeabm_integrate_to(rhs, state, tout, relerr, abserr, y_out, &
            status, tstop, ctx)
        procedure(ode_rhs_t)                 :: rhs
        type(ddeabm_state_t),   intent(inout) :: state
        real(dp), intent(in)                 :: tout
        real(dp), intent(in)                 :: relerr
        real(dp), intent(in)                 :: abserr(:)
        real(dp), allocatable, intent(out)   :: y_out(:)
        type(fortnum_status_t), intent(out)  :: status
        real(dp), intent(in), optional       :: tstop
        class(*), intent(in), optional       :: ctx

        integer  :: neq, l, ltol
        real(dp) :: del, absdel, dt, ha, wt_min
        logical  :: crash

        call status_set(status, FORTNUM_OK, "")
        neq = state%neq

        if (neq < 1 .or. .not. allocated(state%yy)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ddeabm_integrate_to: state not initialised")
            return
        end if
        if (relerr < 0.0_dp .or. any(abserr < 0.0_dp) .or. &
            (relerr == 0.0_dp .and. all(abserr == 0.0_dp))) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "ddeabm_integrate_to: tolerances must be nonnegative, " // &
                "not both 0")
            return
        end if

        if (present(tstop)) call ddeabm_set_tstop(state, tstop)

        allocate(y_out(neq))

        ! tstop / tout consistency: cannot ask to integrate past tstop.
        if (state%has_tstop) then
            if (sign(1.0_dp, tout - state%t) /= &
                    sign(1.0_dp, state%tstop - state%t) .or. &
                abs(tout - state%t) > abs(state%tstop - state%t)) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "ddeabm_integrate_to: tout conflicts with tstop")
                return
            end if
        end if

        ! Continuation guards (DDES): t == tout and direction reversals are
        ! disallowed once started.
        if (state%init /= 0) then
            if (state%t == tout) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "ddeabm_integrate_to: t == tout on continuation")
                return
            end if
            if (state%init == 2 .and. &
                state%delsgn * (tout - state%t) < 0.0_dp) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "ddeabm_integrate_to: direction reversal without restart")
                return
            end if
        end if

        ! INIT=0: evaluate the initial derivative. INIT=1: set the integration
        ! direction and nominal step. INIT=2: fully started, fall through.
        if (state%init == 0) then
            call rhs(state%t, state%yy, state%yp, ctx)
            state%nfev = state%nfev + 1
            state%init = 1
            if (state%t == tout) then
                y_out(:) = state%yy(:)
                state%told = state%t
                return
            end if
        end if
        if (state%init == 1) then
            state%x = state%t
            state%delsgn = sign(1.0_dp, tout - state%t)
            state%h = sign(max(state%fouru * abs(state%x), &
                abs(tout - state%x)), tout - state%x)
            state%init = 2
        end if

        del = tout - state%t
        absdel = abs(del)

        do
            ! Already past tout: interpolate the Adams polynomial and return.
            if (abs(state%x - state%t) >= absdel) then
                call interp(state, tout, y_out, state%ypout)
                state%t = tout
                state%told = state%t
                state%intout = .false.
                return
            end if

            ! Close to tstop and cannot step past it: extrapolate to tout.
            if (state%has_tstop .and. &
                abs(state%tstop - state%x) < state%fouru * abs(state%x)) then
                dt = tout - state%x
                y_out(:) = state%yy(:) + dt * state%yp(:)
                state%t = tout
                state%told = state%t
                return
            end if

            if (state%ksteps > state%max_steps) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "ddeabm_integrate_to: exceeded max_steps")
                state%t = state%x
                state%told = state%t
                y_out(:) = state%yy(:)
                return
            end if

            ! Limit the step to land no further than tstop, set the weights,
            ! then take one Adams step.
            ha = abs(state%h)
            if (state%has_tstop) ha = min(ha, abs(state%tstop - state%x))
            state%h = sign(ha, state%h)
            state%eps = 1.0_dp
            ltol = 1
            wt_min = huge(1.0_dp)
            do l = 1, neq
                if (size(abserr) /= 1) ltol = l
                state%wt(l) = relerr * abs(state%yy(l)) + abserr(ltol)
                wt_min = min(wt_min, state%wt(l))
            end do
            if (wt_min <= 0.0_dp) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "ddeabm_integrate_to: relative error criterion " // &
                    "inappropriate (zero weight)")
                state%t = state%x
                state%told = state%t
                y_out(:) = state%yy(:)
                return
            end if

            call step(rhs, state, crash, ctx)

            if (crash) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "ddeabm_integrate_to: tolerances too small for " // &
                    "machine precision (step crash)")
                state%t = state%x
                state%told = state%t
                y_out(:) = state%yy(:)
                return
            end if

            state%intout = .true.
        end do
    end subroutine ddeabm_integrate_to

    ! Flat one-shot entry matching ode_solve / ode_solve_vode: build state,
    ! integrate t0 -> t1, return the endpoint. For continuation and tstop use
    ! ddeabm_init + ddeabm_integrate_to directly (the re-entrant path KAMEL
    ! needs). rtol/atol default to the ode_problem_t defaults.
    subroutine ode_solve_ddeabm(rhs, t0, t1, y0, t_out, y_out, status, &
            rtol, atol)
        procedure(ode_rhs_t)               :: rhs
        real(dp),               intent(in) :: t0, t1
        real(dp),               intent(in) :: y0(:)
        real(dp), allocatable, intent(out) :: t_out(:)
        real(dp), allocatable, intent(out) :: y_out(:,:)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional     :: rtol, atol

        type(ddeabm_state_t) :: st
        real(dp), allocatable :: yend(:)
        real(dp) :: rt, at(1)

        rt = 1.0e-6_dp
        at = 1.0e-9_dp
        if (present(rtol)) rt = rtol
        if (present(atol)) at = atol

        call ddeabm_init(st, size(y0), t0, y0)
        call ddeabm_integrate_to(rhs, st, t1, rt, at, yend, status)

        allocate(t_out(2))
        allocate(y_out(size(y0), 2))
        t_out(1) = t0
        t_out(2) = t1
        y_out(:, 1) = y0
        y_out(:, 2) = yend
    end subroutine ode_solve_ddeabm

    ! --- DSTEPS-equivalent one-step Adams PECE ---

    ! Integrate one step from state%x to state%x + state%h with the modified
    ! divided-difference Adams formulas, local extrapolation, and order/step
    ! selection. On success state%x advances by state%hold, state%yy holds the
    ! new solution, state%yp the derivative there, and the order/step for the
    ! next step are chosen. crash = .true. if no step is possible (the error
    ! tolerance is below the machine round-off floor); state%eps is then raised
    ! to an acceptable value and nothing is advanced.
    subroutine step(rhs, state, crash, ctx)
        procedure(ode_rhs_t)            :: rhs
        type(ddeabm_state_t), intent(inout) :: state
        logical, intent(out)            :: crash
        class(*), intent(in), optional  :: ctx

        integer  :: neq, k, kp1, kp2, km1, km2, knew, ifail
        integer  :: i, ip1, im1, iq, j, l, jv, nsp1, nsp2, nsm2, limit1, limit2
        real(dp) :: p5eps, round, absh, hnew, err, erk, erkm1, erkm2, erkp1
        real(dp) :: temp1, temp2, temp3, temp4, temp5, temp6, tau, rho, r
        real(dp) :: realns, reali

        neq = state%neq
        crash = .true.

        ! Block 0: round-off / tolerance floor checks; first-step init.
        if (abs(state%h) < state%fouru * abs(state%x)) then
            state%h = sign(state%fouru * abs(state%x), state%h)
            return
        end if
        p5eps = 0.5_dp * state%eps
        round = 0.0_dp
        do l = 1, neq
            round = round + (state%yy(l) / state%wt(l))**2
        end do
        round = state%twou * sqrt(round)
        if (p5eps < round) then
            state%eps = 2.0_dp * round * (1.0_dp + state%fouru)
            return
        end if
        crash = .false.
        state%g(1) = 1.0_dp
        state%g(2) = 0.5_dp
        state%sig(1) = 1.0_dp

        if (state%start) then
            do l = 1, neq
                state%phi(l, 1) = state%yp(l)
                state%phi(l, 2) = 0.0_dp
            end do
            call hstart(rhs, state, ctx)
            state%hold = 0.0_dp
            state%kord = 1
            state%kold = 0
            state%kprev = 0
            state%start = .false.
            state%phase1 = .true.
            state%nornd = .true.
            if (p5eps <= 100.0_dp * round) then
                state%nornd = .false.
                state%phi(:, 15) = 0.0_dp
            end if
        end if
        ifail = 0

        ! The retry loop after an unsuccessful step (block 3 -> block 1).
        do
            k = state%kord
            kp1 = k + 1
            kp2 = k + 2
            km1 = k - 1
            km2 = k - 2

            ! Block 1: coefficients for this step. ns counts steps at the
            ! current size; when k < ns no coefficient changes.
            if (state%h /= state%hold) state%ns = 0
            if (state%ns <= state%kold) state%ns = state%ns + 1
            nsp1 = state%ns + 1

            if (k >= state%ns) then
                state%beta(state%ns) = 1.0_dp
                realns = real(state%ns, dp)
                state%alpha(state%ns) = 1.0_dp / realns
                temp1 = state%h * realns
                state%sig(nsp1) = 1.0_dp
                if (k >= nsp1) then
                    do i = nsp1, k
                        im1 = i - 1
                        temp2 = state%psi(im1)
                        state%psi(im1) = temp1
                        state%beta(i) = state%beta(im1) * state%psi(im1) / temp2
                        temp1 = temp2 + state%h
                        state%alpha(i) = state%h / temp1
                        reali = real(i, dp)
                        state%sig(i + 1) = reali * state%alpha(i) * state%sig(i)
                    end do
                end if
                state%psi(k) = temp1

                ! Coefficients g(*) via the v(*)/w(*) work vectors.
                if (state%ns <= 1) then
                    do iq = 1, k
                        state%v(iq) = 1.0_dp / real(iq * (iq + 1), dp)
                        state%w(iq) = state%v(iq)
                    end do
                    state%ivc = 0
                    state%kgi = 0
                    if (k /= 1) then
                        state%kgi = 1
                        state%gi(1) = state%w(2)
                    end if
                else
                    ! If order was raised, update the diagonal part of v(*).
                    if (k > state%kprev) then
                        if (state%ivc /= 0) then
                            jv = kp1 - state%iv(state%ivc)
                            state%ivc = state%ivc - 1
                        else
                            jv = 1
                            state%v(k) = 1.0_dp / real(k * kp1, dp)
                            state%w(k) = state%v(k)
                            if (k == 2) then
                                state%kgi = 1
                                state%gi(1) = state%w(2)
                            end if
                        end if
                        nsm2 = state%ns - 2
                        if (nsm2 >= jv) then
                            i = 0
                            do j = jv, nsm2
                                i = k - j
                                state%v(i) = state%v(i) &
                                    - state%alpha(j + 1) * state%v(i + 1)
                                state%w(i) = state%v(i)
                            end do
                            if (i == 2) then
                                state%kgi = state%ns - 1
                                state%gi(state%kgi) = state%w(2)
                            end if
                        end if
                    end if

                    ! Update v(*) and set w(*).
                    limit1 = kp1 - state%ns
                    temp5 = state%alpha(state%ns)
                    do iq = 1, limit1
                        state%v(iq) = state%v(iq) - temp5 * state%v(iq + 1)
                        state%w(iq) = state%v(iq)
                    end do
                    state%g(nsp1) = state%w(1)
                    if (limit1 /= 1) then
                        state%kgi = state%ns
                        state%gi(state%kgi) = state%w(2)
                    end if
                    state%w(limit1 + 1) = state%v(limit1 + 1)
                    if (k < state%kold) then
                        state%ivc = state%ivc + 1
                        state%iv(state%ivc) = limit1 + 2
                    end if
                end if

                ! Compute the remaining g(*) in the work vector w(*).
                nsp2 = state%ns + 2
                state%kprev = k
                if (kp1 >= nsp2) then
                    do i = nsp2, kp1
                        limit2 = kp2 - i
                        temp6 = state%alpha(i - 1)
                        do iq = 1, limit2
                            state%w(iq) = state%w(iq) - temp6 * state%w(iq + 1)
                        end do
                        state%g(i) = state%w(1)
                    end do
                end if
            end if

            ! Block 2: predict, evaluate, estimate errors at k, k-1, k-2.
            state%ksteps = state%ksteps + 1

            if (k >= nsp1) then
                do i = nsp1, k
                    temp1 = state%beta(i)
                    state%phi(:, i) = temp1 * state%phi(:, i)
                end do
            end if

            do l = 1, neq
                state%phi(l, kp2) = state%phi(l, kp1)
                state%phi(l, kp1) = 0.0_dp
                state%p(l) = 0.0_dp
            end do
            do j = 1, k
                i = kp1 - j
                ip1 = i + 1
                temp2 = state%g(i)
                do l = 1, neq
                    state%p(l) = state%p(l) + temp2 * state%phi(l, i)
                    state%phi(l, i) = state%phi(l, i) + state%phi(l, ip1)
                end do
            end do
            if (state%nornd) then
                do l = 1, neq
                    state%p(l) = state%yy(l) + state%h * state%p(l)
                end do
            else
                do l = 1, neq
                    tau = state%h * state%p(l) - state%phi(l, 15)
                    state%p(l) = state%yy(l) + tau
                    state%phi(l, 16) = (state%p(l) - state%yy(l)) - tau
                end do
            end if
            state%xold = state%x
            state%x = state%x + state%h
            absh = abs(state%h)
            call rhs(state%x, state%p, state%yp, ctx)
            state%nfev = state%nfev + 1

            erkm2 = 0.0_dp
            erkm1 = 0.0_dp
            erk = 0.0_dp
            do l = 1, neq
                temp3 = 1.0_dp / state%wt(l)
                temp4 = state%yp(l) - state%phi(l, 1)
                if (km2 > 0) erkm2 = erkm2 &
                    + ((state%phi(l, km1) + temp4) * temp3)**2
                if (km2 >= 0) erkm1 = erkm1 &
                    + ((state%phi(l, k) + temp4) * temp3)**2
                erk = erk + (temp4 * temp3)**2
            end do
            if (km2 > 0) erkm2 = absh * state%sig(km1) * GSTR(km2) * sqrt(erkm2)
            if (km2 >= 0) erkm1 = absh * state%sig(k) * GSTR(km1) * sqrt(erkm1)
            temp5 = absh * sqrt(erk)
            err = temp5 * (state%g(k) - state%g(kp1))
            erk = temp5 * state%sig(kp1) * GSTR(k)
            knew = k

            ! Test if the order should be lowered.
            if (km2 > 0) then
                if (max(erkm1, erkm2) <= erk) knew = km1
            else if (km2 == 0) then
                if (erkm1 <= 0.5_dp * erk) knew = km1
            end if

            ! Block 2 exit: accept if err <= eps, else block 3.
            if (err <= state%eps) exit

            ! Block 3: unsuccessful step. Restore x, phi, psi.
            state%phase1 = .false.
            state%x = state%xold
            do i = 1, k
                temp1 = 1.0_dp / state%beta(i)
                ip1 = i + 1
                state%phi(:, i) = temp1 * (state%phi(:, i) - state%phi(:, ip1))
            end do
            if (k >= 2) then
                do i = 2, k
                    state%psi(i - 1) = state%psi(i) - state%h
                end do
            end if

            ! On the third failure set order to one; thereafter use an optimal
            ! step. Double the tolerance and crash if h underflows.
            ifail = ifail + 1
            temp2 = 0.5_dp
            if (ifail > 3) then
                if (p5eps < 0.25_dp * erk) temp2 = sqrt(p5eps / erk)
            end if
            if (ifail >= 3) knew = 1
            state%h = temp2 * state%h
            state%kord = knew
            state%ns = 0
            if (abs(state%h) < state%fouru * abs(state%x)) then
                crash = .true.
                state%h = sign(state%fouru * abs(state%x), state%h)
                state%eps = state%eps + state%eps
                return
            end if
            ! retry: cycle with the new order/step.
        end do

        ! Block 4: successful step. Correct, evaluate, update differences,
        ! then choose order and step for the next step.
        k = state%kord
        kp1 = k + 1
        kp2 = k + 2
        km1 = k - 1
        state%kold = k
        state%hold = state%h

        temp1 = state%h * state%g(kp1)
        if (state%nornd) then
            do l = 1, neq
                temp3 = state%yy(l)
                state%yy(l) = state%p(l) &
                    + temp1 * (state%yp(l) - state%phi(l, 1))
                state%p(l) = temp3
            end do
        else
            do l = 1, neq
                temp3 = state%yy(l)
                rho = temp1 * (state%yp(l) - state%phi(l, 1)) - state%phi(l, 16)
                state%yy(l) = state%p(l) + rho
                state%phi(l, 15) = (state%yy(l) - state%p(l)) - rho
                state%p(l) = temp3
            end do
        end if
        call rhs(state%x, state%yy, state%yp, ctx)
        state%nfev = state%nfev + 1

        do l = 1, neq
            state%phi(l, kp1) = state%yp(l) - state%phi(l, 1)
            state%phi(l, kp2) = state%phi(l, kp1) - state%phi(l, kp2)
        end do
        do i = 1, k
            state%phi(:, i) = state%phi(:, i) + state%phi(:, kp1)
        end do

        ! Estimate the error at order k+1 (unless in phase 1, an order drop is
        ! already chosen, or the step size is not constant), then pick order.
        erkp1 = 0.0_dp
        if (knew == km1 .or. k == DDEABM_MAX_ORDER) state%phase1 = .false.

        if (state%phase1) then
            ! Phase 1: always raise order.
            state%kord = kp1
            erk = erkp1
        else if (knew == km1) then
            ! Order already chosen to drop.
            state%kord = km1
            erk = erkm1
        else if (kp1 <= state%ns) then
            do l = 1, neq
                erkp1 = erkp1 + (state%phi(l, kp2) / state%wt(l))**2
            end do
            erkp1 = absh * GSTR(kp1) * sqrt(erkp1)
            if (k > 1) then
                if (erkm1 <= min(erk, erkp1)) then
                    state%kord = km1
                    erk = erkm1
                else if (erkp1 < erk .and. k /= DDEABM_MAX_ORDER) then
                    state%kord = kp1
                    erk = erkp1
                end if
            else
                if (erkp1 < 0.5_dp * erk) then
                    state%kord = kp1
                    erk = erkp1
                end if
            end if
        end if
        ! else: keep the current order (kp1 > ns, estimate unreliable).

        ! With the new order, determine the step size for the next step.
        hnew = state%h + state%h
        if (.not. state%phase1) then
            if (p5eps < erk * TWO_POW(state%kord + 1)) then
                hnew = state%h
                if (p5eps < erk) then
                    temp2 = real(state%kord + 1, dp)
                    r = (p5eps / erk)**(1.0_dp / temp2)
                    hnew = absh * max(0.5_dp, min(0.9_dp, r))
                    hnew = sign(max(hnew, state%fouru * abs(state%x)), state%h)
                end if
            end if
        end if
        state%h = hnew
        state%nsteps = state%nsteps + 1
    end subroutine step

    ! --- DHSTRT-equivalent starting step ---

    ! Compute a starting step size from the local Lipschitz constant, a bound
    ! on the first derivative, and a bound on the partial w.r.t. the
    ! independent variable, all approximated near (x, yy). Sets state%h.
    ! Uses the max norm (DHVNRM). On entry state%yp holds f(x, yy).
    subroutine hstart(rhs, state, ctx)
        procedure(ode_rhs_t)            :: rhs
        type(ddeabm_state_t), intent(inout) :: state
        class(*), intent(in), optional  :: ctx

        integer  :: neq, j, k, lk, morder
        real(dp) :: a, b, dx, absdx, relper, da, delf, dfdxb, fbnd, dely
        real(dp) :: dfdub, ydpb, tolmin, tolsum, tolexp, tolp, srydpb, h
        real(dp) :: small, big, dy
        real(dp), allocatable :: y(:), yprime(:), spy(:), pv(:), yp(:), sf(:)

        neq = state%neq
        morder = 1
        small = epsilon(1.0_dp)
        big = sqrt(huge(1.0_dp))
        a = state%x
        b = state%x + state%h
        allocate(y(neq), yprime(neq), spy(neq), pv(neq), yp(neq), sf(neq))
        y(:) = state%yy(:)
        yprime(:) = state%yp(:)

        dx = b - a
        absdx = abs(dx)
        relper = small**0.375_dp

        ! Bound (dfdxb) on the partial of f w.r.t. the independent variable,
        ! and a local bound (fbnd) on the first derivative.
        da = sign(max(min(relper * abs(a), absdx), &
            100.0_dp * small * abs(a)), dx)
        if (da == 0.0_dp) da = relper * dx
        call rhs(a + da, y, sf, ctx)
        state%nfev = state%nfev + 1
        do j = 1, neq
            yp(j) = sf(j) - yprime(j)
        end do
        delf = maxnorm(yp)
        dfdxb = big
        if (delf < big * abs(da)) dfdxb = delf / abs(da)
        fbnd = maxnorm(sf)

        ! Estimate (dfdub) the local Lipschitz constant by numerical
        ! differences. The perturbation size is held constant across iterations.
        dely = relper * maxnorm(y)
        if (dely == 0.0_dp) dely = relper
        dely = sign(dely, dx)
        delf = maxnorm(yprime)
        fbnd = max(fbnd, delf)
        if (delf == 0.0_dp) then
            do j = 1, neq
                spy(j) = 0.0_dp
                yp(j) = 1.0_dp
            end do
            delf = maxnorm(yp)
        else
            do j = 1, neq
                spy(j) = yprime(j)
                yp(j) = yprime(j)
            end do
        end if

        dfdub = 0.0_dp
        lk = min(neq + 1, 3)
        do k = 1, lk
            do j = 1, neq
                pv(j) = y(j) + dely * (yp(j) / delf)
            end do
            if (k == 2) then
                call rhs(a + da, pv, yp, ctx)
                state%nfev = state%nfev + 1
                do j = 1, neq
                    pv(j) = yp(j) - sf(j)
                end do
            else
                call rhs(a, pv, yp, ctx)
                state%nfev = state%nfev + 1
                do j = 1, neq
                    pv(j) = yp(j) - yprime(j)
                end do
            end if
            fbnd = max(fbnd, maxnorm(yp))
            delf = maxnorm(pv)
            if (delf >= big * abs(dely)) then
                dfdub = big
                exit
            end if
            dfdub = max(dfdub, delf / abs(dely))
            if (k == lk) exit
            ! Choose the next perturbation vector.
            if (delf == 0.0_dp) delf = 1.0_dp
            do j = 1, neq
                if (k == 2) then
                    dy = y(j)
                    if (dy == 0.0_dp) dy = dely / relper
                else
                    dy = abs(pv(j))
                    if (dy == 0.0_dp) dy = delf
                end if
                if (spy(j) == 0.0_dp) spy(j) = yp(j)
                if (spy(j) /= 0.0_dp) dy = sign(dy, spy(j))
                yp(j) = dy
            end do
            delf = maxnorm(yp)
        end do

        ! Bound (ydpb) on the norm of the second derivative.
        ydpb = dfdxb + dfdub * fbnd

        ! Tolerance parameter for the starting step (mid error-tolerance range).
        tolmin = big
        tolsum = 0.0_dp
        do k = 1, neq
            tolexp = log10(state%wt(k))
            tolmin = min(tolmin, tolexp)
            tolsum = tolsum + tolexp
        end do
        tolp = 10.0_dp**(0.5_dp * (tolsum / real(neq, dp) + tolmin) &
            / real(morder + 1, dp))

        ! Starting step from the first/second derivative information.
        h = absdx
        if (ydpb == 0.0_dp .and. fbnd == 0.0_dp) then
            if (tolp < 1.0_dp) h = absdx * tolp
        else if (ydpb == 0.0_dp) then
            if (tolp < fbnd * absdx) h = tolp / fbnd
        else
            srydpb = sqrt(0.5_dp * ydpb)
            if (tolp < srydpb * absdx) h = tolp / srydpb
        end if

        if (h * dfdub > 1.0_dp) h = 1.0_dp / dfdub
        h = max(h, 100.0_dp * small * abs(a))
        if (h == 0.0_dp) h = small * abs(b)
        state%h = sign(h, dx)
    end subroutine hstart

    ! --- DINTP-equivalent smooth interpolation ---

    ! Approximate the solution (yout) and derivative (ypout) at xout by
    ! evaluating the Adams polynomial built in step(), valid on the last
    ! completed interval [xold, x]. Follows the Watts smoother interpolant.
    subroutine interp(state, xout, yout, ypout)
        type(ddeabm_state_t), intent(in)  :: state
        real(dp),             intent(in)  :: xout
        real(dp),             intent(out) :: yout(:)
        real(dp),             intent(out) :: ypout(:)

        integer  :: neq, kold, kp1, kp2, i, iq, iw, j, jq, l, m
        real(dp) :: hi, h, xi, xim1, xiq, temp1, gdi, alp, gamma, sigma
        real(dp) :: rmu, hmu, gdif, temp2, temp3
        real(dp) :: g(13), c(13), w(13)

        neq = state%neq
        kold = state%kold
        kp1 = kold + 1
        kp2 = kold + 2

        hi = xout - state%xold
        h = state%x - state%xold
        xi = hi / h
        xim1 = xi - 1.0_dp

        ! Initialise w(*) for computing g(*).
        xiq = xi
        temp1 = 1.0_dp
        do iq = 1, kp1
            xiq = xi * xiq
            temp1 = real(iq * (iq + 1), dp)
            w(iq) = xiq / temp1
        end do

        ! Double-integral term gdi.
        if (kold <= state%kgi) then
            gdi = state%gi(kold)
        else
            if (state%ivc > 0) then
                iw = state%iv(state%ivc)
                gdi = state%w(iw)
                m = kold - iw + 3
            else
                gdi = 1.0_dp / temp1
                m = 2
            end if
            if (m <= kold) then
                do i = m, kold
                    gdi = state%w(kp2 - i) - state%alpha(i) * gdi
                end do
            end if
        end if

        ! Compute g(*) and c(*).
        g(1) = xi
        g(2) = 0.5_dp * xi * xi
        c(1) = 1.0_dp
        c(2) = xi
        if (kold >= 2) then
            do i = 2, kold
                alp = state%alpha(i)
                gamma = 1.0_dp + xim1 * alp
                l = kp2 - i
                do jq = 1, l
                    w(jq) = gamma * w(jq) - alp * w(jq + 1)
                end do
                g(i + 1) = w(1)
                c(i + 1) = gamma * c(i)
            end do
        end if

        ! Interpolation parameters.
        sigma = (w(2) - xim1 * w(1)) / gdi
        rmu = xim1 * c(kp1) / gdi
        hmu = rmu / h

        do l = 1, neq
            yout(l) = 0.0_dp
            ypout(l) = 0.0_dp
        end do
        do j = 1, kold
            i = kp2 - j
            gdif = state%g(i) - state%g(i - 1)
            temp2 = (g(i) - g(i - 1)) - sigma * gdif
            temp3 = (c(i) - c(i - 1)) + rmu * gdif
            do l = 1, neq
                yout(l) = yout(l) + temp2 * state%phi(l, i)
                ypout(l) = ypout(l) + temp3 * state%phi(l, i)
            end do
        end do
        ! OY = solution at xold (held in state%p after a successful step),
        ! Y = solution at x (state%yy). DINTP blends the two endpoints.
        do l = 1, neq
            yout(l) = ((1.0_dp - sigma) * state%p(l) + sigma * state%yy(l)) &
                + h * (yout(l) + (g(1) - sigma * state%g(1)) * state%phi(l, 1))
            ypout(l) = hmu * (state%p(l) - state%yy(l)) &
                + (ypout(l) + (c(1) + rmu * state%g(1)) * state%phi(l, 1))
        end do
    end subroutine interp

    ! Max norm over a vector (DHVNRM): the largest component magnitude.
    pure real(dp) function maxnorm(v) result(nrm)
        real(dp), intent(in) :: v(:)
        nrm = maxval(abs(v))
    end function maxnorm

end module fortnum_ode_ddeabm
