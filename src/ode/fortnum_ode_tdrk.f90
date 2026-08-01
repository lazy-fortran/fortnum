module fortnum_ode_tdrk
    ! Two-derivative Runge-Kutta (TDRK) for first-order systems, and general
    ! Runge-Kutta-Nystrom (RKNG) for second-order systems.
    !
    ! ---------------------------------------------------------------------
    ! Why this module exists
    !
    ! Runge-Kutta-Nystrom methods integrate y'' = f(t, y) and are the reason
    ! high-order explicit integrators do so well in celestial mechanics: they
    ! skip the velocity stages, and their order conditions (Nystrom trees) are
    ! far sparser than general Runge-Kutta trees, so a given order needs far
    ! fewer stages. Charged-particle motion -- guiding-centre or full-orbit --
    ! is first order and non-separable, so RKN does not apply to it as written.
    !
    ! It applies after one differentiation. For zdot = F(z),
    !
    !     zddot = F'(z) F(z) =: G(z)
    !
    ! and G depends on z alone, because zdot is itself determined by z. So
    ! zddot = G(z) is in *special* second-order form, the form RKN requires.
    ! RKN carries (z, v) as independent variables while the true flow satisfies
    ! v = F(z); re-synchronising v := F(z) at the start of every step removes
    ! the drift and preserves the order, since the local error bound assumes an
    ! exact (z_n, v_n) to begin with.
    !
    ! The re-synchronised scheme is exactly a special two-derivative Runge-Kutta
    ! method (Chan & Tsai, Numer. Algorithms 53 (2010) 171-194,
    ! doi:10.1007/s11075-009-9349-1):
    !
    !     Z_i   = z + c_i h F(z) + h^2 sum_j abar_ij G(Z_j)
    !     z_new = z +     h F(z) + h^2 sum_i bbar_i  G(Z_i)
    !
    ! F is evaluated once per step, G at each stage. Because the velocity is
    ! discarded and recomputed rather than propagated, only the POSITION order
    ! conditions are needed -- a strictly smaller set than a full RKN pair
    ! requires. That is why the order-4 tableau below needs only two stages.
    !
    ! Whether this pays on a given problem is an empirical question, not a
    ! theoretical one. A TDRK step costs one F plus s evaluations of G; a
    ! classical explicit RK step of the same order costs s_rk evaluations of F.
    ! So TDRK wins iff
    !
    !     1 + s * rho  <  s_rk,     rho = cost(G)/cost(F).
    !
    ! For the order-4 tableau here, s = 2 and s_rk = 4, so break-even is
    ! rho = 1.5.
    !
    ! ---------------------------------------------------------------------
    ! Tableau provenance
    !
    ! The tableaux are CONSTRUCTED from the Nystrom position order conditions
    ! in tdrk_tableau, not transcribed. With the stage-consistency assumption
    ! sum_j abar_ij = c_i^2 / 2, the position conditions through order 4 reduce
    ! to quadrature conditions:
    !
    !     sum b = 1/2,   sum b c = 1/6,   sum b c^2 = 1/12
    !
    ! (the h^4 condition sum_i b_i sum_j abar_ij = 1/24 is then implied). The
    ! conditions themselves are re-checked in the test, so the coefficients
    ! carry no transcription risk.
    !
    ! ---------------------------------------------------------------------
    ! RKNG
    !
    ! Full-orbit motion is genuinely y'' = f(t, y, y') -- the Lorentz force
    ! depends on velocity -- which is *general* RKN form. That needs no
    ! reformulation at all, and is the one place a Nystrom-family method
    ! applies to this project's systems directly.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
                              FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_ode, only: ode_rhs_t

    implicit none
    private

    public :: tdrk_tableau, tdrk_step, tdrk_integrate_fixed
    public :: rkng_step, rkng_integrate_fixed
    public :: tdrk_integrate_adaptive, rkng_integrate_adaptive
    public :: ode_rhs2_t, ode_rhs_rkng_t
    public :: TDRK_MAX_STAGES

    integer, parameter :: TDRK_MAX_STAGES = 4

    ! Step-size controller (Hairer, Norsett, Wanner, "Solving ODEs I", II.4).
    ! The exponent is set per driver from the order of the EMBEDDED solution,
    ! not the propagated one, so it is passed in rather than fixed here.
    real(dp), parameter :: SAFETY  = 0.9_dp
    real(dp), parameter :: FAC_MIN = 0.2_dp
    real(dp), parameter :: FAC_MAX = 5.0_dp
    integer,  parameter :: MAX_STEPS = 10000000

    ! Second derivative of the state along the flow, G = F'(z) F(z), for the
    ! first-order system zdot = F(z). Same shape as ode_rhs_t.
    abstract interface
        subroutine ode_rhs2_t(t, y, g, ctx)
            import :: dp
            real(dp), intent(in)  :: t
            real(dp), intent(in)  :: y(:)
            real(dp), intent(out) :: g(:)
            class(*), intent(in), optional :: ctx
        end subroutine ode_rhs2_t
    end interface

    ! General second-order RHS: y'' = f(t, y, y'), the full-orbit case.
    abstract interface
        subroutine ode_rhs_rkng_t(t, y, yp, ypp, ctx)
            import :: dp
            real(dp), intent(in)  :: t
            real(dp), intent(in)  :: y(:)
            real(dp), intent(in)  :: yp(:)
            real(dp), intent(out) :: ypp(:)
            class(*), intent(in), optional :: ctx
        end subroutine ode_rhs_rkng_t
    end interface

contains

    ! Build a special-RKN tableau of the requested position order.
    !
    ! order 3: c = (0, 2/3),      b = (1/4, 1/4)
    ! order 4: c = (0, 1/2),      b = (1/6, 1/3)
    ! order 4 (3-stage, coupled): c = (0, 1/3, 2/3), b = (1/8, 1/4, 1/8)
    !
    ! The three-stage variant exists so the abar(3,2) coupling path is
    ! exercised; the two-stage variants never touch it.
    subroutine tdrk_tableau(order, three_stage, s, c, abar, bbar, ok)
        integer, intent(in) :: order
        logical, intent(in) :: three_stage
        integer, intent(out) :: s
        real(dp), intent(out) :: c(:), abar(:, :), bbar(:)
        logical, intent(out) :: ok

        c = 0.0_dp
        abar = 0.0_dp
        bbar = 0.0_dp
        ok = .true.

        if (order == 3) then
            s = 2
            c(2) = 2.0_dp/3.0_dp
            abar(2, 1) = 0.5_dp*c(2)**2
            bbar(1) = 0.25_dp
            bbar(2) = 0.25_dp
        else if (order == 4 .and. .not. three_stage) then
            s = 2
            c(2) = 0.5_dp
            abar(2, 1) = 0.5_dp*c(2)**2
            bbar(1) = 1.0_dp/6.0_dp
            bbar(2) = 1.0_dp/3.0_dp
        else if (order == 4) then
            s = 3
            c(2) = 1.0_dp/3.0_dp
            c(3) = 2.0_dp/3.0_dp
            abar(2, 1) = 0.5_dp*c(2)**2
            abar(3, 2) = 0.5_dp*c(3)**2
            bbar(1) = 0.125_dp
            bbar(2) = 0.25_dp
            bbar(3) = 0.125_dp
        else
            s = 0
            ok = .false.
        end if
    end subroutine tdrk_tableau

    ! One TDRK step for zdot = F(z), using G = F'F.
    !
    ! rhs is evaluated once (the velocity re-sync); rhs2 at each stage.
    subroutine tdrk_step(rhs, rhs2, t, y, h, s, c, abar, bbar, gs, ytmp, ynew, &
                         nfev_f, nfev_g, ctx)
        procedure(ode_rhs_t)  :: rhs
        procedure(ode_rhs2_t) :: rhs2
        real(dp), intent(in) :: t, h
        real(dp), intent(in) :: y(:)
        integer,  intent(in) :: s
        real(dp), intent(in) :: c(:), abar(:, :), bbar(:)
        real(dp), intent(inout) :: gs(:, :)      ! (TDRK_MAX_STAGES, neq)
        real(dp), intent(inout) :: ytmp(:)
        real(dp), intent(out) :: ynew(:)
        integer, intent(inout) :: nfev_f, nfev_g
        class(*), intent(in), optional :: ctx

        real(dp) :: fz(size(y)), acc
        integer  :: i, j, k, neq

        neq = size(y)

        call rhs(t, y, fz, ctx)
        nfev_f = nfev_f + 1

        do i = 1, s
            do k = 1, neq
                acc = 0.0_dp
                do j = 1, i - 1
                    acc = acc + abar(i, j)*gs(j, k)
                end do
                ytmp(k) = y(k) + c(i)*h*fz(k) + h*h*acc
            end do
            call rhs2(t + c(i)*h, ytmp(1:neq), gs(i, 1:neq), ctx)
            nfev_g = nfev_g + 1
        end do

        do k = 1, neq
            acc = 0.0_dp
            do i = 1, s
                acc = acc + bbar(i)*gs(i, k)
            end do
            ynew(k) = y(k) + h*fz(k) + h*h*acc
        end do
    end subroutine tdrk_step

    ! Fixed-step TDRK driver. Fixed step is what makes the order measurable, so
    ! this stays even though tdrk_integrate_adaptive exists: the order test
    ! needs a known step sequence, and a caller comparing against fixed-step
    ! symplectic schemes wants resolution, not tolerance, as the knob.
    subroutine tdrk_integrate_fixed(rhs, rhs2, t0, t1, y0, nsteps, order, &
                                    three_stage, yend, nfev_f, nfev_g, status, ctx)
        procedure(ode_rhs_t)  :: rhs
        procedure(ode_rhs2_t) :: rhs2
        real(dp), intent(in) :: t0, t1
        real(dp), intent(in) :: y0(:)
        integer,  intent(in) :: nsteps, order
        logical,  intent(in) :: three_stage
        real(dp), intent(out) :: yend(:)
        integer,  intent(out) :: nfev_f, nfev_g
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: ctx

        real(dp) :: c(TDRK_MAX_STAGES), abar(TDRK_MAX_STAGES, TDRK_MAX_STAGES)
        real(dp) :: bbar(TDRK_MAX_STAGES)
        real(dp), allocatable :: gs(:, :), ytmp(:), y(:), ynew(:)
        real(dp) :: h, t
        integer  :: s, k, neq
        logical  :: ok

        call status_set(status, FORTNUM_OK, "")
        if (nsteps < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "tdrk_integrate_fixed: nsteps must be positive")
            return
        end if

        call tdrk_tableau(order, three_stage, s, c, abar, bbar, ok)
        if (.not. ok) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "tdrk_integrate_fixed: no tableau for that order")
            return
        end if

        neq = size(y0)
        allocate (gs(TDRK_MAX_STAGES, neq), ytmp(neq), y(neq), ynew(neq))
        gs = 0.0_dp
        y = y0
        h = (t1 - t0)/real(nsteps, dp)
        nfev_f = 0
        nfev_g = 0

        t = t0
        do k = 1, nsteps
            call tdrk_step(rhs, rhs2, t, y, h, s, c, abar, bbar, gs, ytmp, ynew, &
                           nfev_f, nfev_g, ctx)
            y = ynew
            t = t0 + real(k, dp)*h
        end do
        yend(1:neq) = y
    end subroutine tdrk_integrate_fixed

    ! One general RKN step for y'' = f(t, y, y'). Unlike TDRK this propagates
    ! the velocity, because for a genuine second-order system the velocity is an
    ! independent state variable rather than a function of position.
    !
    ! Stage positions use the same c/abar as the special tableau. Stage
    ! velocities use avel, whose row sums equal c_i (the consistency condition);
    ! the step velocity uses bvel, from sum bvel = 1, sum bvel c = 1/2,
    ! sum bvel c^2 = 1/3.
    !
    ! The order this reaches is MEASURED rather than claimed: the tableau was
    ! constructed for the special (velocity-independent) Nystrom conditions, and
    ! the general case has strictly more conditions, so it is not assumed to
    ! retain order 4 here. The test reports the observed order.
    subroutine rkng_step(f2, t, y, yp, h, s, c, abar, avel, bbar, bvel, ks, ytmp, &
                         yptmp, ynew, ypnew, nfev, ctx)
        procedure(ode_rhs_rkng_t) :: f2
        real(dp), intent(in) :: t, h
        real(dp), intent(in) :: y(:), yp(:)
        integer,  intent(in) :: s
        real(dp), intent(in) :: c(:), abar(:, :), avel(:, :), bbar(:), bvel(:)
        real(dp), intent(inout) :: ks(:, :)
        real(dp), intent(inout) :: ytmp(:), yptmp(:)
        real(dp), intent(out) :: ynew(:), ypnew(:)
        integer, intent(inout) :: nfev
        class(*), intent(in), optional :: ctx

        real(dp) :: acc
        integer :: i, j, k, neq

        neq = size(y)

        do i = 1, s
            do k = 1, neq
                acc = 0.0_dp
                do j = 1, i - 1
                    acc = acc + abar(i, j)*ks(j, k)
                end do
                ytmp(k) = y(k) + c(i)*h*yp(k) + h*h*acc
                ! Velocity at the stage, needed because f depends on y'.
                acc = 0.0_dp
                do j = 1, i - 1
                    acc = acc + avel(i, j)*ks(j, k)
                end do
                yptmp(k) = yp(k) + h*acc
            end do
            call f2(t + c(i)*h, ytmp(1:neq), yptmp(1:neq), ks(i, 1:neq), ctx)
            nfev = nfev + 1
        end do

        do k = 1, neq
            acc = 0.0_dp
            do i = 1, s
                acc = acc + bbar(i)*ks(i, k)
            end do
            ynew(k) = y(k) + h*yp(k) + h*h*acc

            acc = 0.0_dp
            do i = 1, s
                acc = acc + bvel(i)*ks(i, k)
            end do
            ypnew(k) = yp(k) + h*acc
        end do
    end subroutine rkng_step

    ! Fixed-step RKNG driver for y'' = f(t, y, y').
    subroutine rkng_integrate_fixed(f2, t0, t1, y0, yp0, nsteps, yend, ypend, &
                                    nfev, status, ctx)
        procedure(ode_rhs_rkng_t) :: f2
        real(dp), intent(in) :: t0, t1
        real(dp), intent(in) :: y0(:), yp0(:)
        integer,  intent(in) :: nsteps
        real(dp), intent(out) :: yend(:), ypend(:)
        integer,  intent(out) :: nfev
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: ctx

        real(dp) :: c(TDRK_MAX_STAGES), abar(TDRK_MAX_STAGES, TDRK_MAX_STAGES)
        real(dp) :: bbar(TDRK_MAX_STAGES), bvel(TDRK_MAX_STAGES)
        real(dp) :: avel(TDRK_MAX_STAGES, TDRK_MAX_STAGES)
        real(dp), allocatable :: ks(:, :), ytmp(:), yptmp(:)
        real(dp), allocatable :: y(:), yp(:), ynew(:), ypnew(:)
        real(dp) :: h, t
        integer :: s, k, neq
        logical :: ok

        call status_set(status, FORTNUM_OK, "")
        if (nsteps < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "rkng_integrate_fixed: nsteps must be positive")
            return
        end if

        call tdrk_tableau(4, .true., s, c, abar, bbar, ok)
        if (.not. ok) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "rkng_integrate_fixed: tableau construction failed")
            return
        end if

        ! Velocity weights solve sum bvel = 1, sum bvel c = 1/2,
        ! sum bvel c^2 = 1/3 on c = (0, 1/3, 2/3).
        bvel = 0.0_dp
        bvel(1) = 0.25_dp
        bvel(2) = 0.0_dp
        bvel(3) = 0.75_dp

        ! Stage-velocity coefficients, row sums equal to c_i (consistency).
        avel = 0.0_dp
        avel(2, 1) = c(2)
        avel(3, 2) = c(3)

        neq = size(y0)
        allocate (ks(TDRK_MAX_STAGES, neq), ytmp(neq), yptmp(neq))
        allocate (y(neq), yp(neq), ynew(neq), ypnew(neq))
        ks = 0.0_dp
        y = y0
        yp = yp0
        h = (t1 - t0)/real(nsteps, dp)
        nfev = 0

        t = t0
        do k = 1, nsteps
            call rkng_step(f2, t, y, yp, h, s, c, abar, avel, bbar, bvel, ks, &
                           ytmp, yptmp, ynew, ypnew, nfev, ctx)
            y = ynew
            yp = ypnew
            t = t0 + real(k, dp)*h
        end do
        yend(1:neq) = y
        ypend(1:neq) = yp
    end subroutine rkng_integrate_fixed

    ! ------------------------------------------------------------------ error
    ! control
    !
    ! Both adaptive drivers below use the three-stage order-4 tableau, because
    ! the two-stage one cannot carry an embedded estimate: with c = (0, 1/2) the
    ! order-3 conditions (sum b = 1/2, sum b c = 1/6) already force
    ! b = (1/6, 1/3), which is the order-4 solution itself, so the difference of
    ! the two is identically zero. Three stages leave one free parameter.

    ! Weighted RMS norm of an error vector against the mixed tolerance, the
    ! standard criterion: accepted when the result is <= 1.
    pure function error_norm(err, y, ynew, rtol, atol) result(e)
        real(dp), intent(in) :: err(:), y(:), ynew(:), rtol, atol
        real(dp) :: e, sc
        integer  :: k

        e = 0.0_dp
        do k = 1, size(err)
            sc = atol + rtol*max(abs(y(k)), abs(ynew(k)))
            if (sc <= 0.0_dp) cycle
            e = e + (err(k)/sc)**2
        end do
        e = sqrt(e/real(max(size(err), 1), dp))
    end function error_norm

    ! Elementary step-size controller. phat is the order of the EMBEDDED
    ! solution, so the local error it measures is O(h^(phat+1)).
    pure function step_factor(e, phat, accepted) result(fac)
        real(dp), intent(in) :: e
        integer,  intent(in) :: phat
        logical,  intent(in) :: accepted
        real(dp) :: fac

        if (e <= tiny(1.0_dp)) then
            fac = FAC_MAX
        else
            fac = SAFETY*e**(-1.0_dp/real(phat + 1, dp))
        end if
        fac = min(FAC_MAX, max(FAC_MIN, fac))
        if (.not. accepted) fac = min(fac, 1.0_dp)
    end function step_factor

    ! Adaptive TDRK for zdot = F(z) with G = F'F.
    !
    ! Propagated solution: three-stage order 4, bbar = (1/8, 1/4, 1/8).
    ! Embedded solution:   order 3, bhat = (0, 1/2, 0) -- it satisfies
    ! sum bhat = 1/2 and sum bhat c = 1/6 but gives sum bhat c^2 = 1/18 rather
    ! than the 1/12 order 4 needs, so it is genuinely one order lower and the
    ! difference is a real error estimate rather than round-off.
    !
    ! h0 <= 0 requests a default first step of (t1 - t0)/100; the controller
    ! recovers from a bad guess within a few steps either way.
    subroutine tdrk_integrate_adaptive(rhs, rhs2, t0, t1, y0, rtol, atol, h0, &
                                       yend, hlast, nfev_f, nfev_g, naccept, &
                                       nreject, status, ctx)
        procedure(ode_rhs_t)  :: rhs
        procedure(ode_rhs2_t) :: rhs2
        real(dp), intent(in) :: t0, t1
        real(dp), intent(in) :: y0(:)
        real(dp), intent(in) :: rtol, atol, h0
        real(dp), intent(out) :: yend(:)
        real(dp), intent(out) :: hlast
        integer,  intent(out) :: nfev_f, nfev_g, naccept, nreject
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: ctx

        real(dp) :: c(TDRK_MAX_STAGES), abar(TDRK_MAX_STAGES, TDRK_MAX_STAGES)
        real(dp) :: bbar(TDRK_MAX_STAGES), bhat(TDRK_MAX_STAGES)
        real(dp), allocatable :: gs(:, :), ytmp(:), y(:), ynew(:), err(:)
        real(dp) :: h, t, span, e
        integer  :: s, k, i, neq, nstep
        logical  :: ok, accepted

        call status_set(status, FORTNUM_OK, "")
        nfev_f = 0
        nfev_g = 0
        naccept = 0
        nreject = 0
        neq = size(y0)
        yend(1:neq) = y0
        span = t1 - t0
        hlast = h0
        if (span == 0.0_dp) return

        if (rtol <= 0.0_dp .and. atol <= 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "tdrk_integrate_adaptive: rtol and atol both non-positive")
            return
        end if

        call tdrk_tableau(4, .true., s, c, abar, bbar, ok)
        if (.not. ok) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "tdrk_integrate_adaptive: tableau construction failed")
            return
        end if
        bhat = 0.0_dp
        bhat(2) = 0.5_dp

        allocate (gs(TDRK_MAX_STAGES, neq), ytmp(neq), y(neq), ynew(neq), err(neq))
        gs = 0.0_dp
        y = y0
        t = t0

        if (h0 > 0.0_dp) then
            h = sign(min(h0, abs(span)), span)
        else
            h = span/100.0_dp
        end if

        do nstep = 1, MAX_STEPS
            if (abs(t - t0) >= abs(span)) exit
            if (abs(h) > abs(t1 - t)) h = t1 - t

            call tdrk_step(rhs, rhs2, t, y, h, s, c, abar, bbar, gs, ytmp, ynew, &
                           nfev_f, nfev_g, ctx)
            do k = 1, neq
                err(k) = 0.0_dp
                do i = 1, s
                    err(k) = err(k) + (bbar(i) - bhat(i))*gs(i, k)
                end do
                err(k) = h*h*err(k)
            end do

            e = error_norm(err, y, ynew, rtol, atol)
            accepted = e <= 1.0_dp
            if (accepted) then
                t = t + h
                y = ynew
                naccept = naccept + 1
            else
                nreject = nreject + 1
            end if
            hlast = h
            h = h*step_factor(e, 3, accepted)

            if (abs(h) < 16.0_dp*spacing(abs(t))) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                                "tdrk_integrate_adaptive: step size underflow")
                yend(1:neq) = y
                return
            end if
        end do

        if (abs(t - t0) < abs(span)) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                            "tdrk_integrate_adaptive: step limit reached")
        end if
        yend(1:neq) = y
    end subroutine tdrk_integrate_adaptive

    ! Adaptive RKNG for y'' = f(t, y, y').
    !
    ! Unlike TDRK the velocity is an independent state variable here, so the
    ! error estimate covers both: the position embedded weights satisfy only
    ! sum bhat = 1/2 and the velocity ones only sum bhat_v = 1 and
    ! sum bhat_v c = 1/2, one order below the propagated pair in each case. The
    ! two contributions go into one norm, so neither can be met at the other's
    ! expense.
    !
    ! The propagated method's order on a general (velocity-dependent) f is
    ! measured rather than claimed -- the tableau was built for the special
    ! Nystrom conditions -- so the controller exponent is set from the embedded
    ! order 2, which is safe regardless of which order the propagated pair
    ! actually reaches.
    subroutine rkng_integrate_adaptive(f2, t0, t1, y0, yp0, rtol, atol, h0, &
                                       yend, ypend, hlast, nfev, naccept, &
                                       nreject, status, ctx)
        procedure(ode_rhs_rkng_t) :: f2
        real(dp), intent(in) :: t0, t1
        real(dp), intent(in) :: y0(:), yp0(:)
        real(dp), intent(in) :: rtol, atol, h0
        real(dp), intent(out) :: yend(:), ypend(:)
        real(dp), intent(out) :: hlast
        integer,  intent(out) :: nfev, naccept, nreject
        type(fortnum_status_t), intent(out) :: status
        class(*), intent(in), optional :: ctx

        real(dp) :: c(TDRK_MAX_STAGES), abar(TDRK_MAX_STAGES, TDRK_MAX_STAGES)
        real(dp) :: avel(TDRK_MAX_STAGES, TDRK_MAX_STAGES)
        real(dp) :: bbar(TDRK_MAX_STAGES), bvel(TDRK_MAX_STAGES)
        real(dp) :: bhat(TDRK_MAX_STAGES), bhat_vel(TDRK_MAX_STAGES)
        real(dp), allocatable :: ks(:, :), ytmp(:), yptmp(:)
        real(dp), allocatable :: y(:), yp(:), ynew(:), ypnew(:)
        real(dp), allocatable :: err(:), errv(:)
        real(dp) :: h, t, span, e
        integer  :: s, k, i, neq, nstep
        logical  :: ok, accepted

        call status_set(status, FORTNUM_OK, "")
        nfev = 0
        naccept = 0
        nreject = 0
        neq = size(y0)
        yend(1:neq) = y0
        ypend(1:neq) = yp0
        span = t1 - t0
        hlast = h0
        if (span == 0.0_dp) return

        if (rtol <= 0.0_dp .and. atol <= 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "rkng_integrate_adaptive: rtol and atol both non-positive")
            return
        end if

        call tdrk_tableau(4, .true., s, c, abar, bbar, ok)
        if (.not. ok) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "rkng_integrate_adaptive: tableau construction failed")
            return
        end if

        bvel = 0.0_dp
        bvel(1) = 0.25_dp
        bvel(3) = 0.75_dp
        avel = 0.0_dp
        avel(2, 1) = c(2)
        avel(3, 2) = c(3)

        ! Embedded position: sum bhat = 1/2 only.
        bhat = 0.0_dp
        bhat(1) = 0.5_dp
        ! Embedded velocity: sum = 1, sum c = 1/2, but sum c^2 = 1/6 /= 1/3.
        bhat_vel = 0.0_dp
        bhat_vel(1) = -0.5_dp
        bhat_vel(2) = 1.5_dp

        allocate (ks(TDRK_MAX_STAGES, neq), ytmp(neq), yptmp(neq))
        allocate (y(neq), yp(neq), ynew(neq), ypnew(neq), err(neq), errv(neq))
        ks = 0.0_dp
        y = y0
        yp = yp0
        t = t0

        if (h0 > 0.0_dp) then
            h = sign(min(h0, abs(span)), span)
        else
            h = span/100.0_dp
        end if

        do nstep = 1, MAX_STEPS
            if (abs(t - t0) >= abs(span)) exit
            if (abs(h) > abs(t1 - t)) h = t1 - t

            call rkng_step(f2, t, y, yp, h, s, c, abar, avel, bbar, bvel, ks, &
                           ytmp, yptmp, ynew, ypnew, nfev, ctx)
            do k = 1, neq
                err(k) = 0.0_dp
                errv(k) = 0.0_dp
                do i = 1, s
                    err(k) = err(k) + (bbar(i) - bhat(i))*ks(i, k)
                    errv(k) = errv(k) + (bvel(i) - bhat_vel(i))*ks(i, k)
                end do
                err(k) = h*h*err(k)
                errv(k) = h*errv(k)
            end do

            e = max(error_norm(err, y, ynew, rtol, atol), &
                    error_norm(errv, yp, ypnew, rtol, atol))
            accepted = e <= 1.0_dp
            if (accepted) then
                t = t + h
                y = ynew
                yp = ypnew
                naccept = naccept + 1
            else
                nreject = nreject + 1
            end if
            hlast = h
            h = h*step_factor(e, 2, accepted)

            if (abs(h) < 16.0_dp*spacing(abs(t))) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                                "rkng_integrate_adaptive: step size underflow")
                yend(1:neq) = y
                ypend(1:neq) = yp
                return
            end if
        end do

        if (abs(t - t0) < abs(span)) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                            "rkng_integrate_adaptive: step limit reached")
        end if
        yend(1:neq) = y
        ypend(1:neq) = yp
    end subroutine rkng_integrate_adaptive

end module fortnum_ode_tdrk
