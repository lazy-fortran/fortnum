module fortnum_ode_gauss_radau
    ! Adaptive 15th-order Gauss-Radau integrator for first-order systems, in the
    ! spirit of IAS15 (Rein & Spiegel, MNRAS 446 (2015) 1424, arXiv:1409.4779)
    ! and Everhart's RADAU.
    !
    ! Difference from REBOUND's IAS15: that one is formulated for the
    ! second-order system y'' = f(t, y, y') of gravitational dynamics. This one
    ! is formulated for the FIRST-order system y' = F(t, y), because the systems
    ! it is meant for here -- guiding-centre and full-orbit charged particle
    ! motion -- are first order and non-separable, and cannot be put in the
    ! second-order form REBOUND assumes.
    !
    ! Method. Over one step t = t0 + h*s, s in [0,1], the derivative is
    ! represented as a degree-7 polynomial through the eight Gauss-Radau nodes,
    ! held in the Newton (divided-difference) basis:
    !
    !     F(s) = sum_{k=0..7} g_k N_k(s),   N_k(s) = prod_{m=1..k} (s - c_m)
    !
    ! which integrates termwise to give y(s). The g coefficients are found by
    ! predictor-corrector iteration: evaluate F at the nodes using the current
    ! g, recompute the divided differences, repeat until g stops moving. With
    ! eight Radau nodes the collocation is 15th order.
    !
    ! Node provenance. The nodes are COMPUTED, not transcribed. The free nodes
    ! of a left-Radau rule on [-1,1] are the roots of (P_{n-1}(x) + P_n(x))/(1+x)
    ! for Legendre P; they are found here by bracketing and bisection and then
    ! mapped to [0,1]. Transcribing 25-digit constants from a paper is a silent
    ! failure mode, and a wrong node would still produce a plausible-looking
    ! integrator of reduced order. The order-of-convergence test is what pins
    ! them down.
    !
    ! Derivative policy: trace_rule (ad.md sections 1 and 4), as for the other
    ! steppers here -- the adaptive schedule is the primal's, and a derivative
    ! differentiates the frozen mesh.
    !
    ! No module-level state. The driver allocates its work arrays once, outside
    ! the step loop, and hands them to the stepper.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
                              FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_ode, only: ode_rhs_t, ode_problem_t, ode_solution_t

    implicit none
    private

    public :: ode_integrate_radau, ode_solve_radau
    public :: radau_nodes, radau_newton_basis

    ! Eight nodes, so seven free b coefficients and formal order 15.
    integer, parameter, public :: RADAU_STAGES = 8
    integer, parameter, public :: RADAU_NB = RADAU_STAGES - 1

    ! Predictor-corrector control. IAS15 iterates until the last b coefficient
    ! stops changing at round-off; the cap is a guard, not the usual exit.
    integer,  parameter :: MAX_PC_ITER = 12
    real(dp), parameter :: PC_TOL = 1.0e-16_dp

    ! Step-size control. The local error is estimated from the neglected last
    ! Newton term, h * g_7 * integral(N_7), which is O(h^8) in y, so the
    ! controller exponent is 1/8. Note this embedded estimate behaves like an
    ! order-7 method even though the collocation itself is order 15; that is
    ! deliberate and conservative, and it is what makes the step honest in
    ! y-units rather than a dimensionless IAS15-style epsilon.
    real(dp), parameter :: SAFETY = 0.85_dp
    real(dp), parameter :: FAC_MIN = 0.2_dp
    real(dp), parameter :: FAC_MAX = 8.0_dp
    real(dp), parameter :: ERR_EXPONENT = 1.0_dp/8.0_dp

contains

    ! Legendre P_n(x) by the standard three-term recurrence.
    pure function legendre(n, x) result(p)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: p

        real(dp) :: pm1, pm2
        integer :: k

        if (n == 0) then
            p = 1.0_dp
            return
        end if
        pm2 = 1.0_dp
        pm1 = x
        if (n == 1) then
            p = pm1
            return
        end if
        do k = 2, n
            p = (real(2*k - 1, dp)*x*pm1 - real(k - 1, dp)*pm2)/real(k, dp)
            pm2 = pm1
            pm1 = p
        end do
    end function legendre

    ! Radau defining function, with the known root at x = -1 divided out.
    pure function radau_poly(n, x) result(v)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: v

        v = (legendre(n - 1, x) + legendre(n, x))/(1.0_dp + x)
    end function radau_poly

    ! Left-Radau nodes on [0,1]: c(1) = 0 fixed, the remaining n-1 computed as
    ! roots of the Radau polynomial and mapped from [-1,1].
    subroutine radau_nodes(n, c)
        integer, intent(in) :: n
        real(dp), intent(out) :: c(:)

        integer, parameter :: NSCAN = 20000
        real(dp) :: xa, xb, fa, fb, xm, fm
        integer :: i, iroot, it
        real(dp) :: lo, hi

        c = 0.0_dp
        c(1) = 0.0_dp

        iroot = 1
        lo = -1.0_dp + 1.0e-10_dp
        hi = 1.0_dp
        xa = lo
        fa = radau_poly(n, xa)
        do i = 1, NSCAN
            xb = lo + (hi - lo)*real(i, dp)/real(NSCAN, dp)
            fb = radau_poly(n, xb)
            if (fa == 0.0_dp) then
                iroot = iroot + 1
                c(iroot) = 0.5_dp*(xa + 1.0_dp)
            else if (fa*fb < 0.0_dp) then
                ! Bisection: robust, and the polynomial is cheap.
                do it = 1, 200
                    xm = 0.5_dp*(xa + xb)
                    fm = radau_poly(n, xm)
                    if (fm == 0.0_dp) exit
                    if (fa*fm < 0.0_dp) then
                        xb = xm
                        fb = fm
                    else
                        xa = xm
                        fa = fm
                    end if
                end do
                xm = 0.5_dp*(xa + xb)
                iroot = iroot + 1
                if (iroot <= n) c(iroot) = 0.5_dp*(xm + 1.0_dp)
                xa = lo + (hi - lo)*real(i, dp)/real(NSCAN, dp)
                fa = radau_poly(n, xa)
                cycle
            end if
            xa = xb
            fa = fb
        end do
    end subroutine radau_nodes

    ! Power-basis coefficients of the Newton basis polynomials
    !     N_k(s) = prod_{m=1..k} (s - c_m),   N_0 = 1
    ! built by incremental convolution. Only multiplications and subtractions of
    ! O(1) quantities, so this is stable.
    !
    ! Why not a Vandermonde solve. Fitting F(s) = F0 + sum_j b_j s^j directly by
    ! inverting M(i,j) = c(i+1)**j is algebraically equivalent but numerically
    ! hopeless at tight tolerance: b_7 is recovered from O(h) inputs through
    ! O(h^7) cancellation, so round-off puts a floor under the error indicator,
    ! the controller can never meet the tolerance, and h collapses. Everhart's
    ! divided-difference form avoids the inversion entirely.
    subroutine radau_newton_basis(c, a)
        real(dp), intent(in) :: c(:)
        real(dp), intent(out) :: a(0:, 0:)   ! a(k, p) = coeff of s^p in N_k

        integer :: k, p

        a = 0.0_dp
        a(0, 0) = 1.0_dp
        do k = 1, RADAU_NB
            ! N_k = N_{k-1} * (s - c_k)
            do p = k, 1, -1
                a(k, p) = a(k - 1, p - 1)
            end do
            a(k, 0) = 0.0_dp
            do p = 0, k - 1
                a(k, p) = a(k, p) - c(k)*a(k - 1, p)
            end do
        end do
    end subroutine radau_newton_basis

    ! Newton divided differences of the stage derivatives over the nodes.
    ! g(k, :) is the k-th divided difference, so
    !     F(s) = sum_{k=0..7} g(k,:) N_k(s).
    subroutine radau_divided_differences(c, fs, g, neq)
        real(dp), intent(in) :: c(:)
        real(dp), intent(in) :: fs(0:, :)     ! (0:RADAU_NB, neq) values at nodes
        real(dp), intent(out) :: g(0:, :)     ! (0:RADAU_NB, neq)
        integer, intent(in) :: neq

        integer :: level, k

        g(0:RADAU_NB, 1:neq) = fs(0:RADAU_NB, 1:neq)
        do level = 1, RADAU_NB
            do k = RADAU_NB, level, -1
                g(k, 1:neq) = (g(k, 1:neq) - g(k - 1, 1:neq)) / &
                              (c(k + 1) - c(k + 1 - level))
            end do
        end do
    end subroutine radau_divided_differences

    ! Evaluate y(s) = y0 + h * integral_0^s F(u) du from the Newton coefficients.
    pure subroutine radau_state_at(y0, h, s, g, a, neq, yout)
        real(dp), intent(in) :: y0(:), h, s
        real(dp), intent(in) :: g(0:, :), a(0:, 0:)
        integer, intent(in) :: neq
        real(dp), intent(out) :: yout(:)

        integer :: k, p, m
        real(dp) :: spow(0:RADAU_NB + 1), integ, acc

        spow(0) = 1.0_dp
        do p = 1, RADAU_NB + 1
            spow(p) = spow(p - 1)*s
        end do

        do m = 1, neq
            acc = 0.0_dp
            do k = 0, RADAU_NB
                integ = 0.0_dp
                do p = 0, k
                    integ = integ + a(k, p)*spow(p + 1)/real(p + 1, dp)
                end do
                acc = acc + g(k, m)*integ
            end do
            yout(m) = y0(m) + h*acc
        end do
    end subroutine radau_state_at

    ! One Gauss-Radau step from (t, y) of size h.
    !
    ! g is intent(inout): the previous step's coefficients, rescaled by the step
    ! ratio, predict this step's and cut the corrector to two or three
    ! iterations. g_k carries a factor h^k, hence the (h/h_prev)**k rescaling.
    subroutine gauss_radau_step(rhs, t, y, h, h_prev, c, a, i_last, rtol, atol, &
                                g, fs, ytmp, ynew, err, nfev, niter, ctx)
        procedure(ode_rhs_t) :: rhs
        real(dp), intent(in)    :: t, h, h_prev
        real(dp), intent(in)    :: y(:)
        real(dp), intent(in)    :: c(:), a(0:, 0:)
        real(dp), intent(in)    :: i_last, rtol, atol
        real(dp), intent(inout) :: g(0:, :)
        real(dp), intent(inout) :: fs(0:, :)
        real(dp), intent(out)   :: ytmp(:), ynew(:)
        real(dp), intent(out)   :: err
        integer,  intent(inout) :: nfev
        integer,  intent(out)   :: niter
        class(*), intent(in), optional :: ctx

        integer :: i, k, neq, it
        real(dp) :: ratio, rk, scale_f, glast_prev, delta, sc, e

        neq = size(y)

        if (h_prev /= 0.0_dp) then
            ratio = h/h_prev
            rk = 1.0_dp
            do k = 1, RADAU_NB
                rk = rk*ratio
                g(k, 1:neq) = g(k, 1:neq)*rk
            end do
        else
            g = 0.0_dp
        end if

        call rhs(t, y, fs(0, 1:neq), ctx)
        nfev = nfev + 1
        g(0, 1:neq) = fs(0, 1:neq)

        niter = 0
        do it = 1, MAX_PC_ITER
            niter = it
            glast_prev = maxval(abs(g(RADAU_NB, 1:neq)))

            do i = 1, RADAU_NB
                call radau_state_at(y, h, c(i + 1), g, a, neq, ytmp)
                call rhs(t + h*c(i + 1), ytmp(1:neq), fs(i, 1:neq), ctx)
                nfev = nfev + 1
            end do

            call radau_divided_differences(c, fs, g, neq)

            scale_f = maxval(abs(fs(0, 1:neq)))
            if (scale_f <= 0.0_dp) scale_f = 1.0_dp
            delta = abs(maxval(abs(g(RADAU_NB, 1:neq))) - glast_prev)/scale_f
            if (delta <= PC_TOL) exit
        end do

        call radau_state_at(y, h, 1.0_dp, g, a, neq, ynew)

        ! Local error estimate in y-units: the contribution of the last Newton
        ! term to the step, h * g_7 * integral_0^1 N_7(s) ds, measured against
        ! the usual mixed absolute/relative scale.
        err = 0.0_dp
        do k = 1, neq
            sc = atol + rtol*max(abs(y(k)), abs(ynew(k)))
            if (sc <= 0.0_dp) sc = atol
            if (sc <= 0.0_dp) sc = 1.0_dp
            e = abs(h*g(RADAU_NB, k)*i_last)/sc
            if (e > err) err = e
        end do
    end subroutine gauss_radau_step

    ! Adaptive driver over [problem%t0, problem%t1], recording the trace in
    ! solution the same way ode_integrate and ode_integrate_dop do.
    subroutine ode_integrate_radau(problem, solution, status)
        type(ode_problem_t), intent(in) :: problem
        type(ode_solution_t), intent(inout) :: solution
        type(fortnum_status_t), intent(out) :: status

        real(dp), allocatable :: c(:), a(:, :), g(:, :), fs(:, :)
        real(dp), allocatable :: ytmp(:), ynew(:), yv(:)
        real(dp), allocatable :: t_rec(:), y_rec(:, :), h_rec(:), e_rec(:)
        real(dp) :: t, h, h_prev, hdir, err, fac, span, i_last
        integer  :: neq, nsteps, nrej, nfev, niter, cap, k, rec_cap

        call status_set(status, FORTNUM_OK, "")

        if (.not. allocated(problem%y0)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "ode_integrate_radau: y0 not allocated")
            return
        end if
        if (.not. associated(problem%rhs)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "ode_integrate_radau: rhs not associated")
            return
        end if
        if (problem%rtol <= 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "ode_integrate_radau: rtol must be positive")
            return
        end if

        neq = size(problem%y0)
        span = problem%t1 - problem%t0
        if (span == 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                            "ode_integrate_radau: t1 must differ from t0")
            return
        end if
        hdir = sign(1.0_dp, span)

        allocate (c(RADAU_STAGES), a(0:RADAU_NB, 0:RADAU_NB))
        call radau_nodes(RADAU_STAGES, c)
        call radau_newton_basis(c, a)

        i_last = 0.0_dp
        do k = 0, RADAU_NB
            i_last = i_last + a(RADAU_NB, k)/real(k + 1, dp)
        end do

        allocate (g(0:RADAU_NB, neq), fs(0:RADAU_NB, neq))
        allocate (ytmp(neq), ynew(neq), yv(neq))
        g = 0.0_dp
        h_prev = 0.0_dp

        ! Trace capacity grows by doubling rather than being allocated at
        ! max_steps up front. Callers that drive this one macro-step at a time
        ! (SIMPLE's orbit loop does) would otherwise pay a max_steps-sized
        ! allocation per call, which would dominate the very wall-clock the
        ! integrator benchmark is trying to measure.
        cap = max(problem%max_steps, 1)
        rec_cap = min(cap + 1, 256)
        allocate (t_rec(rec_cap), y_rec(neq, rec_cap), h_rec(rec_cap), e_rec(rec_cap))

        t = problem%t0
        yv = problem%y0
        h = problem%h0
        if (h == 0.0_dp) h = 1.0e-3_dp*abs(span)
        h = hdir*min(abs(h), abs(span))
        if (problem%hmax > 0.0_dp) h = hdir*min(abs(h), problem%hmax)

        nsteps = 0
        nrej = 0
        nfev = 0
        t_rec(1) = t
        y_rec(:, 1) = yv
        h_rec(1) = 0.0_dp
        e_rec(1) = 0.0_dp

        do while ((problem%t1 - t)*hdir > 0.0_dp)
            if (nsteps >= cap) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                                "ode_integrate_radau: step limit reached")
                exit
            end if
            if (abs(h) > abs(problem%t1 - t)) h = problem%t1 - t

            call gauss_radau_step(problem%rhs, t, yv, h, h_prev, c, a, i_last, &
                                  problem%rtol, problem%atol, g, fs, &
                                  ytmp, ynew, err, nfev, niter)

            if (err <= 1.0_dp .or. abs(h) <= abs(problem%hmin)) then
                t = t + h
                yv = ynew
                nsteps = nsteps + 1
                if (nsteps + 1 > rec_cap) then
                    call grow_trace(t_rec, y_rec, h_rec, e_rec, rec_cap, cap + 1)
                end if
                t_rec(nsteps + 1) = t
                y_rec(:, nsteps + 1) = yv
                h_rec(nsteps + 1) = h
                e_rec(nsteps + 1) = err
                h_prev = h
            else
                nrej = nrej + 1
                ! Predictor is no longer trustworthy after a reject.
                g = 0.0_dp
                h_prev = 0.0_dp
            end if

            ! Step-size update, shared by the accept and reject paths.
            if (err > 0.0_dp) then
                fac = SAFETY*(1.0_dp/err)**ERR_EXPONENT
            else
                fac = FAC_MAX
            end if
            fac = max(FAC_MIN, min(FAC_MAX, fac))
            h = h*fac
            if (problem%hmax > 0.0_dp) h = hdir*min(abs(h), problem%hmax)
            if (problem%hmin > 0.0_dp) h = hdir*max(abs(h), problem%hmin)
        end do

        solution%nsteps = nsteps
        solution%nrejected = nrej
        solution%nfev = nfev
        if (allocated(solution%t)) deallocate (solution%t)
        if (allocated(solution%y)) deallocate (solution%y)
        if (allocated(solution%h)) deallocate (solution%h)
        if (allocated(solution%err)) deallocate (solution%err)
        allocate (solution%t(nsteps + 1), solution%y(neq, nsteps + 1))
        allocate (solution%h(nsteps + 1), solution%err(nsteps + 1))
        solution%t = t_rec(1:nsteps + 1)
        solution%y = y_rec(:, 1:nsteps + 1)
        solution%h = h_rec(1:nsteps + 1)
        solution%err = e_rec(1:nsteps + 1)
        solution%status = status
    end subroutine ode_integrate_radau

    ! Flat convenience form, matching ode_solve_dop: the recorded trace comes
    ! back in allocatable arrays rather than just the endpoint.
    subroutine ode_solve_radau(rhs, t0, t1, y0, t_out, y_out, status, rtol, atol)
        procedure(ode_rhs_t) :: rhs
        real(dp), intent(in) :: t0, t1
        real(dp), intent(in) :: y0(:)
        real(dp), allocatable, intent(out) :: t_out(:)
        real(dp), allocatable, intent(out) :: y_out(:, :)
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional :: rtol, atol

        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution

        problem%rhs => rhs
        problem%t0 = t0
        problem%t1 = t1
        problem%y0 = y0
        if (present(rtol)) problem%rtol = rtol
        if (present(atol)) problem%atol = atol

        call ode_integrate_radau(problem, solution, status)

        if (allocated(solution%t)) then
            t_out = solution%t
            y_out = solution%y
        else
            allocate (t_out(1), y_out(size(y0), 1))
            t_out(1) = t0
            y_out(:, 1) = y0
        end if
    end subroutine ode_solve_radau


    ! Double the recorded-trace capacity, never exceeding hard_cap.
    subroutine grow_trace(t_rec, y_rec, h_rec, e_rec, rec_cap, hard_cap)
        real(dp), allocatable, intent(inout) :: t_rec(:), y_rec(:, :)
        real(dp), allocatable, intent(inout) :: h_rec(:), e_rec(:)
        integer, intent(inout) :: rec_cap
        integer, intent(in) :: hard_cap

        real(dp), allocatable :: tt(:), yy(:, :), hh(:), ee(:)
        integer :: newcap, neq

        newcap = min(max(2*rec_cap, rec_cap + 1), hard_cap)
        if (newcap <= rec_cap) return
        neq = size(y_rec, 1)

        allocate (tt(newcap), yy(neq, newcap), hh(newcap), ee(newcap))
        tt(1:rec_cap) = t_rec(1:rec_cap)
        yy(:, 1:rec_cap) = y_rec(:, 1:rec_cap)
        hh(1:rec_cap) = h_rec(1:rec_cap)
        ee(1:rec_cap) = e_rec(1:rec_cap)
        call move_alloc(tt, t_rec)
        call move_alloc(yy, y_rec)
        call move_alloc(hh, h_rec)
        call move_alloc(ee, e_rec)
        rec_cap = newcap
    end subroutine grow_trace

end module fortnum_ode_gauss_radau
