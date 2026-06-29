module fortnum_multiroot_rc
    ! Reverse-communication (RC) mirror of multiroot_hybrid (fortnum_multiroot):
    ! the same Powell hybrid Newton step with a backtracking line search on the
    ! residual 1/2 |F|^2, but with control inverted so it runs on a device.
    !
    ! The host multiroot_hybrid takes a procedure(multiroot_fdf_t) callback and
    ! a class(*) ctx, and the line search re-enters that callback.  Neither
    ! procedure-pointer dispatch nor polymorphic dummies are allowed inside an
    ! !$acc routine seq.  multiroot_step inverts control: the caller drives the
    ! loop, evaluates F and the analytic Jacobian inline at state%x, and calls
    ! multiroot_step repeatedly.  multiroot_step is a pure leaf !$acc routine
    ! seq that resumes from multiroot_rc_t, runs one logical action of the SAME
    ! dogleg / line-search / LU iteration, and returns an action token telling
    ! the caller what to do next (evaluate F and J, stop, or fail).
    !
    ! This is additive: the host multiroot_hybrid / multiroot_hybrids callback
    ! API is unchanged.  The inner LU is shared single-source via
    ! fortnum_linalg%lu_solve; only the phase state machine (genuinely different
    ! control flow from the host do-loop) lives here.
    !
    ! Constraints for routine seq (all met by construction): fixed compile-time
    ! MULTIROOT_RC_MAX_N so every array is a fixed-size component on the stack;
    ! NO allocatable, NO procedure pointers, NO class(*), NO I/O.  Status is a
    ! bare integer action plus an integer fail_code; the host reconstructs the
    ! fortnum_status_t message after the loop (string work is forbidden here).
    !
    ! DERIVATIVE POLICY (ad.md S1, S4): the RC solve is primal_only -- the line
    ! search, the LU steps and the phase machine are iteration, not a
    ! differentiable map.  The differentiable surface is the residual the caller
    ! supplies plus the implicit-rule sensitivity dx*/dp = -J_x^{-1} J_p at the
    ! converged root, identical to multiroot_hybrid's implicit_rule policy.  Do
    ! not differentiate through multiroot_step.
    use fortnum_kinds, only: dp
    use fortnum_linalg, only: lu_solve, LINALG_OK
    implicit none
    private

    ! Largest system the RC stepper accepts (fixed device-stack arrays).
    integer, parameter, public :: MULTIROOT_RC_MAX_N = 8

    ! Action tokens returned by multiroot_step.
    integer, parameter, public :: MULTIROOT_NEED_FJ = 0 ! evaluate F,J at state%x
    integer, parameter, public :: MULTIROOT_DONE = 1 ! converged (root in x)
    integer, parameter, public :: MULTIROOT_FAILED = 2 ! see state%fail_code

    ! fail_code values surfaced on MULTIROOT_FAILED (host maps to fortnum codes).
    integer, parameter, public :: MULTIROOT_RC_SINGULAR = 1 ! singular Jacobian
    integer, parameter, public :: MULTIROOT_RC_MAXITER = 2 ! max_iter reached

    ! Internal phase of the resumable stepper.
    integer, parameter :: PHASE_START = 0 ! ingest F at x0
    integer, parameter :: PHASE_NEWTON = 1 ! have F,J at current x: solve, step
    integer, parameter :: PHASE_LINESEARCH = 2 ! have F,J at a trial point
    integer, parameter :: PHASE_DONE = 3 ! terminal

    ! Caller-owned RC state.  Flat: only scalars and fixed-size arrays, so it
    ! copies cleanly to the device and is valid as a routine-seq dummy.
    type, public :: multiroot_rc_t
        integer :: n = 0 ! system size (1..MULTIROOT_RC_MAX_N)
        integer :: phase = PHASE_START
        integer :: iter = 0 ! Newton steps taken
        integer :: max_iter = 1000
        integer :: ls_count = 0 ! line-search halvings tried
        integer :: fail_code = 0
        real(dp) :: xtol = 1.0e-10_dp
        real(dp) :: ftol = 1.0e-10_dp
        real(dp) :: lambda = 1.0_dp ! current line-search step length
        real(dp) :: g0 = 0.0_dp ! 1/2 |F|^2 at the accepted point
        real(dp) :: x(MULTIROOT_RC_MAX_N) = 0.0_dp ! current iterate / trial point
        real(dp) :: x_base(MULTIROOT_RC_MAX_N) = 0.0_dp ! accepted iterate
        real(dp) :: dx(MULTIROOT_RC_MAX_N) = 0.0_dp ! Newton direction
    end type multiroot_rc_t

    public :: multiroot_rc_init, multiroot_step

contains

    ! Initialise RC state for an n-dimensional system starting at x0.  Optional
    ! tolerances and iteration cap mirror multiroot_hybrid's defaults.  After
    ! this the caller evaluates F (and J) at state%x and calls multiroot_step.
    pure subroutine multiroot_rc_init(state, n, x0, xtol, ftol, max_iter)
        !$acc routine seq
        type(multiroot_rc_t), intent(out) :: state
        integer, intent(in) :: n
        real(dp), intent(in) :: x0(n)
        real(dp), intent(in), optional :: xtol, ftol
        integer, intent(in), optional :: max_iter
        state%n = n
        state%phase = PHASE_START
        state%iter = 0
        state%ls_count = 0
        state%fail_code = 0
        state%lambda = 1.0_dp
        state%g0 = 0.0_dp
        state%xtol = 1.0e-10_dp
        state%ftol = 1.0e-10_dp
        state%max_iter = 1000
        if (present(xtol)) state%xtol = xtol
        if (present(ftol)) state%ftol = ftol
        if (present(max_iter)) state%max_iter = max_iter
        state%x = 0.0_dp
        state%x_base = 0.0_dp
        state%dx = 0.0_dp
        state%x(1:n) = x0
        state%x_base(1:n) = x0
    end subroutine multiroot_rc_init

    ! One logical action of the Powell dogleg iteration, resumable from state.
    !
    ! The caller must have written f(1:n) = F(state%x) and jac(1:n,1:n) =
    ! dF/dx at state%x before each call.  On return, action is:
    !   MULTIROOT_NEED_FJ  evaluate F,J at the (possibly updated) state%x and
    !                      call again,
    !   MULTIROOT_DONE     converged; the root is in state%x(1:n),
    !   MULTIROOT_FAILED   state%fail_code holds MULTIROOT_RC_SINGULAR or
    !                      MULTIROOT_RC_MAXITER.
    !
    ! The math reproduces multiroot_hybrid exactly: stop on |F|_inf <= ftol or
    ! |dx|_inf <= xtol*(|x|_inf + xtol); line search from lambda = 1 with up to
    ! 30 halvings, Armijo slack (1 - 1e-4*lambda)*g0, and accept-last-trial
    ! fallback after 30 halvings.
    pure subroutine multiroot_step(state, f, jac, action)
        !$acc routine seq
        type(multiroot_rc_t), intent(inout) :: state
        real(dp), intent(in) :: f(state%n)
        real(dp), intent(in) :: jac(state%n, state%n)
        integer, intent(out) :: action

        real(dp) :: amat(MULTIROOT_RC_MAX_N, MULTIROOT_RC_MAX_N)
        real(dp) :: rhs(MULTIROOT_RC_MAX_N)
        real(dp) :: gtrial, fmax, dxmax, xmax
        integer :: n

        n = state%n
        fmax = inf_norm(f, n)

        select case (state%phase)
        case (PHASE_START)
            ! F evaluated at x0.  Immediate convergence check.
            if (fmax <= state%ftol) then
                state%phase = PHASE_DONE
                action = MULTIROOT_DONE
                return
            end if
            call begin_newton(state, f, jac, amat, rhs, action)
            return

        case (PHASE_NEWTON)
            call begin_newton(state, f, jac, amat, rhs, action)
            return

        case (PHASE_LINESEARCH)
            ! f,jac are at the trial point state%x = x_base + lambda*dx.
            gtrial = 0.5_dp*dot_n(f, f, n)
            if (gtrial < (1.0_dp - 1.0e-4_dp*state%lambda)*state%g0 &
                .or. state%ls_count >= 30) then
                ! Accept this trial (Armijo, or the 30-halving fallback).
                state%x_base(1:n) = state%x(1:n)
                state%iter = state%iter + 1
                ! Convergence tests at the accepted point.
                if (fmax <= state%ftol) then
                    state%phase = PHASE_DONE
                    action = MULTIROOT_DONE
                    return
                end if
                dxmax = state%lambda*inf_norm(state%dx, n)
                xmax = inf_norm_x(state%x_base, n)
                if (dxmax <= state%xtol*(xmax + state%xtol)) then
                    state%phase = PHASE_DONE
                    action = MULTIROOT_DONE
                    return
                end if
                if (state%iter >= state%max_iter) then
                    state%phase = PHASE_DONE
                    state%fail_code = MULTIROOT_RC_MAXITER
                    action = MULTIROOT_FAILED
                    return
                end if
                ! Next Newton step needs F,J at the accepted point, which the
                ! caller already supplied this call -> reuse it directly.
                call begin_newton(state, f, jac, amat, rhs, action)
                return
            else
                ! Reject: halve lambda and request the next trial point.
                state%lambda = 0.5_dp*state%lambda
                state%ls_count = state%ls_count + 1
                state%x(1:n) = state%x_base(1:n) + state%lambda*state%dx(1:n)
                action = MULTIROOT_NEED_FJ
                return
            end if

        case default
            ! PHASE_DONE or unknown: idempotent.
            if (state%fail_code /= 0) then
                action = MULTIROOT_FAILED
            else
                action = MULTIROOT_DONE
            end if
            return
        end select
    end subroutine multiroot_step

    ! Solve J dx = -F for the Newton direction, set lambda = 1, and request the
    ! first trial point x_base + dx.  On a singular Jacobian, fail.
    pure subroutine begin_newton(state, f, jac, amat, rhs, action)
        !$acc routine seq
        type(multiroot_rc_t), intent(inout) :: state
        real(dp), intent(in) :: f(state%n)
        real(dp), intent(in) :: jac(state%n, state%n)
        real(dp), intent(inout) :: amat(MULTIROOT_RC_MAX_N, MULTIROOT_RC_MAX_N)
        real(dp), intent(inout) :: rhs(MULTIROOT_RC_MAX_N)
        integer, intent(out) :: action
        integer :: n, info, i, j

        n = state%n
        do j = 1, n
            do i = 1, n
                amat(i, j) = jac(i, j)
            end do
        end do
        do i = 1, n
            rhs(i) = -f(i)
        end do
        call lu_solve(n, amat(1:n, 1:n), rhs(1:n), info)
        if (info /= LINALG_OK) then
            state%phase = PHASE_DONE
            state%fail_code = MULTIROOT_RC_SINGULAR
            action = MULTIROOT_FAILED
            return
        end if
        state%dx(1:n) = rhs
        state%g0 = 0.5_dp*dot_n(f, f, n)
        state%lambda = 1.0_dp
        state%ls_count = 0
        state%x_base(1:n) = state%x(1:n)
        state%x(1:n) = state%x_base(1:n) + state%lambda*state%dx(1:n)
        state%phase = PHASE_LINESEARCH
        action = MULTIROOT_NEED_FJ
    end subroutine begin_newton

    ! Infinity norm of a length-n vector (explicit-shape, no array temporaries).
    pure function inf_norm(v, n) result(m)
        !$acc routine seq
        real(dp), intent(in) :: v(n)
        integer, intent(in) :: n
        real(dp) :: m
        integer :: i
        m = 0.0_dp
        do i = 1, n
            if (abs(v(i)) > m) m = abs(v(i))
        end do
    end function inf_norm

    ! Infinity norm of the leading n entries of a fixed-size array.
    pure function inf_norm_x(v, n) result(m)
        !$acc routine seq
        real(dp), intent(in) :: v(MULTIROOT_RC_MAX_N)
        integer, intent(in) :: n
        real(dp) :: m
        integer :: i
        m = 0.0_dp
        do i = 1, n
            if (abs(v(i)) > m) m = abs(v(i))
        end do
    end function inf_norm_x

    ! Dot product of two length-n vectors.
    pure function dot_n(a, b, n) result(s)
        !$acc routine seq
        real(dp), intent(in) :: a(n), b(n)
        integer, intent(in) :: n
        real(dp) :: s
        integer :: i
        s = 0.0_dp
        do i = 1, n
            s = s + a(i)*b(i)
        end do
    end function dot_n

end module fortnum_multiroot_rc
