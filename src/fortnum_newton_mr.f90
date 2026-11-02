module fortnum_newton_mr
    !! Matrix-free inexact Newton-MR for a smooth objective with a symmetric,
    !! possibly indefinite or singular Hessian.
    !!
    !! Each iteration linearizes the problem at x and solves H d = -g by
    !! preconditioned MINRES (`real_minres`) to the forcing eta. Two merits:
    !!
    !! - NEWTON_MR_MERIT_OBJECTIVE (Liu and Roosta, "A Newton-MR algorithm
    !!   with complexity guarantees for nonconvex smooth unconstrained
    !!   optimization", 2022): with npc_exit, MINRES stops at the first
    !!   residual of nonpositive curvature and returns d = M r (a descent
    !!   direction of nonpositive curvature). A solution direction backtracks
    !!   from alpha = 1 under the Armijo rule on f; an NPC direction tracks
    !!   forward (doubling alpha while the rule holds and f decreases). Inside
    !!   the roundoff band of f the trapezoidal change (g(x) + g(x + s)).s/2
    !!   ranks the step instead.
    !! - NEWTON_MR_MERIT_GRADIENT (Roosta, Liu, Xu and Mahoney, "Newton-MR:
    !!   inexact Newton method with minimum residual sub-problem solver",
    !!   EURO J. Comput. Optim. 10, 2022), for stationary points that are not
    !!   minimizers (saddle points of a Lagrangian): the merit is
    !!   phi = g^T M g / 2 with the preconditioner M of the current
    !!   linearization, held fixed across the line search. Preconditioned
    !!   MINRES does not increase |H d + g|_M, so every iterate after the first
    !!   step has the slope (H d)^T M g = (|H d + g|_M^2 - |H d|_M^2
    !!   - |g|_M^2)/2 < 0 without a curvature condition; the NPC exit is not
    !!   used.
    !!
    !! Converged when the Euclidean |g| <= max(tolerance |g0|, absolute_
    !! tolerance). The objective object extends `real_linear_operator_t`
    !! (apply: the Hessian action at the last linearize point) and owns its
    !! state; the driver only reaches it through pointer adapters.
    use fortnum_kinds, only: dp
    use fortnum_krylov, only: real_linear_operator_t, real_minres, &
        KRYLOV_BREAKDOWN, KRYLOV_INVALID_ARGUMENT, KRYLOV_NONPOSITIVE_CURVATURE
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
        ieee_quiet_nan
    implicit none
    private

    integer, parameter, public :: NEWTON_MR_OK = 0
    integer, parameter, public :: NEWTON_MR_MAX_ITERATIONS = 1
    integer, parameter, public :: NEWTON_MR_BREAKDOWN = 2
    integer, parameter, public :: NEWTON_MR_LINE_SEARCH = 4
    integer, parameter, public :: NEWTON_MR_INVALID_ARGUMENT = 5
    integer, parameter, public :: NEWTON_MR_EVALUATION = 6
    !! Forcing eta = min(forcing_max, (|g|/|g0|)^(1/2)) or eta = forcing_max.
    integer, parameter, public :: NEWTON_MR_FORCING_SQRT = 0
    integer, parameter, public :: NEWTON_MR_FORCING_FIXED = 1
    integer, parameter, public :: NEWTON_MR_MERIT_OBJECTIVE = 0
    integer, parameter, public :: NEWTON_MR_MERIT_GRADIENT = 1

    type, abstract, extends(real_linear_operator_t), public :: newton_mr_problem_t
        !! Objective f and gradient g (evaluate), the Hessian point and
        !! preconditioner setup (linearize, default none), the Hessian action
        !! at that point (apply) and an SPD preconditioner z = M r ~ H^{-1} r
        !! (precondition, default identity). A nonzero status of
        !! precondition is a MINRES breakdown.
    contains
        procedure(newton_mr_evaluate_t), deferred :: evaluate
        procedure :: linearize => newton_mr_no_linearize
        procedure :: precondition => newton_mr_identity_precondition
    end type newton_mr_problem_t

    abstract interface
        subroutine newton_mr_evaluate_t(self, x, f, g, status)
            import :: newton_mr_problem_t, dp
            class(newton_mr_problem_t), intent(inout) :: self
            real(dp), intent(in) :: x(:)
            real(dp), intent(out) :: f, g(:)
            integer, intent(out) :: status
        end subroutine newton_mr_evaluate_t
    end interface

    type, public :: newton_mr_options_t
        real(dp) :: tolerance = 1.0e-10_dp
        real(dp) :: absolute_tolerance = 0.0_dp
        integer :: max_iterations = 50
        integer :: max_minres = 1000
        integer :: forcing = NEWTON_MR_FORCING_SQRT
        real(dp) :: forcing_max = 0.5_dp
        real(dp) :: armijo = 1.0e-4_dp
        integer :: max_line_search = 40
        logical :: npc_exit = .true.
        integer :: merit = NEWTON_MR_MERIT_OBJECTIVE
    end type newton_mr_options_t

    type, public :: newton_mr_report_t
        integer :: status = NEWTON_MR_INVALID_ARGUMENT
        integer :: iterations = 0, minres_iterations = 0, minres_max = 0
        integer :: npc_steps = 0, backtracks = 0, evaluations = 0
        real(dp) :: energy = 0.0_dp, initial_gradient = 0.0_dp
        real(dp) :: gradient = 0.0_dp
        !! Wall seconds: linearize, MINRES (operator and preconditioner
        !! applications) and evaluations.
        real(dp) :: linearize_s = 0.0_dp, minres_s = 0.0_dp
        real(dp) :: evaluate_s = 0.0_dp
    end type newton_mr_report_t

    type, extends(real_linear_operator_t) :: hessian_adapter_t
        class(newton_mr_problem_t), pointer :: problem => null()
    contains
        procedure :: apply => apply_hessian_adapter
    end type hessian_adapter_t

    type, extends(real_linear_operator_t) :: preconditioner_adapter_t
        class(newton_mr_problem_t), pointer :: problem => null()
    contains
        procedure :: apply => apply_preconditioner_adapter
    end type preconditioner_adapter_t

    public :: newton_mr_solve

contains

    subroutine newton_mr_no_linearize(self, x, status)
        class(newton_mr_problem_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: status

        status = 0
    end subroutine newton_mr_no_linearize

    subroutine newton_mr_identity_precondition(self, r, z, status)
        class(newton_mr_problem_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status

        z = r
        status = 0
    end subroutine newton_mr_identity_precondition

    subroutine apply_hessian_adapter(self, input, output)
        class(hessian_adapter_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        call self%problem%apply(input, output)
    end subroutine apply_hessian_adapter

    subroutine apply_preconditioner_adapter(self, input, output)
        class(preconditioner_adapter_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        integer :: status

        call self%problem%precondition(input, output, status)
        if (status /= 0) output = ieee_value(1.0_dp, ieee_quiet_nan)
    end subroutine apply_preconditioner_adapter

    subroutine newton_mr_solve(problem, x, options, report)
        !! Newton-MR from x; x is updated in place. report%status:
        !! NEWTON_MR_OK (gradient tolerance met), _MAX_ITERATIONS,
        !! _LINE_SEARCH (no acceptable step), _BREAKDOWN (linearize, MINRES or
        !! preconditioner failure, or no descent direction), _EVALUATION
        !! (failed evaluation at an accepted point or at x) or
        !! _INVALID_ARGUMENT. A MINRES iterate that stopped at its iteration
        !! cap or stagnated is used as an inexact step.
        class(newton_mr_problem_t), intent(inout), target :: problem
        real(dp), intent(inout) :: x(:)
        type(newton_mr_options_t), intent(in) :: options
        type(newton_mr_report_t), intent(out) :: report

        type(hessian_adapter_t) :: hessian
        type(preconditioner_adapter_t) :: preconditioner
        real(dp), allocatable :: g(:), gt(:), g_best(:), d(:), xt(:), mg(:)
        real(dp), allocatable :: work(:)
        real(dp) :: f, ft, gnorm, g0, eta, slope, alpha, gtol, residual, t0
        real(dp) :: f_best, alpha_best, merit, merit_t, merit_best
        integer :: iteration, info, minres_iterations, status, trial, forward, n
        logical :: npc, accepted, gradient_merit

        report%status = NEWTON_MR_INVALID_ARGUMENT
        n = size(x)
        if (n < 1 .or. options%max_iterations < 0 .or. &
            options%max_minres < 1 .or. options%max_line_search < 1) return
        if (.not. (options%forcing_max > 0.0_dp .and. options%forcing_max < 1.0_dp)) &
            return
        if (.not. (options%armijo > 0.0_dp .and. options%armijo < 1.0_dp)) return
        if (options%merit /= NEWTON_MR_MERIT_OBJECTIVE .and. &
            options%merit /= NEWTON_MR_MERIT_GRADIENT) return
        gradient_merit = options%merit == NEWTON_MR_MERIT_GRADIENT
        hessian%problem => problem
        preconditioner%problem => problem
        allocate (g(n), gt(n), g_best(n), d(n), xt(n), mg(n), work(n))

        report%status = NEWTON_MR_EVALUATION
        t0 = clock()
        call problem%evaluate(x, f, g, status)
        report%evaluate_s = report%evaluate_s + clock() - t0
        report%evaluations = 1
        if (status /= 0 .or. .not. ieee_is_finite(f) .or. &
            .not. all(ieee_is_finite(g))) return
        g0 = norm2(g)
        gnorm = g0
        report%initial_gradient = g0
        report%energy = f
        report%gradient = g0
        gtol = max(options%tolerance*g0, options%absolute_tolerance)
        report%status = NEWTON_MR_MAX_ITERATIONS
        do iteration = 1, options%max_iterations
            if (gnorm <= gtol) then
                report%status = NEWTON_MR_OK
                exit
            end if
            report%iterations = iteration
            t0 = clock()
            call problem%linearize(x, status)
            report%linearize_s = report%linearize_s + clock() - t0
            if (status /= 0) then
                report%status = NEWTON_MR_BREAKDOWN
                exit
            end if
            eta = options%forcing_max
            if (options%forcing == NEWTON_MR_FORCING_SQRT) &
                eta = min(eta, sqrt(gnorm/max(g0, tiny(1.0_dp))))
            ! A forcing far below what the stopping test needs only adds
            ! iterations at the roundoff floor.
            eta = max(eta, 0.5_dp*gtol/gnorm)
            d = 0.0_dp
            t0 = clock()
            if (options%npc_exit .and. .not. gradient_merit) then
                call real_minres(hessian, -g, d, eta, options%max_minres, info, &
                    minres_iterations, residual, preconditioner, work)
            else
                call real_minres(hessian, -g, d, eta, options%max_minres, info, &
                    minres_iterations, residual, preconditioner)
            end if
            report%minres_s = report%minres_s + clock() - t0
            report%minres_iterations = report%minres_iterations + minres_iterations
            report%minres_max = max(report%minres_max, minres_iterations)
            if (info == KRYLOV_BREAKDOWN .or. info == KRYLOV_INVALID_ARGUMENT) then
                report%status = NEWTON_MR_BREAKDOWN
                exit
            end if
            npc = info == KRYLOV_NONPOSITIVE_CURVATURE
            if (npc) then
                d = work
                report%npc_steps = report%npc_steps + 1
            end if
            if (gradient_merit) then
                ! phi = g^T M g / 2 at fixed M; slope (H d)^T M g.
                call problem%precondition(g, mg, status)
                if (status /= 0) then
                    report%status = NEWTON_MR_BREAKDOWN
                    exit
                end if
                merit = 0.5_dp*dot_product(g, mg)
                call problem%apply(d, work)
                slope = dot_product(work, mg)
                if (.not. slope < 0.0_dp) then
                    report%status = NEWTON_MR_BREAKDOWN
                    exit
                end if
            else
                merit = f
                slope = dot_product(g, d)
                if (.not. slope < 0.0_dp) then
                    ! Not a descent direction (roundoff, or an unconverged
                    ! first MINRES step): the preconditioned gradient -M g.
                    call problem%precondition(-g, mg, status)
                    slope = dot_product(g, mg)
                    if (status /= 0 .or. .not. slope < 0.0_dp) then
                        report%status = NEWTON_MR_BREAKDOWN
                        exit
                    end if
                    d = mg
                    npc = .true.
                end if
            end if

            accepted = .false.
            alpha = 1.0_dp
            alpha_best = 0.0_dp
            do trial = 1, options%max_line_search
                call trial_point(alpha, ft, merit_t, status)
                if (status == 0) then
                    if (armijo_ok(merit, merit_t, alpha*slope, alpha)) then
                        accepted = .true.
                        call keep(alpha, ft, merit_t)
                        ! NPC: forward tracking while the Armijo rule holds.
                        if (npc .and. trial == 1) then
                            do forward = 1, options%max_line_search
                                alpha = 2.0_dp*alpha
                                call trial_point(alpha, ft, merit_t, status)
                                if (status /= 0) exit
                                if (.not. armijo_ok(merit, merit_t, alpha*slope, &
                                    alpha)) exit
                                if (.not. merit_t < merit_best) exit
                                call keep(alpha, ft, merit_t)
                            end do
                        end if
                        exit
                    end if
                end if
                report%backtracks = report%backtracks + 1
                alpha = 0.5_dp*alpha
            end do
            if (.not. accepted) then
                report%status = NEWTON_MR_LINE_SEARCH
                exit
            end if
            x = x + alpha_best*d
            f = f_best
            g = g_best
            gnorm = norm2(g)
            report%energy = f
            report%gradient = gnorm
        end do
        if (gnorm <= gtol) report%status = NEWTON_MR_OK
        report%energy = f
        report%gradient = gnorm

    contains

        subroutine trial_point(step, value, trial_merit, trial_status)
            real(dp), intent(in) :: step
            real(dp), intent(out) :: value, trial_merit
            integer, intent(out) :: trial_status

            xt = x + step*d
            t0 = clock()
            call problem%evaluate(xt, value, gt, trial_status)
            report%evaluate_s = report%evaluate_s + clock() - t0
            report%evaluations = report%evaluations + 1
            if (trial_status == 0) then
                if (.not. ieee_is_finite(value) .or. &
                    .not. all(ieee_is_finite(gt))) trial_status = 1
            end if
            trial_merit = value
            if (trial_status /= 0 .or. .not. gradient_merit) return
            call problem%precondition(gt, work, trial_status)
            trial_merit = 0.5_dp*dot_product(gt, work)
        end subroutine trial_point

        subroutine keep(step, value, trial_merit)
            real(dp), intent(in) :: step, value, trial_merit

            alpha_best = step
            f_best = value
            merit_best = trial_merit
            g_best = gt
        end subroutine keep

        logical function armijo_ok(merit0, merit1, step_slope, step) result(ok)
            !! Armijo sufficient decrease; inside the roundoff band of the
            !! objective the trapezoidal change ranks the step instead.
            real(dp), intent(in) :: merit0, merit1, step_slope, step
            real(dp) :: band

            ok = .false.
            if (.not. ieee_is_finite(merit1)) return
            band = 64.0_dp*epsilon(1.0_dp)*max(abs(merit0), abs(merit1))
            if (gradient_merit .or. abs(merit1 - merit0) > band) then
                ok = merit1 <= merit0 + options%armijo*step_slope
            else
                ok = 0.5_dp*step*dot_product(g + gt, d) <= options%armijo*step_slope
            end if
        end function armijo_ok
    end subroutine newton_mr_solve

    real(dp) function clock()
        integer(8) :: count, rate

        call system_clock(count, rate)
        clock = real(count, dp)/real(rate, dp)
    end function clock
end module fortnum_newton_mr
