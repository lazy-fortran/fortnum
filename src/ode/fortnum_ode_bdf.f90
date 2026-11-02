module fortnum_ode_bdf
    !! Matrix-free variable-step BDF2 for F(t,y,ydot) = M(t,y) ydot + R(t,y).
    !! See docs/design/ode.md for the implicit callback and error contracts.
    use fortnum_kinds, only: dp
    use fortnum_krylov, only: real_linear_operator_t, real_fgmres, &
        real_fgmres_operator, KRYLOV_OK
    use fortnum_generated_bdf2, only: bdf2_weights, bdf_predict, bdf_error
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    private

    integer, parameter, public :: BDF_SOLVED = 0, BDF_NEWTON_LIMIT = 1
    integer, parameter, public :: BDF_LINEAR_FAILURE = 2, BDF_EVALUATION = 3
    integer, parameter, public :: BDF_LINE_SEARCH = 4, BDF_LINEARIZE = 5

    type, abstract, extends(real_linear_operator_t), public :: bdf_problem_t
    contains
        procedure(bdf_evaluate_t), deferred :: evaluate
        procedure(bdf_linearize_t), deferred :: linearize
        procedure :: precondition => identity_precondition
    end type bdf_problem_t

    abstract interface
        subroutine bdf_evaluate_t(self, t, y, ydot, residual, status)
            import :: bdf_problem_t, dp
            class(bdf_problem_t), intent(inout) :: self
            real(dp), intent(in) :: t, y(:), ydot(:)
            real(dp), intent(out) :: residual(:)
            integer, intent(out) :: status
        end subroutine bdf_evaluate_t

        subroutine bdf_linearize_t(self, t, y, ydot, alpha, status)
            import :: bdf_problem_t, dp
            class(bdf_problem_t), intent(inout) :: self
            real(dp), intent(in) :: t, y(:), ydot(:), alpha
            integer, intent(out) :: status
        end subroutine bdf_linearize_t
    end interface

    type, public :: bdf_options_t
        real(dp) :: rtol = 1.0e-5_dp, atol = 1.0e-8_dp
        real(dp) :: initial_step = 1.0e-3_dp
        real(dp) :: min_step = 0.0_dp, max_step = huge(1.0_dp)
        real(dp) :: newton_tolerance = 0.03_dp
        real(dp) :: forcing_min = 1.0e-4_dp
        integer :: max_steps = 100000, max_newton = 12
        integer :: max_krylov = 200, restart = 30, max_backtracks = 10
        logical :: adaptive = .true.
    end type bdf_options_t

    type, public :: bdf_state_t
        real(dp) :: t = 0.0_dp, next_step = 0.0_dp, previous_step = 0.0_dp
        real(dp), allocatable :: y(:), previous_y(:), derivative(:)
        logical :: derivative_valid = .false.
    end type bdf_state_t

    type, public :: bdf_report_t
        integer :: accepted = 0, rejected_error = 0, rejected_newton = 0
        integer :: newton_iterations = 0, krylov_iterations = 0, evaluations = 0
        integer :: implicit_solves = 0, mass_iterations = 0
        integer :: last_failure = BDF_SOLVED, callback_status = 0
        integer :: linear_status = KRYLOV_OK
        real(dp) :: last_error = 0.0_dp
    end type bdf_report_t

    type, extends(real_linear_operator_t) :: preconditioner_t
        class(bdf_problem_t), pointer :: problem => null()
        integer :: status = 0
    contains
        procedure :: apply => apply_preconditioner
    end type preconditioner_t

    type :: newton_workspace_t
        real(dp), allocatable :: derivative(:), r(:), trial(:), trial_r(:), d(:)
        real(dp), allocatable :: rhs(:)
    end type newton_workspace_t

    public :: bdf_init, bdf_integrate_to

contains

    subroutine identity_precondition(self, r, z, status)
        class(bdf_problem_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status

        z = r
        status = 0
    end subroutine identity_precondition

    subroutine apply_preconditioner(self, input, output)
        class(preconditioner_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        integer :: status

        call self%problem%precondition(input, output, status)
        if (status /= 0) then
            self%status = status
            output = ieee_value(1.0_dp, ieee_quiet_nan)
        end if
    end subroutine apply_preconditioner

    logical function valid_options(options) result(valid)
        type(bdf_options_t), intent(in) :: options

        valid = .false.
        if (.not. ieee_is_finite(options%rtol)) return
        if (.not. ieee_is_finite(options%atol)) return
        if (.not. ieee_is_finite(options%initial_step)) return
        if (.not. ieee_is_finite(options%min_step)) return
        if (.not. ieee_is_finite(options%max_step)) return
        if (.not. ieee_is_finite(options%newton_tolerance)) return
        if (.not. ieee_is_finite(options%forcing_min)) return
        if (options%rtol < 0.0_dp .or. options%atol <= 0.0_dp) return
        if (options%initial_step <= 0.0_dp .or. options%min_step < 0.0_dp) return
        if (options%max_step <= 0.0_dp) return
        if (options%max_step < options%min_step) return
        if (options%newton_tolerance <= 0.0_dp .or. &
            options%newton_tolerance >= 1.0_dp) return
        if (options%forcing_min <= 0.0_dp .or. options%forcing_min > 0.1_dp) return
        if (options%max_steps < 1 .or. options%max_newton < 1) return
        if (options%restart < 1 .or. options%max_krylov < options%restart) return
        if (options%max_backtracks < 0) return
        valid = .true.
    end function valid_options

    subroutine bdf_init(state, t, y, options, status, derivative)
        type(bdf_state_t), intent(out) :: state
        real(dp), intent(in) :: t, y(:)
        type(bdf_options_t), intent(in) :: options
        type(fortnum_status_t), intent(out) :: status
        real(dp), intent(in), optional :: derivative(:)

        call status_set(status, FORTNUM_DOMAIN_ERROR, 'Invalid BDF initial data')
        if (.not. valid_options(options)) return
        if (size(y) < 1 .or. .not. ieee_is_finite(t)) return
        if (.not. all(ieee_is_finite(y))) return
        if (present(derivative)) then
            if (size(derivative) /= size(y)) return
            if (.not. all(ieee_is_finite(derivative))) return
        end if
        state%t = t
        state%y = y
        state%previous_y = y
        allocate (state%derivative(size(y)), source=0.0_dp)
        if (present(derivative)) then
            state%derivative = derivative
            state%derivative_valid = .true.
        end if
        state%next_step = min(options%initial_step, options%max_step)
        call status_set(status, FORTNUM_OK, '')
    end subroutine bdf_init

    subroutine bdf_integrate_to(problem, state, target, options, report, status)
        class(bdf_problem_t), intent(inout), target :: problem
        type(bdf_state_t), intent(inout) :: state
        real(dp), intent(in) :: target
        type(bdf_options_t), intent(in) :: options
        type(bdf_report_t), intent(out) :: report
        type(fortnum_status_t), intent(out) :: status
        type(newton_workspace_t) :: work
        real(dp), allocatable :: high(:), prediction(:), estimate(:)
        real(dp), allocatable :: offset(:), weights(:)
        real(dp) :: h, ratio, alpha, beta, milne, direction, factor, error, power
        integer :: n, attempt, failure

        report = bdf_report_t()
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'Invalid BDF controls or state')
        if (.not. valid_options(options)) return
        if (.not. allocated(state%y)) return
        if (.not. allocated(state%previous_y)) return
        if (.not. allocated(state%derivative)) return
        n = size(state%y)
        if (n < 1 .or. size(state%previous_y) /= n) return
        if (size(state%derivative) /= n) return
        if (.not. all(ieee_is_finite(state%y))) return
        if (.not. all(ieee_is_finite(state%previous_y))) return
        if (.not. all(ieee_is_finite(state%derivative))) return
        if (.not. ieee_is_finite(state%t)) return
        if (.not. ieee_is_finite(target)) return
        if (.not. ieee_is_finite(state%next_step)) return
        if (.not. ieee_is_finite(state%previous_step)) return
        if (state%next_step <= 0.0_dp) return
        call status_set(status, FORTNUM_OK, '')
        if (target == state%t) return
        direction = sign(1.0_dp, target - state%t)
        if (direction*state%previous_step < 0.0_dp) state%previous_step = 0.0_dp
        allocate (high(n), prediction(n), estimate(n), offset(n), weights(n))
        allocate (work%derivative(n), work%r(n), work%trial(n), &
            work%trial_r(n), work%d(n), work%rhs(n))

        do attempt = 1, options%max_steps
            h = direction*min(state%next_step, options%max_step, abs(target - state%t))
            if (state%previous_step /= 0.0_dp) &
                h = direction*min(abs(h), 2.0_dp*abs(state%previous_step))
            if (state%t + h == state%t .or. abs(h) < options%min_step) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    'BDF step underflow or minimum step reached')
                return
            end if
            if (.not. state%derivative_valid) then
                call initial_derivative(problem, state, options, report, failure)
                if (failure /= BDF_SOLVED) then
                    report%last_failure = failure
                    report%rejected_newton = report%rejected_newton + 1
                    state%next_step = 0.25_dp*abs(h)
                    cycle
                end if
            end if
            offset = 0.0_dp
            alpha = 1.0_dp/h
            ratio = 0.0_dp
            milne = 0.5_dp
            power = 0.5_dp
            if (state%previous_step /= 0.0_dp) then
                ratio = h/state%previous_step
                call bdf2_weights(h, ratio, alpha, beta, milne)
                offset = beta*(state%y - state%previous_y)
                power = 1.0_dp/3.0_dp
            end if
            call bdf_predict(h, ratio, state%y, state%previous_y, &
                state%derivative, prediction)
            high = prediction
            call solve_step(state%t + h, state%y, alpha, offset, high, failure)
            report%last_failure = failure
            if (failure /= BDF_SOLVED) then
                report%rejected_newton = report%rejected_newton + 1
                state%next_step = 0.25_dp*abs(h)
                cycle
            end if
            error = 0.0_dp
            if (options%adaptive) then
                weights = options%atol + options%rtol*max(abs(state%y), abs(high))
                call bdf_error(high, prediction, milne, estimate)
                error = norm2(estimate/weights)/sqrt(real(n, dp))
            end if
            report%last_error = error
            if (.not. ieee_is_finite(error)) error = huge(1.0_dp)
            factor = 2.0_dp
            if (.not. options%adaptive) factor = 1.0_dp
            if (error > 0.0_dp) &
                factor = max(0.2_dp, min(2.0_dp, 0.6_dp*error**(-power)))
            state%next_step = min(options%max_step, factor*abs(h))
            if (error > 1.0_dp) then
                report%rejected_error = report%rejected_error + 1
                cycle
            end if
            state%previous_y = state%y
            state%y = high
            state%derivative = work%derivative
            state%previous_step = h
            state%t = state%t + h
            report%accepted = report%accepted + 1
            if (direction*(target - state%t) <= 0.0_dp) then
                state%t = target
                return
            end if
        end do
        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            'BDF attempt budget exhausted')

    contains

        subroutine solve_step(t, base, alpha, offset, x, failure)
            real(dp), intent(in) :: t, base(:), alpha, offset(:)
            real(dp), intent(inout) :: x(:)
            integer, intent(out) :: failure

            report%implicit_solves = report%implicit_solves + 1
            call newton_step(problem, t, base, alpha, offset, x, options, &
                work, report, failure)
        end subroutine solve_step
    end subroutine bdf_integrate_to

    subroutine initial_derivative(problem, state, options, report, failure)
        class(bdf_problem_t), intent(inout) :: problem
        type(bdf_state_t), intent(inout) :: state
        type(bdf_options_t), intent(in) :: options
        type(bdf_report_t), intent(inout) :: report
        integer, intent(out) :: failure
        real(dp) :: zero(size(state%y)), residual(size(state%y))
        real(dp) :: rhs(size(state%y)), slope(size(state%y)), scale, linear_residual
        integer :: iterations, callback_failure

        failure = BDF_EVALUATION
        zero = 0.0_dp
        call problem%evaluate(state%t, state%y, zero, residual, &
            report%callback_status)
        report%evaluations = report%evaluations + 1
        if (report%callback_status /= 0) return
        if (.not. all(ieee_is_finite(residual))) return
        scale = norm2(residual)
        if (.not. ieee_is_finite(scale)) return
        slope = 0.0_dp
        if (scale > 0.0_dp) then
            rhs = -residual/scale
            callback_failure = 0
            call real_fgmres_operator(mass_action, identity, rhs, slope, &
                1.0e-8_dp, options%max_krylov, options%restart, &
                report%linear_status, iterations, linear_residual)
            report%mass_iterations = report%mass_iterations + iterations
            report%callback_status = callback_failure
            failure = BDF_LINEAR_FAILURE
            if (report%linear_status /= KRYLOV_OK .or. callback_failure /= 0) return
            slope = scale*slope
            if (.not. all(ieee_is_finite(slope))) return
        end if
        state%derivative = slope
        state%derivative_valid = .true.
        failure = BDF_SOLVED

    contains

        subroutine mass_action(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)
            integer :: info

            ! F is affine in ydot, so this is M input without a Jacobian setup.
            call problem%evaluate(state%t, state%y, input, output, info)
            report%evaluations = report%evaluations + 1
            if (info /= 0) callback_failure = info
            output = output - residual
            if (callback_failure /= 0) &
                output = ieee_value(1.0_dp, ieee_quiet_nan)
        end subroutine mass_action

        subroutine identity(input, output)
            real(dp), intent(in) :: input(:)
            real(dp), intent(out) :: output(:)

            output = input
        end subroutine identity
    end subroutine initial_derivative

    subroutine newton_step(problem, t, base, alpha, offset, x, options, &
            work, report, failure)
        class(bdf_problem_t), intent(inout), target :: problem
        real(dp), intent(in) :: t, base(:), alpha, offset(:)
        real(dp), intent(inout) :: x(:)
        type(bdf_options_t), intent(in) :: options
        type(newton_workspace_t), intent(inout) :: work
        type(bdf_report_t), intent(inout) :: report
        integer, intent(out) :: failure
        type(preconditioner_t) :: preconditioner
        real(dp) :: rnorm, next_norm, eta, safeguard, linear_residual
        real(dp) :: damping, correction_norm
        integer :: iteration, backtrack, iterations, info
        logical :: small, accepted

        preconditioner%problem => problem
        failure = BDF_EVALUATION
        work%derivative = alpha*(x - base) + offset
        call evaluate(x, work%r)
        if (report%callback_status /= 0 .or. .not. all(ieee_is_finite(work%r))) return
        rnorm = norm2(work%r)
        eta = 0.1_dp
        do iteration = 1, options%max_newton
            if (rnorm == 0.0_dp) then
                failure = BDF_SOLVED
                return
            end if
            failure = BDF_EVALUATION
            if (.not. ieee_is_finite(rnorm)) return
            failure = BDF_LINEARIZE
            call problem%linearize(t, x, work%derivative, alpha, report%callback_status)
            if (report%callback_status /= 0) return
            work%d = 0.0_dp
            ! Normalize so FGMRES's absolute floor still enforces relative forcing.
            work%rhs = -work%r/rnorm
            call real_fgmres(problem, preconditioner, work%rhs, work%d, eta, &
                options%max_krylov, options%restart, info, iterations, &
                linear_residual)
            report%linear_status = info
            report%krylov_iterations = report%krylov_iterations + iterations
            report%callback_status = preconditioner%status
            failure = BDF_LINEAR_FAILURE
            if (info /= KRYLOV_OK .or. preconditioner%status /= 0) return
            work%d = rnorm*work%d
            if (.not. all(ieee_is_finite(work%d))) return
            work%rhs = options%atol + options%rtol*max(abs(base), abs(x))
            correction_norm = norm2(work%d/work%rhs)/sqrt(real(size(x), dp))
            small = correction_norm <= options%newton_tolerance
            damping = 1.0_dp
            accepted = .false.
            do backtrack = 0, options%max_backtracks
                work%trial = x + damping*work%d
                work%derivative = alpha*(work%trial - base) + offset
                call evaluate(work%trial, work%trial_r)
                if (report%callback_status == 0 .and. &
                    all(ieee_is_finite(work%trial_r))) then
                    next_norm = norm2(work%trial_r)
                    if (ieee_is_finite(next_norm)) then
                        accepted = next_norm <= (1.0_dp - 1.0e-4_dp*damping)*rnorm
                        ! A converged correction need not reduce a residual at roundoff.
                        if (small) accepted = .true.
                    end if
                end if
                if (accepted) exit
                damping = damping/2
            end do
            failure = BDF_LINE_SEARCH
            if (.not. accepted) return
            x = work%trial
            work%r = work%trial_r
            report%newton_iterations = report%newton_iterations + 1
            if ((small .and. damping == 1.0_dp) .or. next_norm == 0.0_dp) then
                failure = BDF_SOLVED
                return
            end if
            ! Eisenstat--Walker choice 2, gamma=0.9, exponent=1.5, with safeguard.
            safeguard = 0.9_dp*eta**1.5_dp
            eta = 0.9_dp*(next_norm/rnorm)**1.5_dp
            if (safeguard > 0.1_dp) eta = max(eta, safeguard)
            eta = min(0.9_dp, max(options%forcing_min, eta))
            rnorm = next_norm
        end do
        failure = BDF_NEWTON_LIMIT

    contains

        subroutine evaluate(y, residual)
            real(dp), intent(in) :: y(:)
            real(dp), intent(out) :: residual(:)

            call problem%evaluate(t, y, work%derivative, residual, &
                report%callback_status)
            report%evaluations = report%evaluations + 1
        end subroutine evaluate
    end subroutine newton_step
end module fortnum_ode_bdf
