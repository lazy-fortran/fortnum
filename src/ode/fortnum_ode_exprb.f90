module fortnum_ode_exprb
    !! Autonomous exponential Rosenbrock methods; see docs/design/ode.md.
    use fortnum_kinds, only: dp
    use fortnum_krylov, only: real_linear_operator_t
    use fortnum_phi, only: phi_action, phi_options_t, phi_report_t, PHI_OK
    use fortnum_generated_exprb, only: exprb_controller
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    integer, parameter, public :: EXPRB_EULER = 2, EXPRB32 = 3
    integer, parameter, public :: EXPRB_SOLVED = 0, EXPRB_EVALUATION = 1
    integer, parameter, public :: EXPRB_LINEARIZE = 2, EXPRB_PHI_FAILURE = 3

    type, abstract, extends(real_linear_operator_t), public :: exprb_problem_t
    contains
        procedure(exprb_evaluate_t), deferred :: evaluate
        procedure(exprb_linearize_t), deferred :: linearize
    end type exprb_problem_t

    abstract interface
        subroutine exprb_evaluate_t(self, y, f, status)
            import :: exprb_problem_t, dp
            class(exprb_problem_t), intent(inout) :: self
            real(dp), intent(in) :: y(:)
            real(dp), intent(out) :: f(:)
            integer, intent(out) :: status
        end subroutine exprb_evaluate_t

        subroutine exprb_linearize_t(self, y, status)
            import :: exprb_problem_t, dp
            class(exprb_problem_t), intent(inout) :: self
            real(dp), intent(in) :: y(:)
            integer, intent(out) :: status
        end subroutine exprb_linearize_t
    end interface

    type, public :: exprb_options_t
        real(dp) :: rtol = 1.0e-5_dp, atol = 1.0e-8_dp
        real(dp) :: initial_step = 1.0e-3_dp
        real(dp) :: min_step = 0.0_dp, max_step = huge(1.0_dp)
        integer :: method = EXPRB32, max_steps = 100000
        logical :: adaptive = .true.
        type(phi_options_t) :: phi
    end type exprb_options_t

    type, public :: exprb_state_t
        real(dp) :: t = 0.0_dp, next_step = 0.0_dp
        real(dp), allocatable :: y(:)
    end type exprb_state_t

    type, public :: exprb_report_t
        integer :: accepted = 0, rejected_error = 0, rejected_failure = 0
        integer :: evaluations = 0, linearizations = 0, jacobian_actions = 0
        integer :: phi_calls = 0, phi_substeps = 0, phi_rejections = 0
        integer :: max_krylov_dimension = 0
        integer :: last_failure = EXPRB_SOLVED, callback_status = 0, phi_status = PHI_OK
        real(dp) :: last_error = 0.0_dp
    end type exprb_report_t

    public :: exprb_init, exprb_integrate_to

contains

    logical function valid_options(options) result(valid)
        type(exprb_options_t), intent(in) :: options

        valid = .false.
        if (.not. ieee_is_finite(options%rtol)) return
        if (.not. ieee_is_finite(options%atol)) return
        if (.not. ieee_is_finite(options%initial_step)) return
        if (.not. ieee_is_finite(options%min_step)) return
        if (.not. ieee_is_finite(options%max_step)) return
        if (options%rtol < 0.0_dp .or. options%atol <= 0.0_dp) return
        if (options%initial_step <= 0.0_dp .or. options%min_step < 0.0_dp) return
        if (options%max_step <= 0.0_dp) return
        if (options%max_step < options%min_step .or. options%max_steps < 1) return
        if (options%method /= EXPRB32 .and. options%method /= EXPRB_EULER) return
        if (.not. ieee_is_finite(options%phi%rtol)) return
        if (.not. ieee_is_finite(options%phi%atol)) return
        if (min(options%phi%rtol, options%phi%atol) < 0.0_dp) return
        if (max(options%phi%rtol, options%phi%atol) <= 0.0_dp) return
        if (options%phi%min_dimension < 1) return
        if (options%phi%max_dimension < options%phi%min_dimension) return
        if (options%phi%max_steps < 1 .or. options%phi%max_rejections < 0) return
        valid = .true.
    end function valid_options

    subroutine exprb_init(state, t, y, options, status)
        type(exprb_state_t), intent(out) :: state
        real(dp), intent(in) :: t, y(:)
        type(exprb_options_t), intent(in) :: options
        type(fortnum_status_t), intent(out) :: status

        call status_set(status, FORTNUM_DOMAIN_ERROR, 'Invalid exprb initial data')
        if (.not. valid_options(options)) return
        if (size(y) < 1 .or. .not. ieee_is_finite(t)) return
        if (.not. all(ieee_is_finite(y))) return
        state%t = t
        state%y = y
        state%next_step = min(options%initial_step, options%max_step)
        call status_set(status, FORTNUM_OK, '')
    end subroutine exprb_init

    subroutine exprb_integrate_to(problem, state, target, options, report, status)
        class(exprb_problem_t), intent(inout) :: problem
        type(exprb_state_t), intent(inout) :: state
        real(dp), intent(in) :: target
        type(exprb_options_t), intent(in) :: options
        type(exprb_report_t), intent(out) :: report
        type(fortnum_status_t), intent(out) :: status
        real(dp), allocatable :: f(:), stage(:), delta(:), jf(:), defect(:)
        real(dp), allocatable :: correction(:), high(:), vectors(:, :), weights(:)
        real(dp) :: direction, h, error, factor
        integer :: n, attempt, info

        report = exprb_report_t()
        call status_set(status, FORTNUM_DOMAIN_ERROR, 'Invalid exprb controls or state')
        if (.not. valid_options(options)) return
        if (.not. allocated(state%y)) return
        n = size(state%y)
        if (n < 1) return
        if (.not. all(ieee_is_finite(state%y))) return
        if (.not. ieee_is_finite(state%t)) return
        if (.not. ieee_is_finite(target)) return
        if (.not. ieee_is_finite(state%next_step)) return
        if (state%next_step <= 0.0_dp) return
        call status_set(status, FORTNUM_OK, '')
        if (target == state%t) return
        direction = sign(1.0_dp, target - state%t)
        allocate (f(n), stage(n), delta(n), jf(n), defect(n), correction(n))
        allocate (high(n), vectors(n, 0:3), weights(n))

        do attempt = 1, options%max_steps
            h = direction*min(state%next_step, options%max_step, abs(target - state%t))
            if (state%t + h == state%t .or. abs(h) < options%min_step) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    'exprb step underflow or minimum step reached')
                return
            end if
            call take_step()
            if (report%last_failure /= EXPRB_SOLVED) then
                report%rejected_failure = report%rejected_failure + 1
                state%next_step = 0.25_dp*abs(h)
                cycle
            end if
            error = 0.0_dp
            if (options%adaptive) then
                weights = options%atol + options%rtol*max(abs(state%y), abs(high))
                error = norm2(correction/weights)/sqrt(real(n, dp))
            end if
            if (.not. ieee_is_finite(error)) error = huge(1.0_dp)
            report%last_error = error
            factor = 1.0_dp
            if (options%adaptive) then
                factor = 2.0_dp
                if (error > 0.0_dp) call exprb_controller(error, factor)
                factor = max(0.2_dp, min(2.0_dp, factor))
            end if
            state%next_step = min(options%max_step, factor*abs(h))
            if (error > 1.0_dp) then
                report%rejected_error = report%rejected_error + 1
                cycle
            end if
            state%y = high
            state%t = state%t + h
            report%accepted = report%accepted + 1
            if (direction*(target - state%t) <= 0.0_dp) then
                state%t = target
                return
            end if
        end do
        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            'exprb attempt budget exhausted')

    contains

        subroutine take_step()
            report%last_failure = EXPRB_EVALUATION
            call problem%evaluate(state%y, f, info)
            report%evaluations = report%evaluations + 1
            if (info /= 0) report%callback_status = info
            if (info /= 0 .or. .not. all(ieee_is_finite(f))) return
            report%last_failure = EXPRB_LINEARIZE
            call problem%linearize(state%y, info)
            report%linearizations = report%linearizations + 1
            if (info /= 0) report%callback_status = info
            if (info /= 0) return
            vectors = 0.0_dp
            vectors(:, 1) = h*f
            call action(1, delta)
            if (report%last_failure /= EXPRB_SOLVED) return
            stage = state%y + delta
            report%last_failure = EXPRB_EVALUATION
            if (.not. all(ieee_is_finite(stage))) return
            high = stage
            correction = 0.0_dp
            if (options%adaptive .or. options%method == EXPRB32) then
                call problem%evaluate(stage, defect, info)
                report%evaluations = report%evaluations + 1
                if (info /= 0) report%callback_status = info
                if (info /= 0 .or. .not. all(ieee_is_finite(defect))) return
                ! Use the actual rounded stage displacement in the nonlinear remainder.
                delta = stage - state%y
                call problem%apply(delta, jf)
                report%jacobian_actions = report%jacobian_actions + 1
                report%last_failure = EXPRB_PHI_FAILURE
                if (.not. all(ieee_is_finite(jf))) return
                defect = defect - f - jf
                vectors = 0.0_dp
                vectors(:, 3) = 2*h*defect
                call action(3, correction)
                if (report%last_failure /= EXPRB_SOLVED) return
                if (options%method == EXPRB32) high = stage + correction
            end if
            report%last_failure = EXPRB_EVALUATION
            if (.not. all(ieee_is_finite(high))) return
            report%last_failure = EXPRB_SOLVED
        end subroutine take_step

        subroutine action(order, value)
            integer, intent(in) :: order
            real(dp), intent(out) :: value(:)
            type(phi_report_t) :: diagnostics

            call phi_action(problem, h, vectors(:, 0:order), value, info, &
                options%phi, diagnostics)
            report%phi_calls = report%phi_calls + 1
            report%jacobian_actions = report%jacobian_actions + diagnostics%matvecs
            report%phi_substeps = report%phi_substeps + diagnostics%substeps
            report%phi_rejections = report%phi_rejections + diagnostics%rejections
            report%max_krylov_dimension = &
                max(report%max_krylov_dimension, diagnostics%max_dimension)
            report%phi_status = info
            report%last_failure = EXPRB_SOLVED
            if (info /= PHI_OK) report%last_failure = EXPRB_PHI_FAILURE
        end subroutine action
    end subroutine exprb_integrate_to
end module fortnum_ode_exprb
