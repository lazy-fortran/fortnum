module bdf_fixtures
    use fortnum_kinds, only: dp
    use fortnum_ode_bdf, only: bdf_problem_t
    use fortnum_linalg, only: lu_factor, lu_solve_factored
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    private
    public :: problem_t, robertson_rhs, robertson_jacobian

    type, extends(bdf_problem_t) :: problem_t
        integer :: mode = 1, fail_setup = 0, fail_evaluate = 0
        logical :: fail_precondition = .false., zero_operator = .false.
        logical :: use_preconditioner = .true.
        real(dp) :: stiffness = 1000.0_dp, alpha = 0.0_dp
        real(dp), allocatable :: matrix(:, :), factors(:, :)
        integer, allocatable :: pivots(:)
    contains
        procedure :: evaluate
        procedure :: linearize
        procedure :: apply
        procedure :: precondition
    end type problem_t

contains

    subroutine evaluate(self, t, y, ydot, residual, status)
        class(problem_t), intent(inout) :: self
        real(dp), intent(in) :: t, y(:), ydot(:)
        real(dp), intent(out) :: residual(:)
        integer, intent(out) :: status
        real(dp) :: f(3), v(2), diffusion
        integer :: i, n

        status = 0
        if (self%fail_evaluate > 0) then
            self%fail_evaluate = self%fail_evaluate - 1
            residual = ieee_value(1.0_dp, ieee_quiet_nan)
            return
        end if
        select case (self%mode)
        case (1)
            v(1) = ydot(1) + y(1) - 10*y(2)
            v(2) = ydot(2) + self%stiffness*y(2)
            residual(1) = 2*v(1) + 0.3_dp*v(2)
            residual(2) = 0.1_dp*v(1) + v(2)
        case (2)
            residual = (1 + y**2)*(ydot + y) + &
                self%stiffness*(y - exp(-t)) + (y - exp(-t))**3
        case (3)
            call robertson_rhs(y, f)
            residual = ydot - f
        case (4)
            n = size(y)
            diffusion = 1.0e8_dp/(4*16.0_dp**2)*n**2
            do i = 1, n
                residual(i) = ydot(i) + y(i) + diffusion* &
                    (2*y(i) - y(modulo(i - 2, n) + 1) - y(modulo(i, n) + 1))
            end do
        case (5)
            residual = ydot - (1 + 2*t)
        end select
    end subroutine evaluate

    subroutine linearize(self, t, y, ydot, alpha, status)
        class(problem_t), intent(inout) :: self
        real(dp), intent(in) :: t, y(:), ydot(:), alpha
        integer, intent(out) :: status
        real(dp) :: diffusion
        integer :: i, n

        status = 11
        if (self%fail_setup > 0) then
            self%fail_setup = self%fail_setup - 1
            return
        end if
        n = size(y)
        if (allocated(self%matrix)) then
            if (size(self%matrix, 1) /= n) &
                deallocate (self%matrix, self%factors, self%pivots)
        end if
        if (.not. allocated(self%matrix)) &
            allocate (self%matrix(n, n), self%factors(n, n), self%pivots(n))
        self%matrix = 0.0_dp
        self%alpha = alpha
        select case (self%mode)
        case (1)
            self%matrix(1, 1) = 2*(alpha + 1)
            self%matrix(1, 2) = -20 + 0.3_dp*(alpha + self%stiffness)
            self%matrix(2, 1) = 0.1_dp*(alpha + 1)
            self%matrix(2, 2) = -1 + alpha + self%stiffness
        case (2)
            self%matrix(1, 1) = (1 + y(1)**2)*(alpha + 1) + &
                2*y(1)*(ydot(1) + y(1)) + &
                self%stiffness + 3*(y(1) - exp(-t))**2
        case (3)
            call robertson_jacobian(y, self%matrix)
            self%matrix = -self%matrix
            do i = 1, n
                self%matrix(i, i) = self%matrix(i, i) + alpha
            end do
        case (4)
            diffusion = 1.0e8_dp/(4*16.0_dp**2)*n**2
            do i = 1, n
                self%matrix(i, i) = alpha + 1 + 2*diffusion
                self%matrix(i, modulo(i - 2, n) + 1) = -diffusion
                self%matrix(i, modulo(i, n) + 1) = -diffusion
            end do
        case (5)
            self%matrix(1, 1) = alpha
        end select
        self%factors = self%matrix
        call lu_factor(n, self%factors, self%pivots, status)
    end subroutine linearize

    subroutine apply(self, input, output)
        class(problem_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = matmul(self%matrix, input)
        if (self%zero_operator) output = 0.0_dp
    end subroutine apply

    subroutine precondition(self, r, z, status)
        class(problem_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        real(dp) :: solved(size(r))

        z = r
        status = 0
        if (self%fail_precondition) then
            status = 12
        else if (self%use_preconditioner) then
            solved = r
            call lu_solve_factored(size(r), self%factors, self%pivots, solved, status)
            z = solved
        end if
    end subroutine precondition

    subroutine robertson_rhs(y, f)
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: f(:)

        f(1) = -0.04_dp*y(1) + 1.0e4_dp*y(2)*y(3)
        f(3) = 3.0e7_dp*y(2)**2
        f(2) = -f(1) - f(3)
    end subroutine robertson_rhs

    subroutine robertson_jacobian(y, j)
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: j(:, :)

        j(1, 1) = -0.04_dp
        j(1, 2) = 1.0e4_dp*y(3)
        j(1, 3) = 1.0e4_dp*y(2)
        j(3, 1) = 0.0_dp
        j(3, 2) = 6.0e7_dp*y(2)
        j(3, 3) = 0.0_dp
        j(2, :) = -j(1, :) - j(3, :)
    end subroutine robertson_jacobian
end module bdf_fixtures

program test_fortnum_ode_bdf
    use fortnum_kinds, only: dp
    use fortnum_ode_bdf, only: bdf_options_t, bdf_state_t, bdf_report_t, &
        bdf_init, bdf_integrate_to, BDF_LINEARIZE, &
        BDF_LINEAR_FAILURE, &
        BDF_NEWTON_LIMIT
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_linalg, only: dense_solve
    use bdf_fixtures, only: problem_t, robertson_rhs, robertson_jacobian
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none

    call check_order()
    call check_step_scaling()
    call check_estimator()
    call check_variable_steps()
    call check_error_control()
    call check_robertson_cost()
    call check_diffusion()
    call check_failures()
    print '(a)', 'fortnum_ode_bdf: all tests passed'

contains

    subroutine check_estimator()
        type(problem_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: status
        real(dp), parameter :: ratios(*) = [0.5_dp, 1.0_dp, 2.0_dp]
        real(dp), parameter :: t = 0.4_dp
        real(dp) :: h, k, true_error, estimated, weight, quotient, first
        integer :: i, refinement, direction, order

        problem%mode = 2
        problem%stiffness = 0.0_dp
        options%rtol = 0.01_dp
        options%atol = 1.0e-10_dp
        options%newton_tolerance = 1.0e-10_dp
        options%max_steps = 1
        do order = 1, 2
            do direction = -1, 1, 2
                do i = 1, size(ratios)
                    if (order == 1 .and. i > 1) cycle
                    do refinement = 1, 4
                        h = direction*0.02_dp/2**(refinement - 1)
                        k = h/ratios(i)
                        options%initial_step = abs(h)
                        call bdf_init(state, t, [exp(-t)], options, status)
                        if (order == 2) then
                            call bdf_init(state, t, [exp(-t)], options, status, &
                                derivative=[-exp(-t)])
                            state%previous_step = k
                            state%previous_y = [exp(-(t - k))]
                        end if
                        call bdf_integrate_to(problem, state, t + h, options, &
                            report, status)
                        call require(status%code == FORTNUM_OK, &
                            'exact-start nonlinear single step')
                        call require(report%accepted == 1, 'single accepted step')
                        call require(report%implicit_solves == 1, &
                            'error estimator requires no extra implicit solve')
                        if (order == 1) then
                            call require(report%mass_iterations > 0, &
                                'startup computes the consistent derivative once')
                        else
                            call require(report%mass_iterations == 0, &
                                'supplied derivative avoids the initial mass solve')
                        end if
                        weight = options%atol + options%rtol* &
                            max(exp(-t), abs(state%y(1)))
                        estimated = report%last_error*weight
                        true_error = abs(state%y(1) - exp(-(t + h)))
                        quotient = estimated/true_error
                        if (refinement == 1) first = abs(quotient - 1)
                        print '(a,i0,a,f4.1,a,es10.2,a,f10.6)', &
                            'LTE order=', order, ' omega=', ratios(i), &
                            ' h=', h, ' estimate/true=', quotient
                    end do
                    call require(abs(quotient - 1) < 0.1_dp, &
                        'local error estimate tends to the true local error')
                    call require(abs(quotient - 1) < first/2, &
                        'local error estimate improves with refinement')
                end do
            end do
        end do
    end subroutine check_estimator

    subroutine check_step_scaling()
        type(problem_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: status
        real(dp) :: exponent, tolerances(5), counts(5), x(5), z(5)
        integer :: i

        problem%mode = 2
        problem%stiffness = 0.0_dp
        options%initial_step = 0.01_dp
        do i = 1, 5
            tolerances(i) = 10.0_dp**(-i - 3)
            options%rtol = tolerances(i)
            options%atol = 0.01_dp*options%rtol
            call bdf_init(state, 0.0_dp, [1.0_dp], options, status)
            call bdf_integrate_to(problem, state, 2.0_dp, options, report, status)
            call require(status%code == FORTNUM_OK, 'smooth nonlinear scaling solve')
            counts(i) = real(report%accepted, dp)
            print '(a,es10.2,a,i0)', 'smooth tolerance=', tolerances(i), &
                ' steps=', report%accepted
        end do
        x = log(1/tolerances)
        z = log(counts)
        x = x - sum(x)/size(x)
        z = z - sum(z)/size(z)
        exponent = sum(x*z)/sum(x*x)
        print '(a,f8.4)', 'step-count tolerance exponent=', exponent
        call require(exponent > 0.25_dp .and. exponent < 0.45_dp, &
            'BDF2 step count scales with the cube root of inverse tolerance')
    end subroutine check_step_scaling

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (condition) return
        print '(a)', 'FAIL: '//message
        error stop 1
    end subroutine require

    subroutine check_order()
        type(problem_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: status
        real(dp) :: errors(3), initial(2), exact(2), h, ratio
        integer :: mode, refinement

        options%adaptive = .false.
        options%rtol = 1.0e-10_dp
        options%atol = 1.0e-12_dp
        do mode = 1, 2
            problem%mode = mode
            problem%use_preconditioner = mode == 2
            problem%stiffness = 1000.0_dp
            initial = 1.0_dp
            exact(1) = exp(-1.0_dp)*(1 + 10.0_dp/999)
            exact(2) = exp(-problem%stiffness)
            if (mode == 2) then
                problem%stiffness = 25.0_dp
                exact = exp(-1.0_dp)
            end if
            do refinement = 1, 3
                h = 1.0_dp/(32*2**(refinement - 1))
                options%initial_step = h
                options%max_step = h
                call bdf_init(state, 0.0_dp, initial(:3 - mode), options, status)
                call require(status%code == FORTNUM_OK, 'order initialization')
                call bdf_integrate_to(problem, state, 1.0_dp, options, report, status)
                call require(status%code == FORTNUM_OK, 'fixed-step stiff solve')
                call require(report%rejected_newton == 0, 'order mesh is not changed')
                call require(report%accepted == 32*2**(refinement - 1), &
                    'fixed mesh count')
                errors(refinement) = maxval(abs(state%y - exact(:3 - mode)))
                print '(a,i0,a,es12.4,a,es12.4)', 'order mode=', mode, &
                    ' h=', h, ' error=', errors(refinement)
            end do
            ratio = errors(2)/errors(3)
            call require(ratio > 3.5_dp .and. ratio < 4.5_dp, 'BDF2 order two')
        end do
    end subroutine check_order

    subroutine check_variable_steps()
        type(problem_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state, interleaved
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: status
        real(dp), parameter :: targets(*) = [0.14_dp, 0.21_dp, 0.23_dp, 0.35_dp, 0.5_dp]
        integer :: i

        problem%mode = 5
        options%adaptive = .false.
        options%initial_step = 0.1_dp
        call bdf_init(state, 0.1_dp, [1.11_dp], options, status)
        state%previous_y = [1.0_dp]
        state%previous_step = 0.1_dp
        call bdf_init(interleaved, 0.0_dp, [2.0_dp], options, status)
        do i = 1, size(targets)
            state%next_step = targets(i) - state%t
            call bdf_integrate_to(problem, state, targets(i), options, report, status)
            call require(status%code == FORTNUM_OK, 'nonuniform mesh solve')
            call require(abs(state%y(1) - (1 + state%t + state%t**2)) < 2.0e-12_dp, &
                'BDF2 differentiates a quadratic exactly on a nonuniform mesh')
            call bdf_integrate_to(problem, interleaved, targets(i), options, report, &
                status)
            call require(status%code == FORTNUM_OK, 'interleaved independent state')
        end do
        options%adaptive = .true.
        options%rtol = 1.0e-6_dp
        call bdf_integrate_to(problem, state, 0.0_dp, options, report, status)
        call require(status%code == FORTNUM_OK, 'backward integration restarts history')
        call require(abs(state%y(1) - 1.0_dp) < 1.0e-5_dp, 'backward exact solution')
    end subroutine check_variable_steps

    subroutine dense_reference(final_time, steps, result)
        real(dp), intent(in) :: final_time
        integer, intent(in) :: steps
        real(dp), intent(out) :: result(3)
        real(dp) :: y(3), f0(3), f(3), j(3, 3), r(3), d(3), h
        integer :: step, iteration, info, i

        result = [1.0_dp, 0.0_dp, 0.0_dp]
        h = final_time/steps
        do step = 1, steps
            call robertson_rhs(result, f0)
            y = result
            do iteration = 1, 12
                call robertson_rhs(y, f)
                r = y - result - h/2*(f0 + f)
                call robertson_jacobian(y, j)
                j = -h/2*j
                do i = 1, 3
                    j(i, i) = j(i, i) + 1
                end do
                r = -r
                call dense_solve(j, r, d, info)
                call require(info == 0, 'dense trapezoidal reference solve')
                y = y + d
                if (maxval(abs(d)) < 2.0e-16_dp) exit
            end do
            call require(iteration <= 12, 'dense reference Newton convergence')
            result = y
        end do
    end subroutine dense_reference

    subroutine check_error_control()
        type(problem_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: status
        real(dp) :: coarse(3), fine(3), finer(3), reference(3), error, previous_error
        integer :: decade

        call dense_reference(1.0_dp, 10000, coarse)
        call dense_reference(1.0_dp, 20000, fine)
        call dense_reference(1.0_dp, 40000, finer)
        reference = (4*finer - fine)/3
        call require(maxval(abs(reference - (4*fine - coarse)/3)) < 1.0e-11_dp, &
            'independent dense Robertson reference is converged')
        problem%mode = 3
        options%initial_step = 0.1_dp
        previous_error = huge(1.0_dp)
        do decade = 3, 7
            options%rtol = 10.0_dp**(-decade)
            options%atol = options%rtol*1.0e-4_dp
            call bdf_init(state, 0.0_dp, [1.0_dp, 0.0_dp, 0.0_dp], options, status)
            call bdf_integrate_to(problem, state, 1.0_dp, options, report, status)
            call require(status%code == FORTNUM_OK, 'adaptive Robertson solve')
            error = maxval(abs(state%y - reference))
            print '(a,es10.2,a,es12.4,a,i0,a,i0)', 'Robertson tolerance=', &
                options%rtol, ' error=', error, ' steps=', report%accepted, &
                ' rejected=', report%rejected_error + report%rejected_newton
            call require(error < options%rtol, 'Robertson error follows tolerance')
            call require(error < 0.4_dp*previous_error, &
                'each tolerance decade improves error')
            call require(abs(sum(state%y) - 1) < 1.0e-11_dp, 'Robertson mass balance')
            call require(report%rejected_error > 0, &
                'error test rejects the large initial step')
            previous_error = error
        end do
    end subroutine check_error_control

    subroutine check_robertson_cost()
        type(problem_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: status
        real(dp) :: coarse(3), fine(3), finer(3), reference(3), error, cost
        integer :: decade

        call dense_reference(4.0_dp, 20000, coarse)
        call dense_reference(4.0_dp, 40000, fine)
        call dense_reference(4.0_dp, 80000, finer)
        reference = (4*finer - fine)/3
        call require(maxval(abs(reference - (4*fine - coarse)/3)) < 1.0e-11_dp, &
            'independent long Robertson reference is converged')
        problem%mode = 3
        options%initial_step = 1.0e-6_dp
        do decade = 3, 7
            options%rtol = 10.0_dp**(-decade)
            options%atol = options%rtol*1.0e-4_dp
            call bdf_init(state, 0.0_dp, [1.0_dp, 0.0_dp, 0.0_dp], options, status)
            call bdf_integrate_to(problem, state, 4.0_dp, options, report, status)
            call require(status%code == FORTNUM_OK, 'Robertson cost solve')
            error = maxval(abs(state%y - reference))
            cost = real(report%implicit_solves, dp)/report%accepted
            print '(a,es10.2,a,es12.4,a,i0,a,f8.4)', 'Robertson cost tolerance=', &
                options%rtol, ' error=', error, ' steps=', report%accepted, &
                ' solves/accepted=', cost
            call require(error > 0.1_dp*options%rtol .and. &
                error < 10*options%rtol, 'Robertson achieved accuracy matches cost')
            call require(report%implicit_solves == report%accepted + &
                report%rejected_error + report%rejected_newton, &
                'one implicit solve per Robertson attempt')
            call require(cost < 1.2_dp, 'Robertson retries add little solve cost')
        end do
    end subroutine check_robertson_cost

    subroutine check_diffusion()
        type(problem_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: status
        real(dp), allocatable :: y(:)
        real(dp) :: pi, stiffness, started, elapsed
        integer :: refinement, n, i, counts(3)

        problem%mode = 4
        options%rtol = 1.0e-4_dp
        options%atol = 1.0e-7_dp
        options%initial_step = 0.01_dp
        pi = acos(-1.0_dp)
        do refinement = 1, 3
            n = 16*2**(refinement - 1)
            allocate (y(n))
            do i = 1, n
                y(i) = 1 + 0.01_dp*sin(2*pi*(i - 1)/n) + 0.01_dp*(-1)**i
            end do
            stiffness = 1 + 1.0e8_dp*(real(n, dp)/16)**2
            call cpu_time(started)
            call bdf_init(state, 0.0_dp, y, options, status)
            call bdf_integrate_to(problem, state, 1.0_dp, options, report, status)
            call cpu_time(elapsed)
            elapsed = elapsed - started
            call require(status%code == FORTNUM_OK, 'stiff periodic diffusion solve')
            call require(maxval(abs(state%y - exp(-1.0_dp))) < 2.0e-4_dp, &
                'diffusion Fourier solution after fast modes decay')
            counts(refinement) = report%accepted
            print '(a,i0,a,es10.2,a,i0,a,f8.3)', 'diffusion n=', n, &
                ' stiffness=', stiffness, ' steps=', report%accepted, ' cpu_s=', elapsed
            deallocate (y)
        end do
        call require(maxval(counts) < 300, 'diffusion step count stays bounded')
        call require(maxval(counts) < 2*minval(counts), &
            'refinement does not multiply steps')
    end subroutine check_diffusion

    subroutine check_failures()
        type(problem_t) :: problem
        type(bdf_options_t) :: options
        type(bdf_state_t) :: state
        type(bdf_report_t) :: report
        type(fortnum_status_t) :: status
        real(dp) :: invalid(1)
        integer :: failure

        call bdf_integrate_to(problem, state, 1.0_dp, options, report, status)
        call require(status%code == FORTNUM_DOMAIN_ERROR, &
            'uninitialized state is rejected')
        options%rtol = -1.0_dp
        call bdf_init(state, 0.0_dp, [1.0_dp], options, status)
        call require(status%code == FORTNUM_DOMAIN_ERROR, &
            'negative tolerance is rejected')
        options = bdf_options_t()
        invalid(1) = ieee_value(1.0_dp, ieee_quiet_nan)
        call bdf_init(state, 0.0_dp, invalid, options, status)
        call require(status%code == FORTNUM_DOMAIN_ERROR, &
            'nonfinite initial value is rejected')
        call bdf_init(state, 0.0_dp, [1.0_dp], options, status, derivative=invalid)
        call require(status%code == FORTNUM_DOMAIN_ERROR, &
            'nonfinite initial derivative is rejected')
        call bdf_init(state, 0.0_dp, [1.0_dp], options, status, &
            derivative=[-1.0_dp, -1.0_dp])
        call require(status%code == FORTNUM_DOMAIN_ERROR, &
            'wrong-sized initial derivative is rejected')
        problem%mode = 2
        problem%stiffness = 25.0_dp
        options%initial_step = 0.1_dp
        problem%fail_setup = 1
        problem%fail_evaluate = 1
        call bdf_init(state, 0.0_dp, [1.0_dp], options, status)
        call bdf_integrate_to(problem, state, 0.1_dp, options, report, status)
        call require(status%code == FORTNUM_OK, &
            'retry recovers from transient callback failures')
        call require(report%rejected_newton >= 2, &
            'failed nonlinear attempts are rejected')
        call require(abs(state%y(1) - exp(-0.1_dp)) < 1.0e-5_dp, &
            'retry solution oracle')

        options%max_steps = 1
        do failure = 1, 4
            problem%fail_setup = 0
            problem%fail_precondition = .false.
            problem%zero_operator = .false.
            options%max_newton = 12
            if (failure == 1) problem%fail_setup = 10
            if (failure == 2) problem%fail_precondition = .true.
            if (failure == 3) problem%zero_operator = .true.
            if (failure == 4) options%max_newton = 1
            call bdf_init(state, 0.0_dp, [1.0_dp], options, status)
            call bdf_integrate_to(problem, state, 0.1_dp, options, report, status)
            call require(status%code == FORTNUM_CONVERGENCE_ERROR, &
                'failure exhausts budget')
            call require(state%t == 0 .and. all(state%y == 1), &
                'rejection preserves endpoint')
            call require(state%previous_step == 0 .and. all(state%previous_y == 1), &
                'rejection preserves multistep history')
            call require(report%rejected_newton == 1, 'nonlinear rejection count')
            select case (failure)
            case (1)
                call require(report%last_failure == BDF_LINEARIZE, &
                    'linearize failure status')
                call require(report%callback_status == 11, 'callback failure preserved')
            case (2, 3)
                call require(report%last_failure == BDF_LINEAR_FAILURE, &
                    'linear failure status')
                if (failure == 2) &
                    call require(report%callback_status == 12, &
                    'preconditioner failure preserved')
            case (4)
                call require(report%last_failure == BDF_NEWTON_LIMIT, &
                    'Newton limit status')
            end select
        end do
        problem%zero_operator = .false.
        options = bdf_options_t()
        options%initial_step = 0.01_dp
        options%min_step = 0.1_dp
        call bdf_init(state, 0.0_dp, [1.0_dp], options, status)
        call bdf_integrate_to(problem, state, 0.1_dp, options, report, status)
        call require(status%code == FORTNUM_CONVERGENCE_ERROR, &
            'minimum step is enforced')
    end subroutine check_failures
end program test_fortnum_ode_bdf
