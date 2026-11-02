module exprb_fixtures
    use fortnum_kinds, only: dp
    use fortnum_ode_exprb, only: exprb_problem_t
    use fortnum_ode_bdf, only: bdf_problem_t
    use fortnum_linalg, only: lu_factor, lu_solve_factored
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    private
    public :: problem_t, implicit_problem_t, rhs, jacobian, initial, exact_rd

    type, extends(exprb_problem_t) :: problem_t
        integer :: mode = 1, fail_rhs = 0, fail_setup = 0, fail_action = 0
        integer :: actions = 0, evaluations = 0, linearizations = 0, fail_at = 0
        real(dp) :: diffusion = 1.0_dp, mu = 100.0_dp
        real(dp), allocatable :: point(:)
    contains
        procedure :: evaluate
        procedure :: linearize
        procedure :: apply
    end type problem_t

    type, extends(bdf_problem_t) :: implicit_problem_t
        integer :: mode = 1, actions = 0, evaluations = 0
        real(dp) :: diffusion = 1.0_dp, mu = 100.0_dp
        real(dp), allocatable :: matrix(:, :), factors(:, :)
        integer, allocatable :: pivots(:)
    contains
        procedure :: evaluate => implicit_evaluate
        procedure :: linearize => implicit_linearize
        procedure :: apply => implicit_apply
        procedure :: precondition
    end type implicit_problem_t

contains

    subroutine rhs(mode, diffusion, mu, y, f)
        integer, intent(in) :: mode
        real(dp), intent(in) :: diffusion, mu, y(:)
        real(dp), intent(out) :: f(:)
        real(dp) :: a, s, lambda, d, pi
        integer :: i, n

        pi = acos(-1.0_dp)
        select case (mode)
        case (1)
            f = -y**2
        case (2)
            n = size(y) - 1
            d = diffusion*(n + 1)**2
            lambda = 4*d*sin(pi/(2*(n + 1)))**2
            a = 1/(1 + y(n + 1))
            f(n + 1) = 1.0_dp
            do i = 1, n
                s = sin(pi*i/(n + 1))
                f(i) = -2*d*y(i) + y(i)**2 + &
                    (lambda*a - a**2)*s - a**2*s**2
                if (i > 1) f(i) = f(i) + d*y(i - 1)
                if (i < n) f(i) = f(i) + d*y(i + 1)
            end do
        case (3)
            f(1) = y(2)
            f(2) = mu*(1 - y(1)**2)*y(2) - y(1)
        case (4)
            ! Periodic diffusion retains a slow mode under mesh refinement.
            n = size(y)
            d = diffusion*n**2
            do i = 1, n
                f(i) = -y(i) + d*(y(modulo(i - 2, n) + 1) - &
                    2*y(i) + y(modulo(i, n) + 1))
            end do
        end select
    end subroutine rhs

    subroutine jacobian(mode, diffusion, mu, y, v, jv)
        integer, intent(in) :: mode
        real(dp), intent(in) :: diffusion, mu, y(:), v(:)
        real(dp), intent(out) :: jv(:)
        real(dp) :: a, s, lambda, d, pi
        integer :: i, n

        pi = acos(-1.0_dp)
        select case (mode)
        case (1)
            jv = -2*y*v
        case (2)
            n = size(y) - 1
            d = diffusion*(n + 1)**2
            lambda = 4*d*sin(pi/(2*(n + 1)))**2
            a = 1/(1 + y(n + 1))
            jv(n + 1) = 0.0_dp
            do i = 1, n
                s = sin(pi*i/(n + 1))
                jv(i) = (-2*d + 2*y(i))*v(i) + &
                    ((-lambda*a**2 + 2*a**3)*s + 2*a**3*s**2)*v(n + 1)
                if (i > 1) jv(i) = jv(i) + d*v(i - 1)
                if (i < n) jv(i) = jv(i) + d*v(i + 1)
            end do
        case (3)
            jv(1) = v(2)
            jv(2) = (-2*mu*y(1)*y(2) - 1)*v(1) + &
                mu*(1 - y(1)**2)*v(2)
        case (4)
            n = size(y)
            d = diffusion*n**2
            do i = 1, n
                jv(i) = -v(i) + d*(v(modulo(i - 2, n) + 1) - &
                    2*v(i) + v(modulo(i, n) + 1))
            end do
        end select
    end subroutine jacobian

    function initial(mode, n) result(y)
        integer, intent(in) :: mode, n
        real(dp) :: y(n), pi
        integer :: i

        pi = acos(-1.0_dp)
        y = 1.0_dp
        select case (mode)
        case (2)
            y = exact_rd(n, 0.0_dp)
        case (3)
            y = [2.0_dp, 0.0_dp]
        case (4)
            do i = 1, n
                y(i) = 1 + 0.1_dp*cos(2*pi*(i - 1)/n) + 0.05_dp*(-1)**i
            end do
        end select
    end function initial

    function exact_rd(n, t) result(y)
        integer, intent(in) :: n
        real(dp), intent(in) :: t
        real(dp) :: y(n)
        integer :: i

        do i = 1, n - 1
            y(i) = sin(acos(-1.0_dp)*i/n)/(1 + t)
        end do
        y(n) = t
    end function exact_rd

    subroutine evaluate(self, y, f, status)
        class(problem_t), intent(inout) :: self
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: f(:)
        integer, intent(out) :: status

        self%evaluations = self%evaluations + 1
        status = 0
        if (self%evaluations == self%fail_at) then
            status = 19
            f = 0.0_dp
            return
        end if
        call rhs(self%mode, self%diffusion, self%mu, y, f)
        if (self%fail_rhs > 0) then
            self%fail_rhs = self%fail_rhs - 1
            f = ieee_value(1.0_dp, ieee_quiet_nan)
        end if
    end subroutine evaluate

    subroutine linearize(self, y, status)
        class(problem_t), intent(inout) :: self
        real(dp), intent(in) :: y(:)
        integer, intent(out) :: status

        self%linearizations = self%linearizations + 1
        self%point = y
        status = 0
        if (self%fail_setup > 0) then
            self%fail_setup = self%fail_setup - 1
            status = 17
        end if
    end subroutine linearize

    subroutine apply(self, input, output)
        class(problem_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        self%actions = self%actions + 1
        call jacobian(self%mode, self%diffusion, self%mu, self%point, input, output)
        if (self%fail_action > 0) then
            self%fail_action = self%fail_action - 1
            output = ieee_value(1.0_dp, ieee_quiet_nan)
        end if
    end subroutine apply

    subroutine implicit_evaluate(self, t, y, ydot, residual, status)
        class(implicit_problem_t), intent(inout) :: self
        real(dp), intent(in) :: t, y(:), ydot(:)
        real(dp), intent(out) :: residual(:)
        integer, intent(out) :: status

        self%evaluations = self%evaluations + 1
        call rhs(self%mode, self%diffusion, self%mu, y, residual)
        residual = ydot - residual
        status = 0
    end subroutine implicit_evaluate

    subroutine implicit_linearize(self, t, y, ydot, alpha, status)
        class(implicit_problem_t), intent(inout) :: self
        real(dp), intent(in) :: t, y(:), ydot(:), alpha
        integer, intent(out) :: status
        real(dp) :: v(size(y))
        integer :: i, n

        n = size(y)
        if (.not. allocated(self%matrix)) &
            allocate (self%matrix(n, n), self%factors(n, n), self%pivots(n))
        do i = 1, n
            v = 0.0_dp
            v(i) = 1.0_dp
            call jacobian(self%mode, self%diffusion, self%mu, y, v, self%matrix(:, i))
        end do
        self%actions = self%actions + n
        self%matrix = -self%matrix
        do i = 1, n
            self%matrix(i, i) = self%matrix(i, i) + alpha
        end do
        self%factors = self%matrix
        call lu_factor(n, self%factors, self%pivots, status)
    end subroutine implicit_linearize

    subroutine implicit_apply(self, input, output)
        class(implicit_problem_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        self%actions = self%actions + 1
        output = matmul(self%matrix, input)
    end subroutine implicit_apply

    subroutine precondition(self, r, z, status)
        class(implicit_problem_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status

        z = r
        call lu_solve_factored(size(r), self%factors, self%pivots, z, status)
    end subroutine precondition
end module exprb_fixtures

program test_fortnum_ode_exprb
    use, intrinsic :: iso_fortran_env, only: int64
    use fortnum_kinds, only: dp
    use fortnum_ode_exprb, only: exprb_options_t, exprb_state_t, exprb_report_t, &
        exprb_init, exprb_integrate_to, EXPRB_EULER, EXPRB32
    use fortnum_ode_bdf, only: bdf_options_t, bdf_state_t, bdf_report_t, &
        bdf_init, bdf_integrate_to
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use exprb_fixtures, only: problem_t, implicit_problem_t, rhs, initial, exact_rd
    implicit none

    call orders()
    call tolerance_control()
    call diffusion_refinement()
    call rejection_reuse()
    call failures()
    call comparisons()

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message
        if (.not. condition) then
            print *, 'FAIL: ', message
            error stop 1
        end if
    end subroutine require

    subroutine orders()
        type(problem_t) :: p
        type(exprb_options_t) :: o
        type(exprb_state_t) :: s
        type(exprb_report_t) :: r
        type(fortnum_status_t) :: status
        real(dp) :: errors(4), order, h, target
        integer :: mode, method, k, n

        o%adaptive = .false.
        do mode = 1, 2
            p%mode = mode
            n = 1
            if (mode == 2) n = 16
            target = 0.5_dp
            do method = EXPRB_EULER, EXPRB32
                o%method = method
                do k = 1, 4
                    h = target/(10*2**k)
                    o%initial_step = h
                    o%max_step = h
                    call exprb_init(s, 0.0_dp, initial(mode, n), o, status)
                    call exprb_integrate_to(p, s, target, o, r, status)
                    call require(status%code == FORTNUM_OK, 'fixed-step order solve')
                    if (mode == 1) errors(k) = abs(s%y(1) - 1/(1 + target))
                    if (mode == 2) errors(k) = maxval(abs(s%y - exact_rd(n, target)))
                end do
                order = log(errors(3)/errors(4))/log(2.0_dp)
                print '(a,2i3,5es13.4)', 'order mode/method/errors/slope ', &
                    mode, method, errors, order
                call require(order > method - 0.2_dp, 'manufactured order lower bound')
                call require(order < method + 0.3_dp, 'manufactured order upper bound')
            end do
        end do
    end subroutine orders

    function vdp_reference(steps) result(y)
        integer, intent(in) :: steps
        real(dp) :: y(2), k1(2), k2(2), k3(2), k4(2), h
        integer :: i

        y = initial(3, 2)
        h = 100.0_dp/steps
        do i = 1, steps
            call rhs(3, 1.0_dp, 100.0_dp, y, k1)
            call rhs(3, 1.0_dp, 100.0_dp, y + h*k1/2, k2)
            call rhs(3, 1.0_dp, 100.0_dp, y + h*k2/2, k3)
            call rhs(3, 1.0_dp, 100.0_dp, y + h*k3, k4)
            y = y + h*(k1 + 2*k2 + 2*k3 + k4)/6
        end do
    end function vdp_reference

    subroutine tolerance_control()
        type(problem_t) :: p
        type(exprb_options_t) :: o
        type(exprb_state_t) :: s
        type(exprb_report_t) :: r
        type(fortnum_status_t) :: status
        real(dp) :: reference(2), errors(5), tolerance, reference_difference
        integer :: k, counts(5)

        reference = vdp_reference(1600000)
        reference_difference = maxval(abs(reference - vdp_reference(800000)))
        print '(a,es13.4)', 'RK4 refinement difference ', reference_difference
        call require(reference_difference < 1.0e-9_dp, &
            'independent RK4 reference refinement')
        p%mode = 3
        o%initial_step = 0.05_dp
        do k = 1, 5
            tolerance = 10.0_dp**(-k - 2)
            o%rtol = tolerance
            o%atol = tolerance*0.01_dp
            call exprb_init(s, 0.0_dp, initial(3, 2), o, status)
            call exprb_integrate_to(p, s, 100.0_dp, o, r, status)
            call require(status%code == FORTNUM_OK, 'van der Pol adaptive solve')
            errors(k) = maxval(abs(s%y - reference))
            counts(k) = r%accepted
            print '(a,3es13.4,3i8)', 'tolerance/error/ratio/steps/rejections ', &
                tolerance, errors(k), errors(k)/tolerance, r%accepted, &
                r%rejected_error, r%rejected_failure
            call require(errors(k) < 10*tolerance, 'van der Pol tolerance tracking')
        end do
        call require(errors(5) < errors(1)/100, 'four-decade tolerance improvement')
        call require(counts(5) > counts(1), 'controller responds to tolerance')
    end subroutine tolerance_control

    subroutine diffusion_refinement()
        type(problem_t) :: p
        type(exprb_options_t) :: o
        type(exprb_state_t) :: s
        type(exprb_report_t) :: r
        type(fortnum_status_t) :: status
        integer :: n, k, counts(3)
        real(dp) :: error

        p%mode = 4
        p%diffusion = 1.0e8_dp/(4*16.0_dp**2)
        o%initial_step = 0.01_dp
        o%rtol = 1.0e-5_dp
        o%atol = 1.0e-8_dp
        o%phi%rtol = 1.0e-9_dp
        o%phi%max_dimension = 70
        do k = 1, 3
            n = 8*2**k
            call exprb_init(s, 0.0_dp, initial(4, n), o, status)
            call exprb_integrate_to(p, s, 1.0_dp, o, r, status)
            call require(status%code == FORTNUM_OK, 'stiff diffusion solve')
            error = maxval(abs(s%y - exp(-1.0_dp)))
            counts(k) = r%accepted
            print '(a,i4,es13.4,5i9)', &
                'diffusion n/error/steps/J/dim/substeps/rejects ', &
                n, error, r%accepted, r%jacobian_actions, r%max_krylov_dimension, &
                r%phi_substeps, r%phi_rejections
            call require(error < 2.0e-6_dp, 'diffusion exact Fourier oracle')
            call require(r%accepted < 40, 'diffusion bounded step count')
            call require(r%jacobian_actions == p%actions, 'Jacobian action accounting')
            p%actions = 0
        end do
        call require(maxval(counts) <= 2*minval(counts), 'refinement step independence')
    end subroutine diffusion_refinement

    subroutine rejection_reuse()
        type(problem_t) :: p
        type(exprb_options_t) :: o
        type(exprb_state_t) :: s
        type(exprb_report_t) :: r
        type(fortnum_status_t) :: status
        integer :: method, evaluations, linearizations

        o%initial_step = 0.5_dp
        o%rtol = 1.0e-7_dp
        o%atol = 1.0e-9_dp
        do method = EXPRB_EULER, EXPRB32
            p = problem_t()
            o%method = method
            call exprb_init(s, 0.0_dp, [1.0_dp], o, status)
            call exprb_integrate_to(p, s, 0.5_dp, o, r, status)
            call require(status%code == FORTNUM_OK, 'rejected-step solve')
            call require(r%rejected_error > 0, 'large initial step forces rejection')
            call require(r%rejected_failure == 0, 'only error rejections')
            call require(abs(s%y(1) - 1/1.5_dp) < 5.0e-5_dp, 'reuse exact oracle')
            call require(p%linearizations == r%accepted, &
                'one linearization per accepted step despite rejections')
            call require(p%evaluations == 2*r%accepted + r%rejected_error, &
                'rejected steps reuse base evaluation')
            call require(r%linearizations == p%linearizations, &
                'rejection linearization accounting')
            call require(r%evaluations == p%evaluations, &
                'rejection evaluation accounting')
            print '(a,4i8)', 'reuse method/steps/rejections/linearizations ', &
                method, r%accepted, r%rejected_error, r%linearizations
            evaluations = p%evaluations
            linearizations = p%linearizations
            call exprb_integrate_to(p, s, 0.6_dp, o, r, status)
            call require(status%code == FORTNUM_OK, 'reuse continuation')
            call require(abs(s%y(1) - 1/1.6_dp) < 5.0e-5_dp, &
                'reuse continuation exact oracle')
            call require(p%linearizations - linearizations == r%accepted, &
                'continuation refreshes accepted endpoint')
            call require(r%linearizations == p%linearizations - linearizations, &
                'continuation linearization accounting')
            call require(r%evaluations == p%evaluations - evaluations, &
                'continuation evaluation accounting')
        end do
    end subroutine rejection_reuse

    subroutine failures()
        type(problem_t) :: p
        type(exprb_options_t) :: o
        type(exprb_state_t) :: s
        type(exprb_report_t) :: r
        type(fortnum_status_t) :: status
        integer :: kind

        o%initial_step = 0.05_dp
        do kind = 1, 4
            p = problem_t()
            if (kind == 1) p%fail_rhs = 1
            if (kind == 2) p%fail_setup = 1
            if (kind == 3) p%fail_action = 1
            if (kind == 4) p%fail_at = 2
            call exprb_init(s, 0.0_dp, [1.0_dp], o, status)
            call exprb_integrate_to(p, s, 0.2_dp, o, r, status)
            call require(status%code == FORTNUM_OK, 'retry recovers')
            call require(r%rejected_failure == 1, 'one callback/phi rejection')
            if (kind == 4) &
                call require(r%callback_status == 19, 'callback status retained')
            call require(abs(s%y(1) - 1/1.2_dp) < 1.0e-5_dp, 'retry exact oracle')
            call require(r%jacobian_actions == p%actions, 'retry action accounting')
            call require(r%evaluations == p%evaluations, 'evaluation accounting')
            call require(r%linearizations == p%linearizations, &
                'linearization accounting')
            if (kind == 2) then
                call require(p%linearizations == r%accepted + 1, &
                    'retry failed linearization only')
            else
                call require(p%linearizations == r%accepted, &
                    'reuse linearization after RHS or phi failure')
            end if
            if (kind == 1 .or. kind == 4) then
                call require(p%evaluations == 2*r%accepted + r%rejected_error + 1, &
                    'count failed RHS evaluation')
            else
                call require(p%evaluations == 2*r%accepted + r%rejected_error, &
                    'reuse base evaluation after setup or phi failure')
            end if
        end do
        call exprb_integrate_to(p, s, 0.4_dp, o, r, status)
        call require(status%code == FORTNUM_OK, 'continuation')
        call require(abs(s%y(1) - 1/1.4_dp) < 2.0e-5_dp, 'continuation oracle')
        call exprb_integrate_to(p, s, 0.0_dp, o, r, status)
        call require(status%code == FORTNUM_OK, 'reverse integration')
        call require(abs(s%y(1) - 1) < 3.0e-5_dp, 'reverse exact oracle')
        call exprb_integrate_to(p, s, s%t, o, r, status)
        call require(r%accepted == 0, 'equal endpoint no-op')
        o%method = EXPRB_EULER
        call exprb_init(s, 0.0_dp, [1.0_dp], o, status)
        call exprb_integrate_to(p, s, 0.2_dp, o, r, status)
        call require(status%code == FORTNUM_OK, 'adaptive Euler')
        call require(abs(s%y(1) - 1/1.2_dp) < 5.0e-4_dp, 'adaptive Euler oracle')
        o%max_steps = 1
        p%fail_rhs = 3
        call exprb_init(s, 0.0_dp, [1.0_dp], o, status)
        call exprb_integrate_to(p, s, 0.2_dp, o, r, status)
        call require(status%code == FORTNUM_CONVERGENCE_ERROR, 'attempt exhaustion')
        call require(s%t == 0 .and. s%y(1) == 1, 'rejection preserves accepted state')
        o%max_steps = 100
        o%min_step = 0.1_dp
        call exprb_integrate_to(p, s, 0.02_dp, o, r, status)
        call require(status%code == FORTNUM_CONVERGENCE_ERROR, 'minimum step')
        o%rtol = -1.0_dp
        call exprb_init(s, 0.0_dp, [1.0_dp], o, status)
        call require(status%code == FORTNUM_DOMAIN_ERROR, 'invalid tolerance')
        o = exprb_options_t()
        o%phi%max_dimension = 0
        call exprb_init(s, 0.0_dp, [1.0_dp], o, status)
        call require(status%code == FORTNUM_DOMAIN_ERROR, 'invalid phi controls')
    end subroutine failures

    subroutine comparisons()
        real(dp), allocatable :: y0(:), reference(:)
        real(dp), parameter :: error_target = 1.0e-6_dp
        real(dp) :: target, error, tolerance, power, samples(3), median
        integer :: mode, method, n, calibration, repeat, steps, evaluations, actions
        integer(int64) :: clock_start, clock_end, rate

        do mode = 1, 4
            n = 1
            target = 0.5_dp
            if (mode == 2) n = 16
            if (mode == 3) then
                n = 2
                target = 100.0_dp
            end if
            if (mode == 4) then
                n = 16
                target = 1.0_dp
            end if
            y0 = initial(mode, n)
            reference = y0
            if (mode == 1) reference = 1/(1 + target)
            if (mode == 2) reference = exact_rd(n, target)
            if (mode == 3) reference = vdp_reference(1600000)
            if (mode == 4) reference = exp(-target)
            do method = 2, 4
                tolerance = 1.0e-5_dp
                power = 1.5_dp
                if (method == EXPRB32) power = 1.0_dp
                ! Calibration is excluded from the measured samples.
                do calibration = 1, 10
                    call comparison_solve(mode, method, tolerance, y0, &
                        reference, target, error, steps, evaluations, actions)
                    if (error <= error_target) then
                        if (error >= 0.1_dp*error_target .or. mode == 4) exit
                    end if
                    tolerance = tolerance*max(0.02_dp, min(10.0_dp, &
                        (0.6_dp*error_target/max(error, tiny(error)))**power))
                end do
                call require(error <= error_target, 'matched-error comparison')
                if (mode /= 4) call require(error >= 0.1_dp*error_target, &
                    'comparison errors within one decade')
                do repeat = 1, 3
                    call system_clock(clock_start, rate)
                    call comparison_solve(mode, method, tolerance, y0, &
                        reference, target, error, steps, evaluations, actions)
                    call system_clock(clock_end)
                    samples(repeat) = real(clock_end - clock_start, dp)/rate
                end do
                median = sum(samples) - minval(samples) - maxval(samples)
                print '(a,2i3,2es13.4,3i9,3es13.4)', &
                    'compare mode/method/tol/error/steps/f/J/'// &
                    'seconds(median,min,max) ', &
                    mode, method, tolerance, error, steps, evaluations, actions, &
                    median, minval(samples), maxval(samples)
            end do
        end do

    end subroutine comparisons

    subroutine comparison_solve(mode, method, tolerance, y0, reference, target, &
            error, steps, evaluations, actions)
        integer, intent(in) :: mode, method
        real(dp), intent(in) :: tolerance, y0(:), reference(:), target
        real(dp), intent(out) :: error
        integer, intent(out) :: steps, evaluations, actions
        type(problem_t) :: p
        type(implicit_problem_t) :: b
        type(exprb_options_t) :: o
        type(bdf_options_t) :: bo
        type(exprb_state_t) :: s
        type(bdf_state_t) :: bs
        type(exprb_report_t) :: r
        type(bdf_report_t) :: br
        type(fortnum_status_t) :: status
        real(dp) :: derivative(size(y0))

        if (method <= 3) then
            p = problem_t()
            p%mode = mode
            if (mode == 4) p%diffusion = 1.0e8_dp/(4*16.0_dp**2)
            o = exprb_options_t()
            o%method = method
            o%rtol = tolerance
            o%atol = tolerance*0.01_dp
            o%initial_step = 0.001_dp
            o%phi%rtol = 1.0e-9_dp
            call exprb_init(s, 0.0_dp, y0, o, status)
            call exprb_integrate_to(p, s, target, o, r, status)
            error = maxval(abs(s%y - reference))
            steps = r%accepted
            evaluations = r%evaluations
            actions = r%jacobian_actions
        else
            b = implicit_problem_t()
            b%mode = mode
            if (mode == 4) b%diffusion = 1.0e8_dp/(4*16.0_dp**2)
            bo = bdf_options_t()
            bo%rtol = tolerance
            bo%atol = tolerance*0.01_dp
            bo%initial_step = 1.0e-5_dp
            derivative = y0
            call rhs(mode, b%diffusion, b%mu, y0, derivative)
            b%evaluations = 1
            call bdf_init(bs, 0.0_dp, y0, bo, status, derivative)
            call bdf_integrate_to(b, bs, target, bo, br, status)
            error = maxval(abs(bs%y - reference))
            steps = br%accepted
            evaluations = b%evaluations
            actions = b%actions
        end if
        call require(status%code == FORTNUM_OK, 'comparison solve')
    end subroutine comparison_solve
end program test_fortnum_ode_exprb
