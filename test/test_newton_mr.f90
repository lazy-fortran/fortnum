module newton_mr_fixtures
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_newton_mr, only: newton_mr_problem_t
    implicit none
    private

    public :: rosenbrock_t, lagrangian_t, double_well_t

    type, extends(newton_mr_problem_t) :: rosenbrock_t
        !! Extended Rosenbrock: sum over pairs of 100 (x2 - x1^2)^2 + (1 - x1)^2,
        !! minimizer x = 1, indefinite Hessian away from it.
        real(dp), allocatable :: h(:, :)
        logical :: jacobi = .false., fail = .false.
    contains
        procedure :: evaluate => rosenbrock_evaluate
        procedure :: linearize => rosenbrock_linearize
        procedure :: apply => rosenbrock_apply
        procedure :: precondition => rosenbrock_precondition
    end type rosenbrock_t

    type, extends(newton_mr_problem_t) :: double_well_t
        !! sum(x^4/4 - x^2/2): minimizers x_i = +-1, Hessian 3 x^2 - 1 < 0 for
        !! |x| < 1/sqrt(3).
        real(dp), allocatable :: h(:)
    contains
        procedure :: evaluate => double_well_evaluate
        procedure :: linearize => double_well_linearize
        procedure :: apply => double_well_apply
    end type double_well_t

    type, extends(newton_mr_problem_t) :: lagrangian_t
        !! L(x, l) = sum(x^4/4 + x^2/2) - c.x + l (|x|^2 - 1)/2: the
        !! stationary point is a saddle (bordered Hessian).
        real(dp), allocatable :: c(:), d(:), xl(:)
        real(dp) :: schur = 1.0_dp
    contains
        procedure :: evaluate => lagrangian_evaluate
        procedure :: linearize => lagrangian_linearize
        procedure :: apply => lagrangian_apply
        procedure :: precondition => lagrangian_precondition
    end type lagrangian_t

contains

    subroutine rosenbrock_evaluate(self, x, f, g, status)
        class(rosenbrock_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f, g(:)
        integer, intent(out) :: status
        integer :: i

        f = 0.0_dp
        do i = 1, size(x), 2
            f = f + 100.0_dp*(x(i + 1) - x(i)**2)**2 + (1.0_dp - x(i))**2
            g(i) = -400.0_dp*x(i)*(x(i + 1) - x(i)**2) - 2.0_dp*(1.0_dp - x(i))
            g(i + 1) = 200.0_dp*(x(i + 1) - x(i)**2)
        end do
        status = 0
    end subroutine rosenbrock_evaluate

    subroutine rosenbrock_linearize(self, x, status)
        class(rosenbrock_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: status
        integer :: i

        if (.not. allocated(self%h)) allocate (self%h(size(x), size(x)))
        self%h = 0.0_dp
        do i = 1, size(x), 2
            self%h(i, i) = 1200.0_dp*x(i)**2 - 400.0_dp*x(i + 1) + 2.0_dp
            self%h(i, i + 1) = -400.0_dp*x(i)
            self%h(i + 1, i) = -400.0_dp*x(i)
            self%h(i + 1, i + 1) = 200.0_dp
        end do
        status = 0
    end subroutine rosenbrock_linearize

    subroutine rosenbrock_apply(self, input, output)
        class(rosenbrock_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = matmul(self%h, input)
    end subroutine rosenbrock_apply

    subroutine rosenbrock_precondition(self, r, z, status)
        class(rosenbrock_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        integer :: i

        status = 0
        if (self%fail) then
            z = 0.0_dp
            status = 1
            return
        end if
        if (.not. self%jacobi) then
            z = r
            return
        end if
        do i = 1, size(r)
            z(i) = r(i)/max(abs(self%h(i, i)), 1.0_dp)
        end do
    end subroutine rosenbrock_precondition

    subroutine double_well_evaluate(self, x, f, g, status)
        class(double_well_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f, g(:)
        integer, intent(out) :: status

        f = sum(x**4/4 - x**2/2)
        g = x**3 - x
        status = 0
    end subroutine double_well_evaluate

    subroutine double_well_linearize(self, x, status)
        class(double_well_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: status

        self%h = 3*x**2 - 1
        status = 0
    end subroutine double_well_linearize

    subroutine double_well_apply(self, input, output)
        class(double_well_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)

        output = self%h*input
    end subroutine double_well_apply

    subroutine lagrangian_evaluate(self, x, f, g, status)
        class(lagrangian_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f, g(:)
        integer, intent(out) :: status
        integer :: n

        n = size(self%c)
        associate (y => x(1:n), l => x(n + 1))
            f = sum(y**4/4 + y**2/2) - dot_product(self%c, y) &
                + l*(dot_product(y, y) - 1.0_dp)/2
            g(1:n) = y**3 + y - self%c + l*y
            g(n + 1) = (dot_product(y, y) - 1.0_dp)/2
        end associate
        status = 0
    end subroutine lagrangian_evaluate

    subroutine lagrangian_linearize(self, x, status)
        class(lagrangian_t), intent(inout) :: self
        real(dp), intent(in) :: x(:)
        integer, intent(out) :: status
        integer :: n

        n = size(self%c)
        self%xl = x
        self%d = 3*x(1:n)**2 + 1 + x(n + 1)
        ! SPD block preconditioner: |K|^-1 and the Schur approximation.
        self%schur = 1.0_dp/sum(x(1:n)**2/abs(self%d))
        status = 0
    end subroutine lagrangian_linearize

    subroutine lagrangian_apply(self, input, output)
        class(lagrangian_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        integer :: n

        n = size(self%c)
        output(1:n) = self%d*input(1:n) + self%xl(1:n)*input(n + 1)
        output(n + 1) = dot_product(self%xl(1:n), input(1:n))
    end subroutine lagrangian_apply

    subroutine lagrangian_precondition(self, r, z, status)
        class(lagrangian_t), intent(inout) :: self
        real(dp), intent(in) :: r(:)
        real(dp), intent(out) :: z(:)
        integer, intent(out) :: status
        integer :: n

        n = size(self%c)
        z(1:n) = r(1:n)/abs(self%d)
        z(n + 1) = self%schur*r(n + 1)
        status = 0
    end subroutine lagrangian_precondition
end module newton_mr_fixtures

program test_newton_mr
    !! Newton-MR against independent oracles: the known minimizer of the
    !! extended Rosenbrock function and of a double well from a start of
    !! negative curvature (objective merit, curvature exits and forward
    !! tracking taken), and the saddle point of a
    !! constrained Lagrangian from a scalar bisection on the multiplier
    !! (preconditioned gradient merit).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_newton_mr, only: newton_mr_solve, newton_mr_options_t, &
        newton_mr_report_t, NEWTON_MR_OK, NEWTON_MR_BREAKDOWN, &
        NEWTON_MR_INVALID_ARGUMENT, NEWTON_MR_MERIT_GRADIENT
    use newton_mr_fixtures, only: rosenbrock_t, lagrangian_t, double_well_t
    implicit none

    integer, parameter :: n = 10, m = 8
    type(rosenbrock_t) :: rosenbrock
    type(lagrangian_t) :: lagrangian
    type(double_well_t) :: well
    type(newton_mr_options_t) :: options
    type(newton_mr_report_t) :: report
    real(dp) :: x(n), y(m + 1), exact(m), lambda
    integer :: i

    ! 1. Rosenbrock from (-1.2, 1) pairs.
    options%tolerance = 1.0e-12_dp
    options%max_iterations = 200
    x(1:n:2) = -1.2_dp
    x(2:n:2) = 1.0_dp
    call newton_mr_solve(rosenbrock, x, options, report)
    print '(a,4i6)', 'rosenbrock: iterations, minres, npc, evaluations ', &
        report%iterations, report%minres_iterations, report%npc_steps, &
        report%evaluations
    call require(report%status == NEWTON_MR_OK, "Newton-MR converges on Rosenbrock")
    call require(maxval(abs(x - 1.0_dp)) < 1.0e-8_dp, &
        "Newton-MR finds the Rosenbrock minimizer")
    call require(report%gradient <= 1.0e-12_dp*report%initial_gradient, &
        "reported gradient meets the tolerance")

    rosenbrock%jacobi = .true.
    x(1:n:2) = -1.2_dp
    x(2:n:2) = 1.0_dp
    call newton_mr_solve(rosenbrock, x, options, report)
    call require(report%status == NEWTON_MR_OK .and. &
        maxval(abs(x - 1.0_dp)) < 1.0e-8_dp, &
        "preconditioned Newton-MR finds the Rosenbrock minimizer")

    rosenbrock%fail = .true.
    x = 0.0_dp
    call newton_mr_solve(rosenbrock, x, options, report)
    call require(report%status == NEWTON_MR_BREAKDOWN, &
        "a failed preconditioner is a breakdown")
    rosenbrock%fail = .false.

    options%armijo = 0.0_dp
    call newton_mr_solve(rosenbrock, x, options, report)
    call require(report%status == NEWTON_MR_INVALID_ARGUMENT, &
        "a zero Armijo constant is rejected")
    options = newton_mr_options_t()

    ! 2. Double well from 0 < x < 1/sqrt(3): negative curvature everywhere.
    do i = 1, n
        x(i) = 0.05_dp*i
    end do
    options%tolerance = 1.0e-12_dp
    call newton_mr_solve(well, x, options, report)
    print '(a,4i6)', 'double well: iterations, minres, npc, evaluations ', &
        report%iterations, report%minres_iterations, report%npc_steps, &
        report%evaluations
    call require(report%status == NEWTON_MR_OK .and. &
        maxval(abs(x - 1.0_dp)) < 1.0e-10_dp, &
        "Newton-MR leaves the negative-curvature region to the minimizer")
    call require(report%npc_steps >= 1, "the concave start takes a curvature step")
    options = newton_mr_options_t()

    ! 3. Lagrangian saddle with the preconditioned gradient merit.
    allocate (lagrangian%c(m))
    do i = 1, m
        lagrangian%c(i) = 1.0_dp + 0.3_dp*i
    end do
    y(1:m) = lagrangian%c/norm2(lagrangian%c)
    y(m + 1) = 0.0_dp
    options%merit = NEWTON_MR_MERIT_GRADIENT
    options%tolerance = 1.0e-12_dp
    options%max_iterations = 100
    call newton_mr_solve(lagrangian, y, options, report)
    print '(a,4i6)', 'lagrangian: iterations, minres, npc, evaluations ', &
        report%iterations, report%minres_iterations, report%npc_steps, &
        report%evaluations
    call require(report%status == NEWTON_MR_OK, "Newton-MR converges on the saddle")
    call require(report%npc_steps == 0, "the gradient merit takes no curvature exit")
    call multiplier_oracle(lagrangian%c, exact, lambda)
    call require(maxval(abs(y(1:m) - exact)) < 1.0e-10_dp .and. &
        abs(y(m + 1) - lambda) < 1.0e-9_dp*max(1.0_dp, abs(lambda)), &
        "Newton-MR finds the bisection saddle point")

    print '(a)', 'PASS'

contains

    subroutine multiplier_oracle(c, xs, l)
        !! x_i(l) solves x^3 + (1 + l) x = c_i (monotone, unique for l > -1);
        !! |x(l)|^2 decreases in l, so bisection finds |x|^2 = 1.
        real(dp), intent(in) :: c(:)
        real(dp), intent(out) :: xs(:), l
        real(dp) :: lo, hi
        integer :: it

        lo = -1.0_dp
        hi = 1.0e3_dp
        do it = 1, 200
            l = 0.5_dp*(lo + hi)
            call components(c, l, xs)
            if (dot_product(xs, xs) > 1.0_dp) then
                lo = l
            else
                hi = l
            end if
        end do
        l = 0.5_dp*(lo + hi)
        call components(c, l, xs)
    end subroutine multiplier_oracle

    subroutine components(c, l, xs)
        real(dp), intent(in) :: c(:), l
        real(dp), intent(out) :: xs(:)
        real(dp) :: lo, hi, t
        integer :: i, it

        do i = 1, size(c)
            lo = 0.0_dp
            hi = c(i)
            do it = 1, 200
                t = 0.5_dp*(lo + hi)
                if (t**3 + (1.0_dp + l)*t > c(i)) then
                    hi = t
                else
                    lo = t
                end if
            end do
            xs(i) = 0.5_dp*(lo + hi)
        end do
    end subroutine components

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            write (error_unit, '(a)') "FAIL: "//message
            error stop 1
        end if
    end subroutine require
end program test_newton_mr
