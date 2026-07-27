module vector_root_vjp_hybrid_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use fortnum_kinds, only: dp
    use fortnum_multiroot, only: multiroot_implicit_vjp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    type, bind(c) :: residual_gradient_t
        real(c_double) :: values(4)
    end type residual_gradient_t

    public :: analytical_root_vjp, hybrid_root_vjp
    public :: autodiff_iterations_root_vjp, finite_difference_root_vjp
    public :: residual_one, residual_two

    interface
        function enzyme_autodiff(f, x1, x2, p1, p2) result(gradient) &
                bind(c, name="__enzyme_autodiff")
            import :: c_double, c_funptr, residual_gradient_t
            type(c_funptr), value :: f
            real(c_double), value :: x1, x2, p1, p2
            type(residual_gradient_t) :: gradient
        end function enzyme_autodiff
    end interface

contains

    pure function residual_one(x1, x2, p1, p2) result(f) &
            bind(c, name="fortnum_vector_root_vjp_residual_one")
        real(c_double), value :: x1, x2, p1, p2
        real(c_double) :: f

        f = x1*x1 + x2 - p1 + 0.0_c_double*p2
    end function residual_one

    pure function residual_two(x1, x2, p1, p2) result(f) &
            bind(c, name="fortnum_vector_root_vjp_residual_two")
        real(c_double), value :: x1, x2, p1, p2
        real(c_double) :: f

        f = x1 + x2*x2 + 0.0_c_double*p1 - p2
    end function residual_two

    pure function newton_objective(p1, p2, u1, u2) result(objective) &
            bind(c, name="fortnum_vector_newton_objective")
        real(c_double), value :: p1, p2, u1, u2
        real(c_double) :: objective
        real(c_double) :: x1, x2

        call fixed_newton_solve(p1, p2, x1, x2)
        objective = u1*x1 + u2*x2
    end function newton_objective

    pure subroutine fixed_newton_solve(p1, p2, x1, x2)
        real(c_double), intent(in) :: p1, p2
        real(c_double), intent(out) :: x1, x2
        real(c_double) :: determinant, f1, f2, step1, step2
        integer :: iteration

        x1 = 0.8_c_double
        x2 = 1.1_c_double
        do iteration = 1, 12
            f1 = residual_one(x1, x2, p1, p2)
            f2 = residual_two(x1, x2, p1, p2)
            determinant = 4.0_c_double*x1*x2 - 1.0_c_double
            step1 = (-2.0_c_double*x2*f1 + f2)/determinant
            step2 = (-2.0_c_double*x1*f2 + f1)/determinant
            x1 = x1 + step1
            x2 = x2 + step2
        end do
    end subroutine fixed_newton_solve

    function analytical_root_vjp(x, p, u) result(jtu)
        real(dp), intent(in) :: x(2), p(2), u(2)
        real(dp) :: jtu(2)
        type(fortnum_status_t) :: status

        call multiroot_implicit_vjp(analytical_state_jacobian, &
            analytical_parameter_vjp, x, p, u, jtu, status)
        if (.not. status_ok(status)) error stop "analytical vector-root VJP failed"
    end function analytical_root_vjp

    function hybrid_root_vjp(x, p, u) result(jtu)
        real(dp), intent(in) :: x(2), p(2), u(2)
        real(dp) :: jtu(2)
        type(fortnum_status_t) :: status

        call multiroot_implicit_vjp(autodiff_state_jacobian, &
            autodiff_parameter_vjp, x, p, u, jtu, status)
        if (.not. status_ok(status)) error stop "hybrid vector-root VJP failed"
    end function hybrid_root_vjp

    function autodiff_iterations_root_vjp(p, u) result(jtu)
        real(dp), intent(in) :: p(2), u(2)
        real(dp) :: jtu(2)
        type(residual_gradient_t) :: gradient

        gradient = enzyme_autodiff(c_funloc(newton_objective), &
            p(1), p(2), u(1), u(2))
        jtu = gradient%values(1:2)
    end function autodiff_iterations_root_vjp

    function finite_difference_root_vjp(p, u) result(jtu)
        real(dp), intent(in) :: p(2), u(2)
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: jtu(2), pm(2), pp(2)
        integer :: parameter

        do parameter = 1, 2
            pp = p
            pm = p
            pp(parameter) = pp(parameter) + h
            pm(parameter) = pm(parameter) - h
            jtu(parameter) = (newton_objective( &
                pp(1), pp(2), u(1), u(2)) - newton_objective( &
                pm(1), pm(2), u(1), u(2)))/(2.0_dp*h)
        end do
    end function finite_difference_root_vjp

    subroutine analytical_state_jacobian(x, p, jac_x, context)
        real(dp), intent(in) :: x(:), p(:)
        real(dp), intent(out) :: jac_x(size(x), size(x))
        class(*), intent(inout), optional :: context

        if (size(p) /= 2) error stop "vector residual expects two parameters"
        jac_x(1, 1) = 2.0_dp*x(1)
        jac_x(1, 2) = 1.0_dp
        jac_x(2, 1) = 1.0_dp
        jac_x(2, 2) = 2.0_dp*x(2)
    end subroutine analytical_state_jacobian

    subroutine analytical_parameter_vjp(x, p, u, f_p_t_u, context)
        real(dp), intent(in) :: x(:), p(:), u(:)
        real(dp), intent(out) :: f_p_t_u(size(p))
        class(*), intent(inout), optional :: context

        if (size(x) /= 2) error stop "vector residual expects two states"
        f_p_t_u = -u
    end subroutine analytical_parameter_vjp

    subroutine autodiff_state_jacobian(x, p, jac_x, context)
        real(dp), intent(in) :: x(:), p(:)
        real(dp), intent(out) :: jac_x(size(x), size(x))
        class(*), intent(inout), optional :: context
        type(residual_gradient_t) :: gradient_one, gradient_two

        gradient_one = enzyme_autodiff(c_funloc(residual_one), &
            x(1), x(2), p(1), p(2))
        gradient_two = enzyme_autodiff(c_funloc(residual_two), &
            x(1), x(2), p(1), p(2))
        jac_x(1, :) = gradient_one%values(1:2)
        jac_x(2, :) = gradient_two%values(1:2)
    end subroutine autodiff_state_jacobian

    subroutine autodiff_parameter_vjp(x, p, u, f_p_t_u, context)
        real(dp), intent(in) :: x(:), p(:), u(:)
        real(dp), intent(out) :: f_p_t_u(size(p))
        class(*), intent(inout), optional :: context
        type(residual_gradient_t) :: gradient_one, gradient_two

        gradient_one = enzyme_autodiff(c_funloc(residual_one), &
            x(1), x(2), p(1), p(2))
        gradient_two = enzyme_autodiff(c_funloc(residual_two), &
            x(1), x(2), p(1), p(2))
        f_p_t_u = u(1)*gradient_one%values(3:4) &
            + u(2)*gradient_two%values(3:4)
    end subroutine autodiff_parameter_vjp

end module vector_root_vjp_hybrid_kernel

program enzyme_vector_root_vjp_hybrid
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use vector_root_vjp_hybrid_kernel, only: analytical_root_vjp, &
        hybrid_root_vjp, autodiff_iterations_root_vjp, &
        finite_difference_root_vjp
    implicit none

    interface
        function peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function peak_rss_bytes
    end interface

    real(dp), parameter :: h = 1.0e-5_dp
    real(dp) :: p(2), pp(2), pm(2), u(2), xstar(2), xp(2), xm(2)
    real(dp) :: analytical(2), got(2), reference(2)
    real(dp) :: autodiff_iterations(2), diagnostic(2), errors(4)
    character(32) :: argument, candidate
    integer :: parameter

    call get_command_argument(1, argument)
    call get_command_argument(2, candidate)
    if (trim(argument) == "--benchmark") then
        call run_benchmark()
    else if (trim(argument) == "--peak-rss") then
        call run_peak_rss(trim(candidate))
    else
        p = [2.0_dp, 2.0_dp]
        u = [1.3_dp, -0.4_dp]
        xstar = solve_root(p)
        do parameter = 1, size(p)
            pp = p
            pm = p
            pp(parameter) = pp(parameter) + h
            pm(parameter) = pm(parameter) - h
            xp = solve_root(pp)
            xm = solve_root(pm)
            reference(parameter) = dot_product(u, xp - xm)/(2.0_dp*h)
        end do
        analytical = analytical_root_vjp(xstar, p, u)
        got = hybrid_root_vjp(xstar, p, u)
        autodiff_iterations = autodiff_iterations_root_vjp(p, u)
        diagnostic = finite_difference_root_vjp(p, u)
        errors(1) = maxval(abs(analytical - reference))
        errors(2) = maxval(abs(got - reference))
        errors(3) = maxval(abs(autodiff_iterations - reference))
        errors(4) = maxval(abs(diagnostic - reference))

        if (maxval(errors) > 1.0e-9_dp) then
            print *, "vector-root VJP mismatch", errors
            error stop 1
        end if
        print *, "PASS vector-root VJP tournament errors", errors
    end if

contains

    pure function solve_root(parameters) result(x)
        real(dp), intent(in) :: parameters(2)
        real(dp) :: x(2), f(2), jacobian(2, 2), step(2), determinant
        integer :: iteration

        x = [0.8_dp, 1.1_dp]
        do iteration = 1, 12
            f(1) = x(1)*x(1) + x(2) - parameters(1)
            f(2) = x(1) + x(2)*x(2) - parameters(2)
            jacobian(1, 1) = 2.0_dp*x(1)
            jacobian(1, 2) = 1.0_dp
            jacobian(2, 1) = 1.0_dp
            jacobian(2, 2) = 2.0_dp*x(2)
            determinant = jacobian(1, 1)*jacobian(2, 2) &
                - jacobian(1, 2)*jacobian(2, 1)
            step(1) = (-f(1)*jacobian(2, 2) + jacobian(1, 2)*f(2)) &
                /determinant
            step(2) = (-jacobian(1, 1)*f(2) + f(1)*jacobian(2, 1)) &
                /determinant
            x = x + step
        end do
    end function solve_root

    subroutine run_benchmark()
        integer, parameter :: samples = 15
        integer(int64), parameter :: reps = 100000_int64
        real(dp) :: analytical_times(samples), hybrid_times(samples)
        real(dp) :: autodiff_times(samples), diagnostic_times(samples), sink
        integer :: sample

        do sample = 1, 3
            call time_candidate("analytical", reps/100_int64, sink)
            call time_candidate("hybrid", reps/100_int64, sink)
            call time_candidate("autodiff", reps/100_int64, sink)
            call time_candidate("diagnostic", reps/100_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate("analytical", reps, analytical_times(sample))
            call time_candidate("hybrid", reps, hybrid_times(sample))
            call time_candidate("autodiff", reps, autodiff_times(sample))
            call time_candidate("diagnostic", reps, diagnostic_times(sample))
        end do
        call report("analytical", analytical_times, reps)
        call report("hybrid", hybrid_times, reps)
        call report("autodiff", autodiff_times, reps)
        call report("diagnostic", diagnostic_times, reps)
    end subroutine run_benchmark

    subroutine run_peak_rss(name)
        character(*), intent(in) :: name
        integer(int64), parameter :: reps = 100000_int64
        real(dp) :: elapsed

        if ((name /= "analytical") .and. (name /= "hybrid") .and. &
                (name /= "autodiff") .and. (name /= "diagnostic")) then
            error stop "usage: --peak-rss analytical|hybrid|autodiff|diagnostic"
        end if
        call time_candidate(name, reps, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine time_candidate(name, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        integer(int64) :: i, start, finish, rate
        real(dp) :: parameters(2), cotangent(2), x(2), jtu(2), sink

        sink = 0.0_dp
        call system_clock(start, rate)
        do i = 1_int64, reps
            x(1) = 1.05_dp + 0.001_dp*real(mod(i, 101_int64), dp)
            x(2) = 0.85_dp + 0.001_dp*real(mod(i, 79_int64), dp)
            parameters(1) = x(1)*x(1) + x(2)
            parameters(2) = x(1) + x(2)*x(2)
            cotangent(1) = 0.02_dp*real(mod(i, 17_int64) - 8_int64, dp)
            cotangent(2) = -0.6_dp
            select case (name)
            case ("analytical")
                jtu = analytical_root_vjp(x, parameters, cotangent)
            case ("hybrid")
                jtu = hybrid_root_vjp(x, parameters, cotangent)
            case ("autodiff")
                jtu = autodiff_iterations_root_vjp(parameters, cotangent)
            case ("diagnostic")
                jtu = finite_difference_root_vjp(parameters, cotangent)
            end select
            sink = sink + sum(jtu)
        end do
        call system_clock(finish)
        elapsed_ns = 1.0e9_dp*real(finish - start, dp) &
            / (real(rate, dp)*real(reps, dp))
        if (sink /= sink) error stop "benchmark failed"
    end subroutine time_candidate

    subroutine report(name, values, reps)
        character(*), intent(in) :: name
        real(dp), intent(in) :: values(:)
        integer(int64), intent(in) :: reps
        real(dp) :: ordered(size(values)), deviations(size(values))
        real(dp) :: median, mad

        ordered = values
        call sort_values(ordered)
        median = ordered((size(ordered) + 1)/2)
        deviations = abs(values - median)
        call sort_values(deviations)
        mad = deviations((size(deviations) + 1)/2)
        write (*, "(a,',',i0,',',f12.4,',',f12.4)") &
            name, reps, median, mad
    end subroutine report

    subroutine sort_values(values)
        real(dp), intent(inout) :: values(:)
        real(dp) :: temporary
        integer :: i, j

        do i = 1, size(values) - 1
            do j = i + 1, size(values)
                if (values(j) < values(i)) then
                    temporary = values(i)
                    values(i) = values(j)
                    values(j) = temporary
                end if
            end do
        end do
    end subroutine sort_values

end program enzyme_vector_root_vjp_hybrid
