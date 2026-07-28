module vector_root_vjp_hybrid_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_generated_enzyme_vector_root_objective, only: &
        fortnum_enzyme_vector_root_objective_vjp
    use fortnum_generated_enzyme_vector_root_vjp_residual_one, only: &
        fortnum_enzyme_vector_root_vjp_residual_one_vjp
    use fortnum_generated_enzyme_vector_root_vjp_residual_two, only: &
        fortnum_enzyme_vector_root_vjp_residual_two_vjp
    use fortnum_generated_vector_root_residual_jacobian, only: &
        fortnum_vector_root_residual_jacobian_kernel
    use fortnum_generated_vector_root_residual_vjp, only: &
        fortnum_vector_root_residual_vjp_kernel
    use fortnum_kinds, only: dp
    use fortnum_multiroot, only: multiroot_implicit_vjp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    public :: analytical_root_vjp, hybrid_root_vjp
    public :: autodiff_iterations_root_vjp, finite_difference_root_vjp
    public :: residual_one, residual_two

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
        real(dp) :: ignored_u1bar, ignored_u2bar

        call fortnum_enzyme_vector_root_objective_vjp( &
            p(1), p(2), u(1), u(2), 1.0_dp, jtu(1), jtu(2), &
            ignored_u1bar, ignored_u2bar)
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
        call fortnum_vector_root_residual_jacobian_kernel(x, jac_x)
    end subroutine analytical_state_jacobian

    subroutine analytical_parameter_vjp(x, p, u, f_p_t_u, context)
        real(dp), intent(in) :: x(:), p(:), u(:)
        real(dp), intent(out) :: f_p_t_u(size(p))
        class(*), intent(inout), optional :: context

        if (size(x) /= 2) error stop "vector residual expects two states"
        call fortnum_vector_root_residual_vjp_kernel(u, f_p_t_u)
    end subroutine analytical_parameter_vjp

    subroutine autodiff_state_jacobian(x, p, jac_x, context)
        real(dp), intent(in) :: x(:), p(:)
        real(dp), intent(out) :: jac_x(size(x), size(x))
        class(*), intent(inout), optional :: context
        real(dp) :: ignored_p1bar, ignored_p2bar

        call fortnum_enzyme_vector_root_vjp_residual_one_vjp( &
            x(1), x(2), p(1), p(2), 1.0_dp, jac_x(1, 1), &
            jac_x(1, 2), ignored_p1bar, ignored_p2bar)
        call fortnum_enzyme_vector_root_vjp_residual_two_vjp( &
            x(1), x(2), p(1), p(2), 1.0_dp, jac_x(2, 1), &
            jac_x(2, 2), ignored_p1bar, ignored_p2bar)
    end subroutine autodiff_state_jacobian

    subroutine autodiff_parameter_vjp(x, p, u, f_p_t_u, context)
        real(dp), intent(in) :: x(:), p(:), u(:)
        real(dp), intent(out) :: f_p_t_u(size(p))
        class(*), intent(inout), optional :: context
        real(dp) :: ignored_x1bar, ignored_x2bar
        real(dp) :: pbar_one(2), pbar_two(2)

        call fortnum_enzyme_vector_root_vjp_residual_one_vjp( &
            x(1), x(2), p(1), p(2), u(1), ignored_x1bar, &
            ignored_x2bar, pbar_one(1), pbar_one(2))
        call fortnum_enzyme_vector_root_vjp_residual_two_vjp( &
            x(1), x(2), p(1), p(2), u(2), ignored_x1bar, &
            ignored_x2bar, pbar_two(1), pbar_two(2))
        f_p_t_u = pbar_one + pbar_two
    end subroutine autodiff_parameter_vjp

end module vector_root_vjp_hybrid_kernel

program enzyme_vector_root_vjp_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, write_fixture_result
    use vector_root_vjp_hybrid_kernel, only: analytical_root_vjp, &
        hybrid_root_vjp, autodiff_iterations_root_vjp, &
        finite_difference_root_vjp
    implicit none

    integer, parameter :: repetitions = 100000
    real(dp), parameter :: h = 1.0e-5_dp
    real(dp) :: p(2), pp(2), pm(2), u(2), xstar(2), xp(2), xm(2)
    real(dp) :: analytical(2), got(2), reference(2)
    real(dp) :: autodiff_iterations(2), diagnostic(2), errors(4)
    real(dp) :: samples(fixture_sample_count), sink
    character(32) :: argument, candidate
    integer :: parameter, candidate_kind

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
        call benchmark_one("analytical", 1)
        call benchmark_one("hybrid", 2)
        call benchmark_one("autodiff", 3)
        call benchmark_one("diagnostic", 4)
    end subroutine run_benchmark

    subroutine run_peak_rss(name)
        character(*), intent(in) :: name

        call parse_candidate(name)
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine benchmark_one(name, kind)
        character(*), intent(in) :: name
        integer, intent(in) :: kind
        real(dp) :: median, mad

        candidate_kind = kind
        call collect_fixture_samples(measure_candidate, samples)
        call median_mad(samples, median, mad)
        call write_fixture_result(name, repetitions, median, mad)
    end subroutine benchmark_one

    subroutine parse_candidate(name)
        character(*), intent(in) :: name

        select case (name)
        case ("analytical")
            candidate_kind = 1
        case ("hybrid")
            candidate_kind = 2
        case ("autodiff")
            candidate_kind = 3
        case ("diagnostic")
            candidate_kind = 4
        case default
            error stop "usage: --peak-rss analytical|hybrid|autodiff|diagnostic"
        end select
    end subroutine parse_candidate

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(dp) :: elapsed_ns
        integer :: i
        real(dp) :: parameters(2), cotangent(2), x(2), jtu(2), local_sink

        cotangent(2) = -0.6_dp
        local_sink = 0.0_dp
        call timer%start()
        do i = 1, repetitions
            x(1) = 1.05_dp + 0.001_dp*real(mod(i, 101), dp)
            x(2) = 0.85_dp + 0.001_dp*real(mod(i, 79), dp)
            parameters(1) = x(1)*x(1) + x(2)
            parameters(2) = x(1) + x(2)*x(2)
            cotangent(1) = 0.02_dp*real(mod(i, 17) - 8, dp)
            select case (candidate_kind)
            case (1)
                jtu = analytical_root_vjp(x, parameters, cotangent)
            case (2)
                jtu = hybrid_root_vjp(x, parameters, cotangent)
            case (3)
                jtu = autodiff_iterations_root_vjp(parameters, cotangent)
            case default
                jtu = finite_difference_root_vjp(parameters, cotangent)
            end select
            local_sink = local_sink + sum(jtu)
        end do
        elapsed_ns = timer%elapsed_ns()/real(repetitions, dp)
        if (local_sink /= local_sink) error stop "vector-root VJP benchmark failed"
        sink = local_sink
    end function measure_candidate

end program enzyme_vector_root_vjp_hybrid
