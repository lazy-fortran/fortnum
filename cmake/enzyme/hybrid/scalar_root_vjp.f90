module scalar_root_vjp_hybrid_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_generated_enzyme_scalar_root_vjp_newton, only: &
        fortnum_enzyme_scalar_root_vjp_newton_vjp
    use fortnum_generated_enzyme_scalar_root_vjp_residual, only: &
        fortnum_enzyme_scalar_root_vjp_residual_vjp
    use fortnum_generated_scalar_root_residual_vjp, only: &
        fortnum_scalar_root_residual_vjp_kernel
    use fortnum_kinds, only: dp
    use fortnum_roots, only: root_implicit_vjp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    public :: analytical_root_vjp, hybrid_root_vjp, &
        autodiff_iterations_root_vjp, finite_difference_root_vjp, residual

contains

    pure function residual(x, p1, p2) result(f) &
            bind(c, name="fortnum_scalar_root_vjp_residual")
        real(c_double), value :: x, p1, p2
        real(c_double) :: f

        f = x*x*x + p1*x - p2
    end function residual

    pure function newton_solve(x0, p1, p2) result(x) &
            bind(c, name="fortnum_scalar_root_vjp_newton_solve")
        real(c_double), value :: x0, p1, p2
        real(c_double) :: x
        integer :: iteration

        x = 1.0_c_double + 0.0_c_double*x0
        do iteration = 1, 12
            x = x - residual(x, p1, p2)/(3.0_c_double*x*x + p1)
        end do
    end function newton_solve

    function analytical_root_vjp(x, p, u) result(jtu)
        real(dp), intent(in) :: x, p(2), u
        real(dp) :: jtu(2)
        type(fortnum_status_t) :: status

        call root_implicit_vjp(analytical_residual_vjp, x, p, u, jtu, status)
        if (.not. status_ok(status)) error stop "analytical root VJP failed"
    end function analytical_root_vjp

    function hybrid_root_vjp(x, p, u) result(jtu)
        real(dp), intent(in) :: x, p(2), u
        real(dp) :: jtu(2)
        type(fortnum_status_t) :: status

        call root_implicit_vjp(autodiff_residual_vjp, x, p, u, jtu, status)
        if (.not. status_ok(status)) error stop "hybrid root VJP failed"
    end function hybrid_root_vjp

    function autodiff_iterations_root_vjp(p, u) result(jtu)
        real(dp), intent(in) :: p(2), u
        real(dp) :: jtu(2)
        real(dp) :: ignored_xbar

        call fortnum_enzyme_scalar_root_vjp_newton_vjp( &
            0.0_dp, p(1), p(2), u, ignored_xbar, jtu(1), jtu(2))
    end function autodiff_iterations_root_vjp

    function finite_difference_root_vjp(p, u) result(jtu)
        real(dp), intent(in) :: p(2), u
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: jtu(2), pm(2), pp(2)
        integer :: parameter

        do parameter = 1, 2
            pp = p
            pm = p
            pp(parameter) = pp(parameter) + h
            pm(parameter) = pm(parameter) - h
            jtu(parameter) = u*( &
                newton_solve(0.0_dp, pp(1), pp(2)) &
                - newton_solve(0.0_dp, pm(1), pm(2)))/(2.0_dp*h)
        end do
    end function finite_difference_root_vjp

    subroutine analytical_residual_vjp(x, p, u, f_x, f_p_t_u, context)
        real(dp), intent(in) :: x, p(:), u
        real(dp), intent(out) :: f_x, f_p_t_u(size(p))
        class(*), intent(inout), optional :: context

        call fortnum_scalar_root_residual_vjp_kernel( &
            x, p(1), u, f_x, f_p_t_u(1), f_p_t_u(2))
    end subroutine analytical_residual_vjp

    subroutine autodiff_residual_vjp(x, p, u, f_x, f_p_t_u, context)
        real(dp), intent(in) :: x, p(:), u
        real(dp), intent(out) :: f_x, f_p_t_u(size(p))
        class(*), intent(inout), optional :: context
        real(dp) :: p1bar, p2bar

        call fortnum_enzyme_scalar_root_vjp_residual_vjp( &
            x, p(1), p(2), 1.0_dp, f_x, p1bar, p2bar)
        f_p_t_u = u*[p1bar, p2bar]
    end subroutine autodiff_residual_vjp

end module scalar_root_vjp_hybrid_kernel

program enzyme_scalar_root_vjp_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, write_fixture_result
    use scalar_root_vjp_hybrid_kernel, only: analytical_root_vjp, &
        hybrid_root_vjp, autodiff_iterations_root_vjp, &
        finite_difference_root_vjp, residual
    implicit none

    integer, parameter :: repetitions = 100000
    real(dp), parameter :: h = 1.0e-5_dp
    real(dp), parameter :: u = 1.3_dp
    real(dp) :: p(2), pp(2), pm(2), xstar, got(2), reference(2)
    real(dp) :: analytical(2), analytical_error, hybrid_error
    real(dp) :: autodiff_iterations(2), autodiff_iterations_error
    real(dp) :: diagnostic(2), diagnostic_error
    real(dp) :: samples(fixture_sample_count), sink
    character(32) :: argument, candidate
    integer :: candidate_kind, parameter

    call get_command_argument(1, argument)
    call get_command_argument(2, candidate)
    if (trim(argument) == "--benchmark") then
        call run_benchmark()
    else if (trim(argument) == "--peak-rss") then
        call run_peak_rss(trim(candidate))
    else
        p = [0.7_dp, 3.0_dp]
        xstar = solve_root(p)
        do parameter = 1, size(p)
            pp = p
            pm = p
            pp(parameter) = pp(parameter) + h
            pm(parameter) = pm(parameter) - h
            reference(parameter) = u*(solve_root(pp) - solve_root(pm)) &
                /(2.0_dp*h)
        end do
        got = hybrid_root_vjp(xstar, p, u)
        analytical = analytical_root_vjp(xstar, p, u)
        autodiff_iterations = autodiff_iterations_root_vjp(p, u)
        diagnostic = finite_difference_root_vjp(p, u)
        analytical_error = maxval(abs(analytical - reference))
        hybrid_error = maxval(abs(got - reference))
        autodiff_iterations_error = maxval(abs(autodiff_iterations - reference))
        diagnostic_error = maxval(abs(diagnostic - reference))

        if (max(analytical_error, hybrid_error, autodiff_iterations_error, &
            diagnostic_error) > 1.0e-9_dp) then
            print *, "scalar-root VJP mismatch", analytical_error, &
                hybrid_error, autodiff_iterations_error, diagnostic_error
            error stop 1
        end if
        print *, "PASS scalar-root VJP tournament", analytical_error, &
            hybrid_error, autodiff_iterations_error, diagnostic_error
    end if

contains

    pure function solve_root(parameters) result(x)
        real(dp), intent(in) :: parameters(2)
        real(dp) :: x, a, b, middle, fm
        integer :: iteration

        a = 0.0_dp
        b = 3.0_dp
        do iteration = 1, 100
            middle = 0.5_dp*(a + b)
            fm = residual(middle, parameters(1), parameters(2))
            if (fm > 0.0_dp) then
                b = middle
            else
                a = middle
            end if
        end do
        x = 0.5_dp*(a + b)
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
        real(dp) :: parameters(2), cotangent, x, jtu(2), local_sink

        parameters(1) = 0.7_dp
        local_sink = 0.0_dp
        call timer%start()
        do i = 1, repetitions
            x = 1.15_dp + 0.001_dp*real(mod(i, 101), dp)
            parameters(2) = x*x*x + parameters(1)*x
            cotangent = 0.8_dp + 0.001_dp*real(mod(i, 17), dp)
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
        if (local_sink /= local_sink) error stop "scalar-root VJP benchmark failed"
        sink = local_sink
    end function measure_candidate

end program enzyme_scalar_root_vjp_hybrid
