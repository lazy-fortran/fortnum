module scalar_root_hybrid_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_generated_enzyme_scalar_root_newton, only: &
        fortnum_enzyme_scalar_root_newton_jvp
    use fortnum_generated_enzyme_scalar_root_residual, only: &
        fortnum_enzyme_scalar_root_residual_jvp
    use fortnum_generated_scalar_root_residual_jvp, only: &
        fortnum_scalar_root_residual_jvp_kernel
    use fortnum_kinds, only: dp
    use fortnum_roots, only: root_implicit_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    public :: analytical_root_jvp, hybrid_root_jvp, &
        autodiff_iterations_root_jvp, finite_difference_root_jvp, residual

contains

    pure function residual(x, p1, p2) result(f) &
            bind(c, name="fortnum_scalar_root_residual")
        real(c_double), value :: x, p1, p2
        real(c_double) :: f

        f = x*x*x + p1*x - p2
    end function residual

    pure function newton_solve(x0, p1, p2) result(x) &
            bind(c, name="fortnum_scalar_root_newton_solve")
        real(c_double), value :: x0, p1, p2
        real(c_double) :: x
        integer :: iteration

        x = 1.0_c_double + 0.0_c_double*x0
        do iteration = 1, 12
            x = x - residual(x, p1, p2)/(3.0_c_double*x*x + p1)
        end do
    end function newton_solve

    function analytical_root_jvp(x, p, tp) result(dx)
        real(dp), intent(in) :: x, p(2), tp(2)
        real(dp) :: dx
        type(fortnum_status_t) :: status

        call root_implicit_jvp(analytical_residual_jvp, x, p, tp, dx, status)
        if (.not. status_ok(status)) error stop "analytical root JVP failed"
    end function analytical_root_jvp

    function hybrid_root_jvp(x, p, tp) result(dx)
        real(dp), intent(in) :: x, p(2), tp(2)
        real(dp) :: dx
        type(fortnum_status_t) :: status

        call root_implicit_jvp(autodiff_residual_jvp, x, p, tp, dx, status)
        if (.not. status_ok(status)) error stop "hybrid root JVP failed"
    end function hybrid_root_jvp

    function autodiff_iterations_root_jvp(p, tp) result(dx)
        real(dp), intent(in) :: p(2), tp(2)
        real(dp) :: dx

        dx = fortnum_enzyme_scalar_root_newton_jvp(0.0_dp, 0.0_dp, &
            p(1), tp(1), p(2), tp(2))
    end function autodiff_iterations_root_jvp

    function finite_difference_root_jvp(p, tp) result(dx)
        real(dp), intent(in) :: p(2), tp(2)
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: dx

        dx = (newton_solve(0.0_dp, p(1) + h*tp(1), p(2) + h*tp(2)) &
            - newton_solve(0.0_dp, p(1) - h*tp(1), p(2) - h*tp(2))) &
            /(2.0_dp*h)
    end function finite_difference_root_jvp

    subroutine analytical_residual_jvp(x, p, tp, f_x, f_p_tp, context)
        real(dp), intent(in) :: x, p(:), tp(:)
        real(dp), intent(out) :: f_x, f_p_tp
        class(*), intent(inout), optional :: context

        call fortnum_scalar_root_residual_jvp_kernel( &
            x, p(1), tp(1), tp(2), f_x, f_p_tp)
    end subroutine analytical_residual_jvp

    subroutine autodiff_residual_jvp(x, p, tp, f_x, f_p_tp, context)
        real(dp), intent(in) :: x, p(:), tp(:)
        real(dp), intent(out) :: f_x, f_p_tp
        class(*), intent(inout), optional :: context

        f_x = fortnum_enzyme_scalar_root_residual_jvp(x, 1.0_dp, &
            p(1), 0.0_dp, p(2), 0.0_dp)
        f_p_tp = fortnum_enzyme_scalar_root_residual_jvp(x, 0.0_dp, &
            p(1), tp(1), p(2), tp(2))
    end subroutine autodiff_residual_jvp

end module scalar_root_hybrid_kernel

program enzyme_scalar_root_hybrid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, write_fixture_result
    use scalar_root_hybrid_kernel, only: analytical_root_jvp, &
        hybrid_root_jvp, autodiff_iterations_root_jvp, &
        finite_difference_root_jvp, residual
    implicit none

    integer, parameter :: repetitions = 100000
    real(dp), parameter :: h = 1.0e-5_dp
    real(dp) :: p(2), pp(2), pm(2), tp(2), xstar, got, reference
    real(dp) :: analytical, autodiff_iterations, diagnostic
    real(dp) :: errors(4)
    real(dp) :: samples(fixture_sample_count), sink
    character(32) :: argument, candidate
    integer :: candidate_kind

    call get_command_argument(1, argument)
    call get_command_argument(2, candidate)
    if (trim(argument) == "--benchmark") then
        call run_benchmark()
    else if (trim(argument) == "--peak-rss") then
        call run_peak_rss(trim(candidate))
    else
        p = [0.7_dp, 3.0_dp]
        tp = [0.4_dp, -0.6_dp]
        xstar = solve_root(p)
        pp = p + h*tp
        pm = p - h*tp
        reference = (solve_root(pp) - solve_root(pm))/(2.0_dp*h)
        got = hybrid_root_jvp(xstar, p, tp)
        analytical = analytical_root_jvp(xstar, p, tp)
        autodiff_iterations = autodiff_iterations_root_jvp(p, tp)
        diagnostic = finite_difference_root_jvp(p, tp)
        errors(1) = abs(analytical - reference)
        errors(2) = abs(got - reference)
        errors(3) = abs(autodiff_iterations - reference)
        errors(4) = abs(diagnostic - reference)

        if (maxval(errors) > 1.0e-9_dp) then
            print *, "scalar-root JVP mismatch", analytical, got, &
                autodiff_iterations, diagnostic, reference
            error stop 1
        end if
        print *, "PASS scalar-root JVP tournament errors", errors
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
        real(dp) :: parameters(2), direction(2), x, local_sink

        parameters(1) = 0.7_dp
        direction = [0.4_dp, -0.6_dp]
        local_sink = 0.0_dp
        call timer%start()
        do i = 1, repetitions
            x = 1.15_dp + 0.001_dp*real(mod(i, 101), dp)
            parameters(2) = x*x*x + parameters(1)*x
            select case (candidate_kind)
            case (1)
                local_sink = local_sink + analytical_root_jvp( &
                    x, parameters, direction)
            case (2)
                local_sink = local_sink + hybrid_root_jvp( &
                    x, parameters, direction)
            case (3)
                local_sink = local_sink + autodiff_iterations_root_jvp( &
                    parameters, direction)
            case default
                local_sink = local_sink + finite_difference_root_jvp( &
                    parameters, direction)
            end select
        end do
        elapsed_ns = timer%elapsed_ns()/real(repetitions, dp)
        if (local_sink /= local_sink) error stop "scalar-root JVP benchmark failed"
        sink = local_sink
    end function measure_candidate

end program enzyme_scalar_root_hybrid
