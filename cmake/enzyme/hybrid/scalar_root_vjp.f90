module scalar_root_vjp_hybrid_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use fortnum_kinds, only: dp
    use fortnum_roots, only: root_implicit_vjp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    type, bind(c) :: residual_gradient_t
        real(c_double) :: values(3)
    end type residual_gradient_t

    public :: analytical_root_vjp, hybrid_root_vjp, &
        autodiff_iterations_root_vjp, finite_difference_root_vjp, residual

    interface
        function enzyme_autodiff(f, x, p1, p2) result(gradient) &
                bind(c, name="__enzyme_autodiff")
            import :: c_double, c_funptr, residual_gradient_t
            type(c_funptr), value :: f
            real(c_double), value :: x, p1, p2
            type(residual_gradient_t) :: gradient
        end function enzyme_autodiff
    end interface

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
        type(residual_gradient_t) :: gradient

        gradient = enzyme_autodiff(c_funloc(newton_solve), &
            0.0_dp, p(1), p(2))
        jtu = u*gradient%values(2:3)
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

        f_x = 3.0_dp*x*x + p(1)
        f_p_t_u(1) = u*x
        f_p_t_u(2) = -u
    end subroutine analytical_residual_vjp

    subroutine autodiff_residual_vjp(x, p, u, f_x, f_p_t_u, context)
        real(dp), intent(in) :: x, p(:), u
        real(dp), intent(out) :: f_x, f_p_t_u(size(p))
        class(*), intent(inout), optional :: context
        type(residual_gradient_t) :: gradient

        gradient = enzyme_autodiff(c_funloc(residual), x, p(1), p(2))
        f_x = gradient%values(1)
        f_p_t_u = u*gradient%values(2:3)
    end subroutine autodiff_residual_vjp

end module scalar_root_vjp_hybrid_kernel

program enzyme_scalar_root_vjp_hybrid
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use scalar_root_vjp_hybrid_kernel, only: analytical_root_vjp, &
        hybrid_root_vjp, autodiff_iterations_root_vjp, &
        finite_difference_root_vjp, residual
    implicit none

    interface
        function peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function peak_rss_bytes
    end interface

    real(dp), parameter :: h = 1.0e-5_dp
    real(dp), parameter :: u = 1.3_dp
    real(dp) :: p(2), pp(2), pm(2), xstar, got(2), reference(2)
    real(dp) :: analytical(2), analytical_error, hybrid_error
    real(dp) :: autodiff_iterations(2), autodiff_iterations_error
    real(dp) :: diagnostic(2), diagnostic_error
    character(32) :: argument, candidate
    integer :: parameter

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
        integer, parameter :: samples = 15
        integer(int64), parameter :: reps = 100000_int64
        real(dp) :: analytical(samples), hybrid(samples)
        real(dp) :: autodiff_iterations(samples), diagnostic(samples), sink
        integer :: sample

        do sample = 1, 3
            call time_candidate("analytical", reps/100_int64, sink)
            call time_candidate("hybrid", reps/100_int64, sink)
            call time_candidate("autodiff", reps/100_int64, sink)
            call time_candidate("diagnostic", reps/100_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate("analytical", reps, analytical(sample))
            call time_candidate("hybrid", reps, hybrid(sample))
            call time_candidate("autodiff", reps, autodiff_iterations(sample))
            call time_candidate("diagnostic", reps, diagnostic(sample))
        end do
        call report("analytical", analytical, reps)
        call report("hybrid", hybrid, reps)
        call report("autodiff", autodiff_iterations, reps)
        call report("diagnostic", diagnostic, reps)
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
        real(dp) :: parameters(2), cotangent, x, jtu(2), sink

        parameters(1) = 0.7_dp
        sink = 0.0_dp
        call system_clock(start, rate)
        do i = 1_int64, reps
            x = 1.15_dp + 0.001_dp*real(mod(i, 101_int64), dp)
            parameters(2) = x*x*x + parameters(1)*x
            cotangent = 0.8_dp + 0.001_dp*real(mod(i, 17_int64), dp)
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

end program enzyme_scalar_root_vjp_hybrid
