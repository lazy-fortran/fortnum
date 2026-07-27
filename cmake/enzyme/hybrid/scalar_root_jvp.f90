module scalar_root_hybrid_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use fortnum_kinds, only: dp
    use fortnum_roots, only: root_implicit_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    public :: analytical_root_jvp, hybrid_root_jvp, residual

    interface
        function enzyme_fwddiff(f, x, dx, p1, dp1, p2, dp2) result(df) &
                bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), value :: x, dx, p1, dp1, p2, dp2
            real(c_double) :: df
        end function enzyme_fwddiff
    end interface

contains

    pure function residual(x, p1, p2) result(f) &
            bind(c, name="fortnum_scalar_root_residual")
        real(c_double), value :: x, p1, p2
        real(c_double) :: f

        f = x*x*x + p1*x - p2
    end function residual

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

    subroutine analytical_residual_jvp(x, p, tp, f_x, f_p_tp, context)
        real(dp), intent(in) :: x, p(:), tp(:)
        real(dp), intent(out) :: f_x, f_p_tp
        class(*), intent(inout), optional :: context

        f_x = 3.0_dp*x*x + p(1)
        f_p_tp = x*tp(1) - tp(2)
    end subroutine analytical_residual_jvp

    subroutine autodiff_residual_jvp(x, p, tp, f_x, f_p_tp, context)
        real(dp), intent(in) :: x, p(:), tp(:)
        real(dp), intent(out) :: f_x, f_p_tp
        class(*), intent(inout), optional :: context

        f_x = enzyme_fwddiff(c_funloc(residual), x, 1.0_dp, &
            p(1), 0.0_dp, p(2), 0.0_dp)
        f_p_tp = enzyme_fwddiff(c_funloc(residual), x, 0.0_dp, &
            p(1), tp(1), p(2), tp(2))
    end subroutine autodiff_residual_jvp

end module scalar_root_hybrid_kernel

program enzyme_scalar_root_hybrid
    use, intrinsic :: iso_c_binding, only: c_double, c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use scalar_root_hybrid_kernel, only: analytical_root_jvp, &
        hybrid_root_jvp, residual
    implicit none

    interface
        function peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function peak_rss_bytes
    end interface

    real(dp), parameter :: h = 1.0e-5_dp
    real(dp) :: p(2), pp(2), pm(2), tp(2), xstar, got, reference
    character(32) :: argument, candidate

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

        if (abs(got - reference) > 1.0e-9_dp) then
            print *, "hybrid scalar-root JVP mismatch", got, reference
            error stop 1
        end if
        print *, "PASS hybrid scalar-root JVP", got
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
        integer(int64), parameter :: reps = 2000000_int64
        real(dp) :: analytical(samples), hybrid(samples), sink
        integer :: sample

        do sample = 1, 3
            call time_candidate("analytical", reps/100_int64, sink)
            call time_candidate("hybrid", reps/100_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate("analytical", reps, analytical(sample))
            call time_candidate("hybrid", reps, hybrid(sample))
        end do
        call report("analytical", analytical, reps)
        call report("hybrid", hybrid, reps)
    end subroutine run_benchmark

    subroutine run_peak_rss(name)
        character(*), intent(in) :: name
        integer(int64), parameter :: reps = 2000000_int64
        real(dp) :: elapsed

        if ((name /= "analytical") .and. (name /= "hybrid")) then
            error stop "usage: --peak-rss analytical|hybrid"
        end if
        call time_candidate(name, reps, elapsed)
        write (*, "(i0)") peak_rss_bytes()
    end subroutine run_peak_rss

    subroutine time_candidate(name, reps, elapsed_ns)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: reps
        real(dp), intent(out) :: elapsed_ns
        integer(int64) :: i, start, finish, rate
        real(dp) :: parameters(2), direction(2), x, sink

        parameters(1) = 0.7_dp
        direction = [0.4_dp, -0.6_dp]
        sink = 0.0_dp
        call system_clock(start, rate)
        do i = 1_int64, reps
            x = 1.15_dp + 0.001_dp*real(mod(i, 101_int64), dp)
            parameters(2) = x*x*x + parameters(1)*x
            if (name == "analytical") then
                sink = sink + analytical_root_jvp(x, parameters, direction)
            else
                sink = sink + hybrid_root_jvp(x, parameters, direction)
            end if
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

end program enzyme_scalar_root_hybrid
