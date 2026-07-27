module vector_root_hybrid_kernel
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr
    use fortnum_kinds, only: dp
    use fortnum_multiroot, only: multiroot_implicit_jvp
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    private

    public :: analytical_root_jvp, hybrid_root_jvp
    public :: residual_one, residual_two

    interface
        function enzyme_fwddiff(f, x1, dx1, x2, dx2, p1, dp1, p2, dp2) &
                result(df) bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), value :: x1, dx1, x2, dx2
            real(c_double), value :: p1, dp1, p2, dp2
            real(c_double) :: df
        end function enzyme_fwddiff
    end interface

contains

    pure function residual_one(x1, x2, p1, p2) result(f) &
            bind(c, name="fortnum_vector_root_residual_one")
        real(c_double), value :: x1, x2, p1, p2
        real(c_double) :: f

        f = x1*x1 + x2 - p1 + 0.0_c_double*p2
    end function residual_one

    pure function residual_two(x1, x2, p1, p2) result(f) &
            bind(c, name="fortnum_vector_root_residual_two")
        real(c_double), value :: x1, x2, p1, p2
        real(c_double) :: f

        f = x1 + x2*x2 + 0.0_c_double*p1 - p2
    end function residual_two

    function analytical_root_jvp(x, p, tp) result(dx)
        real(dp), intent(in) :: x(2), p(2), tp(2)
        real(dp) :: dx(2)
        type(fortnum_status_t) :: status

        call multiroot_implicit_jvp(analytical_residual_jvp, x, p, tp, &
            dx, status)
        if (.not. status_ok(status)) error stop "analytical vector-root JVP failed"
    end function analytical_root_jvp

    function hybrid_root_jvp(x, p, tp) result(dx)
        real(dp), intent(in) :: x(2), p(2), tp(2)
        real(dp) :: dx(2)
        type(fortnum_status_t) :: status

        call multiroot_implicit_jvp(autodiff_residual_jvp, x, p, tp, &
            dx, status)
        if (.not. status_ok(status)) error stop "hybrid vector-root JVP failed"
    end function hybrid_root_jvp

    subroutine analytical_residual_jvp(x, p, tp, jac_x, f_p_tp, context)
        real(dp), intent(in) :: x(:), p(:), tp(:)
        real(dp), intent(out) :: jac_x(size(x), size(x))
        real(dp), intent(out) :: f_p_tp(size(x))
        class(*), intent(inout), optional :: context

        if (size(p) /= 2) error stop "vector residual expects two parameters"
        jac_x(1, 1) = 2.0_dp*x(1)
        jac_x(1, 2) = 1.0_dp
        jac_x(2, 1) = 1.0_dp
        jac_x(2, 2) = 2.0_dp*x(2)
        f_p_tp = -tp
    end subroutine analytical_residual_jvp

    subroutine autodiff_residual_jvp(x, p, tp, jac_x, f_p_tp, context)
        real(dp), intent(in) :: x(:), p(:), tp(:)
        real(dp), intent(out) :: jac_x(size(x), size(x))
        real(dp), intent(out) :: f_p_tp(size(x))
        class(*), intent(inout), optional :: context

        jac_x(1, 1) = enzyme_fwddiff(c_funloc(residual_one), &
            x(1), 1.0_dp, x(2), 0.0_dp, p(1), 0.0_dp, p(2), 0.0_dp)
        jac_x(1, 2) = enzyme_fwddiff(c_funloc(residual_one), &
            x(1), 0.0_dp, x(2), 1.0_dp, p(1), 0.0_dp, p(2), 0.0_dp)
        jac_x(2, 1) = enzyme_fwddiff(c_funloc(residual_two), &
            x(1), 1.0_dp, x(2), 0.0_dp, p(1), 0.0_dp, p(2), 0.0_dp)
        jac_x(2, 2) = enzyme_fwddiff(c_funloc(residual_two), &
            x(1), 0.0_dp, x(2), 1.0_dp, p(1), 0.0_dp, p(2), 0.0_dp)
        f_p_tp(1) = enzyme_fwddiff(c_funloc(residual_one), &
            x(1), 0.0_dp, x(2), 0.0_dp, p(1), tp(1), p(2), tp(2))
        f_p_tp(2) = enzyme_fwddiff(c_funloc(residual_two), &
            x(1), 0.0_dp, x(2), 0.0_dp, p(1), tp(1), p(2), tp(2))
    end subroutine autodiff_residual_jvp

end module vector_root_hybrid_kernel

program enzyme_vector_root_hybrid
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use vector_root_hybrid_kernel, only: analytical_root_jvp, hybrid_root_jvp
    implicit none

    interface
        function peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function peak_rss_bytes
    end interface

    real(dp), parameter :: h = 1.0e-5_dp
    real(dp) :: p(2), pp(2), pm(2), tp(2), xstar(2), xp(2), xm(2)
    real(dp) :: analytical(2), got(2), reference(2)
    real(dp) :: analytical_error, hybrid_error
    character(32) :: argument, candidate

    call get_command_argument(1, argument)
    call get_command_argument(2, candidate)
    if (trim(argument) == "--benchmark") then
        call run_benchmark()
    else if (trim(argument) == "--peak-rss") then
        call run_peak_rss(trim(candidate))
    else
        p = [2.0_dp, 2.0_dp]
        tp = [0.4_dp, -0.6_dp]
        xstar = solve_root(p)
        pp = p + h*tp
        pm = p - h*tp
        xp = solve_root(pp)
        xm = solve_root(pm)
        reference = (xp - xm)/(2.0_dp*h)
        analytical = analytical_root_jvp(xstar, p, tp)
        got = hybrid_root_jvp(xstar, p, tp)
        analytical_error = maxval(abs(analytical - reference))
        hybrid_error = maxval(abs(got - reference))

        if (analytical_error > 1.0e-9_dp) then
            print *, "analytical vector-root JVP mismatch", analytical, reference
            error stop 1
        end if
        if (hybrid_error > 1.0e-9_dp) then
            print *, "hybrid vector-root JVP mismatch", got, reference
            error stop 1
        end if
        print *, "PASS vector-root JVPs", analytical_error, hybrid_error
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
        integer(int64), parameter :: reps = 2000000_int64
        real(dp) :: analytical_times(samples), hybrid_times(samples), sink
        integer :: sample

        do sample = 1, 3
            call time_candidate("analytical", reps/100_int64, sink)
            call time_candidate("hybrid", reps/100_int64, sink)
        end do
        do sample = 1, samples
            call time_candidate("analytical", reps, analytical_times(sample))
            call time_candidate("hybrid", reps, hybrid_times(sample))
        end do
        call report("analytical", analytical_times, reps)
        call report("hybrid", hybrid_times, reps)
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
        real(dp) :: parameters(2), direction(2), x(2), dx(2), sink

        direction(2) = -0.6_dp
        sink = 0.0_dp
        call system_clock(start, rate)
        do i = 1_int64, reps
            x(1) = 1.05_dp + 0.001_dp*real(mod(i, 101_int64), dp)
            x(2) = 0.85_dp + 0.001_dp*real(mod(i, 79_int64), dp)
            parameters(1) = x(1)*x(1) + x(2)
            parameters(2) = x(1) + x(2)*x(2)
            direction(1) = 0.02_dp*real(mod(i, 17_int64) - 8_int64, dp)
            if (name == "analytical") then
                dx = analytical_root_jvp(x, parameters, direction)
            else
                dx = hybrid_root_jvp(x, parameters, direction)
            end if
            sink = sink + sum(dx)
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

end program enzyme_vector_root_hybrid
