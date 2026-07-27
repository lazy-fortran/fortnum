module dawson_outer_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    private
    public :: outer, outer_autodiff

    interface
        function fortnum_dawson_kernel(x) result(f) &
                bind(c, name="fortnum_dawson_kernel")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: f
        end function fortnum_dawson_kernel

        function fortnum_dawson_kernel_autodiff(x) result(f) &
                bind(c, name="fortnum_dawson_kernel_autodiff")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: f
        end function fortnum_dawson_kernel_autodiff
    end interface

contains

    function outer(x) result(y) bind(c, name="fortnum_dawson_outer")
        real(c_double), value :: x
        real(c_double) :: y, f
        f = fortnum_dawson_kernel(x)
        y = sin(f) + f*f
    end function outer

    function outer_autodiff(x) result(y) bind(c, name="fortnum_dawson_outer_autodiff")
        real(c_double), value :: x
        real(c_double) :: y, f
        f = fortnum_dawson_kernel_autodiff(x)
        y = sin(f) + f*f
    end function outer_autodiff

end module dawson_outer_kernel

program enzyme_dawson_hybrid
    use, intrinsic :: iso_c_binding, only: c_double, c_funloc, c_funptr, c_int
    use, intrinsic :: iso_fortran_env, only: int64
    use dawson_outer_kernel, only: outer, outer_autodiff
    use fortnum_generated_dawson_outer, only: fortnum_dawson_outer_kernel
    implicit none

    interface
        function enzyme_fwddiff(f, x, dx) result(dy) &
                bind(c, name="__enzyme_fwddiff")
            import :: c_double, c_funptr
            type(c_funptr), value :: f
            real(c_double), value :: x, dx
            real(c_double) :: dy
        end function enzyme_fwddiff

        subroutine rule_reset() bind(c, name="fortnum_dawson_rule_reset")
        end subroutine rule_reset

        function rule_calls() result(n) bind(c, name="fortnum_dawson_rule_calls")
            import :: c_int
            integer(c_int) :: n
        end function rule_calls

        subroutine rule_disable_count() &
                bind(c, name="fortnum_dawson_rule_disable_count")
        end subroutine rule_disable_count

        function fortnum_dawson_kernel(x) result(f) &
                bind(c, name="fortnum_dawson_kernel")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: f
        end function fortnum_dawson_kernel
    end interface

    real(c_double), parameter :: x = 0.7_c_double
    real(c_double), parameter :: dx = -0.4_c_double
    real(c_double), parameter :: h = 1.0e-6_c_double
    real(c_double) :: got, reference, scale
    character(32) :: argument

    call get_command_argument(1, argument)
    if (trim(argument) == "--benchmark") then
        call run_benchmark()
        stop
    end if
    call rule_reset()
    got = enzyme_fwddiff(c_funloc(outer), x, dx)
    reference = (outer(x + h*dx) - outer(x - h*dx))/(2.0_c_double*h)
    scale = max(1.0_c_double, abs(reference))

    if (abs(got - reference)/scale > 2.0e-9_c_double) then
        print *, "hybrid Dawson JVP mismatch", got, reference
        error stop 1
    end if
    if (rule_calls() /= 1_c_int) then
        print *, "analytical Dawson rule was not selected", rule_calls()
        error stop 1
    end if
    print *, "PASS hybrid Dawson JVP", got

contains

    function analytical_jvp(x, dx) result(dy)
        real(c_double), intent(in) :: x, dx
        real(c_double) :: dy, f, value
        f = fortnum_dawson_kernel(x)
        call fortnum_dawson_outer_kernel(x, f, dx, value, dy)
    end function analytical_jvp

    subroutine run_benchmark()
        integer, parameter :: samples = 15
        integer, parameter :: reps = 200000
        real(c_double) :: analytical(samples), autodiff(samples), hybrid(samples)
        real(c_double) :: sink
        integer :: i

        call rule_disable_count()
        do i = 1, 3
            call time_candidate(i, 2000, sink)
        end do
        do i = 1, samples
            call time_candidate(1, reps, analytical(i))
            call time_candidate(2, reps, autodiff(i))
            call time_candidate(3, reps, hybrid(i))
        end do
        call report("analytical", analytical, reps)
        call report("autodiff", autodiff, reps)
        call report("hybrid", hybrid, reps)
        if (sink == huge(sink)) error stop 1
    end subroutine run_benchmark

    subroutine time_candidate(which, reps, elapsed_ns)
        integer, intent(in) :: which, reps
        real(c_double), intent(out) :: elapsed_ns
        integer(int64) :: start, finish, rate
        integer :: i
        real(c_double) :: xi, sink

        sink = 0.0_c_double
        call system_clock(start, rate)
        do i = 1, reps
            xi = 0.65_c_double + 0.001_c_double*real(mod(i, 101), c_double)
            select case (which)
            case (1)
                sink = sink + analytical_jvp(xi, dx)
            case (2)
                sink = sink + enzyme_fwddiff(c_funloc(outer_autodiff), xi, dx)
            case (3)
                sink = sink + enzyme_fwddiff(c_funloc(outer), xi, dx)
            end select
        end do
        call system_clock(finish)
        elapsed_ns = 1.0e9_c_double*real(finish - start, c_double)/ &
            (real(rate, c_double)*real(reps, c_double))
        if (sink == huge(sink)) print *, sink
    end subroutine time_candidate

    subroutine report(name, values, reps)
        character(*), intent(in) :: name
        real(c_double), intent(in) :: values(:)
        integer, intent(in) :: reps
        real(c_double) :: ordered(size(values)), deviations(size(values))
        real(c_double) :: median, mad

        ordered = values
        call sort_values(ordered)
        median = ordered((size(ordered) + 1)/2)
        deviations = abs(values - median)
        call sort_values(deviations)
        mad = deviations((size(deviations) + 1)/2)
        write (*, "(a,',',i0,',',f12.4,',',f12.4)") name, reps, median, mad
    end subroutine report

    subroutine sort_values(values)
        real(c_double), intent(inout) :: values(:)
        real(c_double) :: tmp
        integer :: i, j
        do i = 1, size(values) - 1
            do j = i + 1, size(values)
                if (values(j) < values(i)) then
                    tmp = values(i)
                    values(i) = values(j)
                    values(j) = tmp
                end if
            end do
        end do
    end subroutine sort_values

end program enzyme_dawson_hybrid
