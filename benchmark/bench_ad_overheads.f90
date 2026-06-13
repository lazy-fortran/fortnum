program bench_ad_overheads
    ! Measures derivative-product overhead ratios JVP/primal, VJP/primal,
    ! grad/primal, HVP/primal for fortnum kernels. No real products exist yet
    ! (#40); a placeholder kernel and its hand-coded JVP seed the table so the
    ! harness builds and runs, and the output format is ready for #40 to add
    ! real (primal, derivative) pairs via registry%add(..., deriv=...).
    !
    ! Reuses fortnum_benchmark so the ns/call timing, the JSON schema, and
    ! gate.py stay shared with bench_main. The seed pair bounds the harness
    ! overhead itself: the JVP/primal ratio of two arithmetically similar
    ! kernels should sit near one.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark, only: bench_registry_t
    implicit none

    type(bench_registry_t) :: registry
    integer(int64)         :: reps
    logical                :: as_json
    integer                :: iarg
    character(64)          :: arg

    reps = 20000000_int64

    as_json = .false.
    do iarg = 1, command_argument_count()
        call get_command_argument(iarg, arg)
        if (trim(arg) == "--json") as_json = .true.
    end do

    ! Seed: placeholder primal and its forward product. The deriv= kernel
    ! makes registry%run report the JVP/primal ratio column. #40 replaces this
    ! with real (foo, foo_jvp) pairs and adds vjp/grad/hvp pairs alongside.
    call registry%add("placeholder_quadratic", primal_quadratic, &
        deriv=jvp_quadratic)

    call registry%run(reps, as_json=as_json)

contains

    ! Placeholder primal: a scalar reduction of a fixed quadratic so the timed
    ! work is real and the optimiser cannot fold it away (the result is
    ! consumed by the timing accumulator).
    function primal_quadratic() result(sink)
        real(dp) :: sink
        real(dp), parameter :: x(4) = &
            [0.5_dp, 1.5_dp, -3.0_dp, 2.0_dp]
        sink = sum(x*x)
    end function primal_quadratic

    ! Forward product of the same kernel at a fixed direction; stands in for a
    ! real foo_jvp so the ratio column has a value.
    function jvp_quadratic() result(sink)
        real(dp) :: sink
        real(dp), parameter :: x(4) = &
            [0.5_dp, 1.5_dp, -3.0_dp, 2.0_dp]
        real(dp), parameter :: v(4) = &
            [1.0_dp, -1.0_dp, 0.25_dp, 0.5_dp]
        sink = sum(2.0_dp*x*v)
    end function jvp_quadratic

end program bench_ad_overheads
