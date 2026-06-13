program bench_main
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark, only: bench_registry_t
    implicit none

    type(bench_registry_t) :: registry
    integer(int64)         :: reps

    ! Repetition count sized so even a sub-nanosecond kernel runs long enough
    ! for system_clock resolution to give a stable ns/call figure.
    reps = 50000000_int64

    ! Seed cases.  Both run today; they bound the harness overhead itself.
    call registry%add("version_string_len", kernel_version_len)
    call registry%add("arith_kernel_fma", kernel_arith)

    ! M6 (#36) wires derivative kernels here by passing the optional deriv=
    ! argument to registry%add, e.g.
    !   call registry%add("rosenbrock", primal_rosenbrock, deriv=vjp_rosenbrock)
    ! The run() table then fills the "deriv ns" and "ratio" columns.

    call registry%run(reps)

contains

    ! Touches a fortnum public symbol so the seed benchmark exercises the
    ! library, not just intrinsics.
    function kernel_version_len() result(sink)
        use fortnum_version, only: fortnum_version_string
        real(dp) :: sink
        sink = real(len_trim(fortnum_version_string), dp)
    end function kernel_version_len

    ! Trivial fused multiply-add kernel as a fixed arithmetic baseline.  Uses
    ! a static counter so successive calls differ and cannot be folded.
    function kernel_arith() result(sink)
        real(dp)       :: sink
        integer(int64), save :: n = 0_int64
        real(dp)       :: x
        n = n + 1_int64
        x = real(mod(n, 1024_int64), dp)
        sink = x * x + 1.0_dp
    end function kernel_arith

end program bench_main
