program bench_static_registry
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_build_selection, only: fortnum_lookup_build_selection
    implicit none

    integer, parameter :: samples = 15
    integer, parameter :: query_count = 4
    integer(int64), parameter :: reps = 500000_int64
    character(48), parameter :: operators(query_count) = [character(48) :: &
        "dawson_outer", "linear_solve", "linear_solve", "multiroot_implicit"]
    character(16), parameter :: products(query_count) = [character(16) :: &
        "jvp", "jvp", "vjp", "jvp"]
    character(128), parameter :: workloads(query_count) = [character(128) :: &
        "scalar calls, x cycled over [0.65, 0.75], one direction", &
        "16x16 dense system, 200000 directions per sample, primal LU reusable", &
        "16x16 dense system, 200000 cotangents per sample, transposed primal LU reusable", &
        "16x16 dense residual Jacobian, 200000 parameter directions per sample"]
    character(64), parameter :: candidates(query_count) = [character(64) :: &
        "analytical", "reuse_primal_lu", "reuse_transposed_primal_lu", &
        "default_dense_solve"]
    real(dp) :: warmup
    integer :: sample
    character(16) :: candidate

    call get_command_argument(1, candidate)
    if ((trim(candidate) /= "registry") .and. (trim(candidate) /= "direct")) then
        error stop "usage: bench_static_registry registry|direct"
    end if

    do sample = 1, 3
        warmup = run_sample(candidate, reps/10_int64)
    end do
    do sample = 1, samples
        write (*, "(f0.4)") run_sample(candidate, reps)
    end do

contains

    function run_sample(name, count) result(ns_per_call)
        character(*), intent(in) :: name
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call
        integer(int64) :: k, tick0, tick1, rate
        integer :: index
        character(64) :: selected
        logical :: found
        real(dp) :: sink

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do k = 1_int64, count
            index = int(mod(k - 1_int64, int(query_count, int64))) + 1
            if (trim(name) == "registry") then
                call fortnum_lookup_build_selection(operators(index), &
                    products(index), workloads(index), selected, found)
            else
                selected = candidates(index)
                found = .true.
            end if
            if (found) sink = sink + real(len_trim(selected), dp)
        end do
        call system_clock(tick1)
        if (sink <= 0.0_dp) error stop "benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp &
            / (real(rate, dp)*real(count, dp))
    end function run_sample

end program bench_static_registry
