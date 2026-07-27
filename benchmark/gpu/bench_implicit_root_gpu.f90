program bench_implicit_root_gpu
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_gpu_implicit_root_wrapper, only: implicit_root_jvp_batch
    implicit none

    integer, parameter :: samples = 31
    real(dp), allocatable :: parameters(:), directions(:), roots(:), tangents(:)
    logical, allocatable :: reliable(:)
    character(16) :: argument, residency, option
    integer :: i, n, sample
    logical :: peak_only
    real(dp) :: sink

    call get_command_argument(1, argument)
    read (argument, *) n
    call get_command_argument(2, residency)
    call get_command_argument(3, option)
    if (n < 1) error stop "batch size must be positive"
    if (trim(residency) /= "transfer" .and. &
        trim(residency) /= "resident") then
        error stop "residency must be transfer or resident"
    end if
    peak_only = trim(option) == "--peak-rss"

    allocate (parameters(n), directions(n), roots(n), tangents(n), reliable(n))
    do i = 1, n
        parameters(i) = 0.5_dp + 1.5_dp* &
            real(mod(17*i, 4093), dp)/4092.0_dp
        directions(i) = -0.8_dp + 1.6_dp* &
            real(mod(23*i, 4091), dp)/4090.0_dp
    end do

    if (trim(residency) == "resident") then
        !$acc data copyin(parameters, directions) create(roots, tangents, reliable)
        !$omp target data map(to: parameters, directions) &
        !$omp& map(alloc: roots, tangents, reliable)
        call run_benchmark()
        !$omp target update from(roots, tangents, reliable)
        !$omp end target data
        !$acc update self(roots, tangents, reliable)
        !$acc end data
    else
        call run_benchmark()
    end if
    sink = roots(1) + tangents(n)
    if (sink /= sink .or. .not. all(reliable)) then
        error stop "implicit root benchmark failed"
    end if

contains

    subroutine run_benchmark()
        sink = run_sample()
        if (peak_only) then
            write (*, "(i0)") peak_rss_bytes()
            return
        end if
        do sample = 1, 3
            sink = run_sample()
        end do
        do sample = 1, samples
            write (*, "(f0.6)") run_sample()
        end do
    end subroutine run_benchmark

    function run_sample() result(milliseconds)
        real(dp) :: milliseconds
        integer(int64) :: tick0, tick1, rate

        call system_clock(tick0, rate)
        if (trim(residency) == "transfer") then
            !$acc data copyin(parameters, directions) copyout(roots, tangents, reliable)
            !$omp target data map(to: parameters, directions) &
            !$omp& map(from: roots, tangents, reliable)
            call implicit_root_jvp_batch( &
                n, parameters, directions, roots, tangents, reliable)
            !$omp end target data
            !$acc end data
        else
            call implicit_root_jvp_batch( &
                n, parameters, directions, roots, tangents, reliable)
        end if
        call system_clock(tick1)
        milliseconds = real(tick1 - tick0, dp)*1.0e3_dp/real(rate, dp)
    end function run_sample

end program bench_implicit_root_gpu
