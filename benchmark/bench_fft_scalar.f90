program bench_fft_scalar
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_fft, only: fft_c2c_jvp, fft_c2c_vjp
    implicit none

    integer, parameter :: sample_count = 15
    character(16) :: action, argument, product
    complex(dp), allocatable :: values(:)
    integer(int64) :: repetitions
    integer :: n, sample
    logical :: peak_only
    real(dp) :: warmup

    call read_setting("FORTNUM_FFT_SCALAR_ACTION", "validate", action)
    if (trim(action) == "validate") then
        call validate_length(8)
        call validate_length(30)
        call validate_length(31)
        write (*, "(a)") "PASS"
        stop
    end if
    if (trim(action) /= "benchmark" .and. trim(action) /= "peak-rss") &
        error stop "action must be validate, benchmark, or peak-rss"
    peak_only = trim(action) == "peak-rss"

    call read_setting("FORTNUM_FFT_SCALAR_PRODUCT", "jvp", product)
    if (trim(product) /= "jvp" .and. trim(product) /= "vjp") &
        error stop "product must be jvp or vjp"
    call read_setting("FORTNUM_FFT_SCALAR_LENGTH", "8", argument)
    read (argument, *) n
    if (n < 1) error stop "length must be positive"

    allocate (values(n))
    call fill_values(values, 1_int64)
    repetitions = max(16_int64, 1048576_int64/int(n, int64))
    call read_repetitions(repetitions)
    warmup = run_sample(max(1_int64, repetitions/10_int64))
    if (peak_only) then
        warmup = run_sample(repetitions)
        write (*, "(i0)") peak_rss_bytes()
        stop
    end if
    do sample = 1, 3
        warmup = run_sample(max(1_int64, repetitions/10_int64))
    end do
    do sample = 1, sample_count
        write (*, "(f0.4)") run_sample(repetitions)
    end do

contains

    subroutine read_setting(name, default_value, value)
        character(*), intent(in) :: name, default_value
        character(*), intent(out) :: value
        integer :: length, status

        value = default_value
        call get_environment_variable(name, value, length=length, status=status)
        if (status /= 0 .or. length < 1) value = default_value
    end subroutine read_setting

    subroutine read_repetitions(value)
        integer(int64), intent(inout) :: value
        character(32) :: text
        integer :: ios, length, status

        call get_environment_variable( &
            "FORTNUM_FFT_SCALAR_REPETITIONS", text, length=length, status=status)
        if (status /= 0 .or. length < 1) return
        read (text(:length), *, iostat=ios) value
        if (ios /= 0 .or. value < 1_int64) &
            error stop "repetitions must be positive"
    end subroutine read_repetitions

    function run_sample(count) result(ns_per_call)
        integer(int64), intent(in) :: count
        integer(int64) :: iteration, tick0, tick1, rate
        real(dp) :: ns_per_call, sink

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            call fill_values(values, iteration)
            if (trim(product) == "jvp") then
                call fft_c2c_jvp(values, -1)
            else
                call fft_c2c_vjp(values, -1)
            end if
            sink = sink + real(values(1), dp)
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "scalar FFT benchmark produced NaN"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine validate_length(length)
        integer, intent(in) :: length
        complex(dp), allocatable :: actual(:), expected(:), input(:)
        real(dp), parameter :: tolerance = 4.0e-12_dp
        real(dp) :: scale

        allocate (actual(length), expected(length), input(length))
        call fill_values(input, 3_int64)

        actual = input
        call fft_c2c_jvp(actual, -1)
        call direct_dft(input, -1, expected)
        scale = max(1.0_dp, maxval(abs(expected)))
        if (maxval(abs(actual - expected)) > tolerance*scale) &
            error stop "scalar FFT JVP disagrees with direct DFT"

        actual = input
        call fft_c2c_vjp(actual, -1)
        call direct_dft(input, +1, expected)
        scale = max(1.0_dp, maxval(abs(expected)))
        if (maxval(abs(actual - expected)) > tolerance*scale) &
            error stop "scalar FFT VJP disagrees with direct adjoint DFT"
    end subroutine validate_length

    pure subroutine fill_values(output, iteration)
        complex(dp), intent(out) :: output(:)
        integer(int64), intent(in) :: iteration
        integer :: i

        do i = 1, size(output)
            output(i) = cmplx( &
                sin(0.013_dp*real(i + iteration, dp)), &
                cos(0.017_dp*real(3*i + iteration, dp)), dp)
        end do
    end subroutine fill_values

    pure subroutine direct_dft(input, sign, output)
        complex(dp), intent(in) :: input(:)
        integer, intent(in) :: sign
        complex(dp), intent(out) :: output(size(input))
        real(dp), parameter :: pi = acos(-1.0_dp)
        complex(dp) :: phase
        real(dp) :: angle
        integer :: frequency, sample, length

        length = size(input)
        do frequency = 0, length - 1
            output(frequency + 1) = cmplx(0.0_dp, 0.0_dp, dp)
            do sample = 0, length - 1
                angle = real(sign, dp)*2.0_dp*pi* &
                    real(frequency*sample, dp)/real(length, dp)
                phase = cmplx(cos(angle), sin(angle), dp)
                output(frequency + 1) = output(frequency + 1) + &
                    input(sample + 1)*phase
            end do
        end do
    end subroutine direct_dft

end program bench_fft_scalar
