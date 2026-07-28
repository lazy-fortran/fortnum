program test_enzyme_fixture_support
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_enzyme_fixture_support, only: fixture_peak_rss_bytes, &
        fixture_sample_count, fixture_timer_t, fixture_warmup_count, &
        collect_fixture_samples, median_mad, read_fixture_environment, &
        read_fixture_integer, write_fixture_result, write_fixture_scaling_result
    implicit none

    type(fixture_timer_t) :: timer
    real(dp) :: median, mad
    real(dp) :: samples(4)
    character(32) :: text
    integer :: integer_value, calls, unit
    logical :: valid
    character(256) :: line

    if (fixture_warmup_count /= 3) error stop "unexpected warmup default"
    if (fixture_sample_count /= 15) error stop "unexpected sample default"

    call read_fixture_environment( &
        "FORTNUM_FIXTURE_TEST_TEXT", "missing", text)
    if (trim(text) /= "configured") error stop "environment text mismatch"
    call read_fixture_environment( &
        "FORTNUM_FIXTURE_TEST_MISSING", "fallback", text)
    if (trim(text) /= "fallback") error stop "environment fallback mismatch"

    call read_fixture_integer( &
        "FORTNUM_FIXTURE_TEST_INTEGER", 7, integer_value, valid)
    if (.not. valid .or. integer_value /= 42) then
        error stop "environment integer mismatch"
    end if
    call read_fixture_integer( &
        "FORTNUM_FIXTURE_TEST_INVALID", 7, integer_value, valid)
    if (valid .or. integer_value /= 7) then
        error stop "invalid integer fallback mismatch"
    end if

    calls = 0
    call collect_fixture_samples(measure, samples, 2)
    if (calls /= 6) error stop "warmup calls were not discarded"
    if (any(samples /= [3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp])) then
        error stop "sample collection mismatch"
    end if

    call median_mad([9.0_dp, 1.0_dp, 5.0_dp, 3.0_dp], median, mad)
    if (median /= 4.0_dp .or. mad /= 2.0_dp) then
        error stop "even median/MAD mismatch"
    end if
    call median_mad([9.0_dp, 1.0_dp, 5.0_dp, 3.0_dp, 7.0_dp], median, mad)
    if (median /= 5.0_dp .or. mad /= 2.0_dp) then
        error stop "odd median/MAD mismatch"
    end if

    open (newunit=unit, status="scratch", action="readwrite")
    call write_fixture_result("candidate", 11, 2.5_dp, 0.25_dp, unit)
    rewind (unit)
    read (unit, "(a)") line
    close (unit)
    if (trim(line) /= &
        "candidate,11,  2.5000000000000000E+00,  2.5000000000000000E-01") then
        error stop "standard result format mismatch"
    end if

    open (newunit=unit, status="scratch", action="readwrite")
    call write_fixture_scaling_result( &
        "candidate", 4, 11, 2.5_dp, 0.25_dp, unit)
    rewind (unit)
    read (unit, "(a)") line
    close (unit)
    if (trim(line) /= &
        "candidate,4,11,  2.5000000000000000E+00,  2.5000000000000000E-01") then
        error stop "scaling result format mismatch"
    end if

    call timer%start()
    if (timer%elapsed_ns() < 0.0_dp) error stop "timer moved backwards"
    if (fixture_peak_rss_bytes() <= 0_int64) error stop "peak RSS unavailable"

contains

    function measure() result(value)
        real(dp) :: value

        calls = calls + 1
        value = real(calls, dp)
    end function measure

end program test_enzyme_fixture_support
