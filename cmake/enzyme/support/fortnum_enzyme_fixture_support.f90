module fortnum_enzyme_fixture_support
    use, intrinsic :: iso_c_binding, only: c_int64_t
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, output_unit
    implicit none
    private

    integer, parameter, public :: fixture_warmup_count = 3
    integer, parameter, public :: fixture_sample_count = 15

    type, public :: fixture_timer_t
        integer(int64), private :: started = 0_int64
        integer(int64), private :: rate = 0_int64
    contains
        procedure :: start => fixture_timer_start
        procedure :: elapsed_ns => fixture_timer_elapsed_ns
    end type fixture_timer_t

    abstract interface
        function fixture_measurement() result(value)
            import dp
            real(dp) :: value
        end function fixture_measurement
    end interface

    public :: collect_fixture_samples
    public :: fixture_peak_rss_bytes
    public :: median_mad
    public :: read_fixture_environment
    public :: read_fixture_integer
    public :: write_fixture_result

    interface
        function c_peak_rss_bytes() bind(c, name="fortnum_peak_rss_bytes") &
                result(bytes)
            import c_int64_t
            integer(c_int64_t) :: bytes
        end function c_peak_rss_bytes
    end interface

contains

    subroutine read_fixture_environment(name, default_value, value)
        character(*), intent(in) :: name, default_value
        character(*), intent(out) :: value
        integer :: status

        value = default_value
        call get_environment_variable(name, value, status=status)
        if (status /= 0) value = default_value
    end subroutine read_fixture_environment

    subroutine read_fixture_integer(name, default_value, value, valid)
        character(*), intent(in) :: name
        integer, intent(in) :: default_value
        integer, intent(out) :: value
        logical, intent(out) :: valid
        character(64) :: text
        integer :: ios, status

        value = default_value
        valid = .true.
        call get_environment_variable(name, text, status=status)
        if (status /= 0) return
        read (text, *, iostat=ios) value
        if (ios /= 0) then
            value = default_value
            valid = .false.
        end if
    end subroutine read_fixture_integer

    subroutine fixture_timer_start(self)
        class(fixture_timer_t), intent(inout) :: self

        call system_clock(self%started, self%rate)
    end subroutine fixture_timer_start

    function fixture_timer_elapsed_ns(self) result(nanoseconds)
        class(fixture_timer_t), intent(in) :: self
        real(dp) :: nanoseconds
        integer(int64) :: finished

        call system_clock(finished)
        if (self%rate <= 0_int64) error stop "fixture timer was not started"
        nanoseconds = real(finished - self%started, dp)*1.0e9_dp/ &
            real(self%rate, dp)
    end function fixture_timer_elapsed_ns

    subroutine collect_fixture_samples(measure, samples, warmup_count)
        procedure(fixture_measurement) :: measure
        real(dp), intent(out) :: samples(:)
        integer, intent(in), optional :: warmup_count
        integer :: count, i
        real(dp) :: discarded

        count = fixture_warmup_count
        if (present(warmup_count)) count = warmup_count
        if (count < 0) error stop "warmup count must be nonnegative"
        discarded = 0.0_dp
        do i = 1, count
            discarded = measure()
        end do
        do i = 1, size(samples)
            samples(i) = measure()
        end do
        if (discarded /= discarded) error stop "fixture warmup produced NaN"
    end subroutine collect_fixture_samples

    subroutine median_mad(values, median, mad)
        real(dp), intent(in) :: values(:)
        real(dp), intent(out) :: median, mad
        real(dp), allocatable :: ordered(:), deviations(:)

        if (size(values) < 1) error stop "median requires at least one value"
        ordered = values
        call sort_values(ordered)
        median = middle_value(ordered)
        deviations = abs(values - median)
        call sort_values(deviations)
        mad = middle_value(deviations)
    end subroutine median_mad

    subroutine write_fixture_result(name, repetitions, median, mad, unit)
        character(*), intent(in) :: name
        integer, intent(in) :: repetitions
        real(dp), intent(in) :: median, mad
        integer, intent(in), optional :: unit
        integer :: destination

        destination = output_unit
        if (present(unit)) destination = unit
        write (destination, "(a,',',i0,',',es24.16,',',es24.16)") &
            trim(name), repetitions, median, mad
    end subroutine write_fixture_result

    function fixture_peak_rss_bytes() result(bytes)
        integer(int64) :: bytes

        bytes = int(c_peak_rss_bytes(), int64)
    end function fixture_peak_rss_bytes

    function middle_value(ordered) result(value)
        real(dp), intent(in) :: ordered(:)
        real(dp) :: value
        integer :: midpoint

        midpoint = size(ordered)/2
        if (mod(size(ordered), 2) == 0) then
            value = 0.5_dp*(ordered(midpoint) + ordered(midpoint + 1))
        else
            value = ordered(midpoint + 1)
        end if
    end function middle_value

    subroutine sort_values(values)
        real(dp), intent(inout) :: values(:)
        real(dp) :: key
        integer :: i, j

        do i = 2, size(values)
            key = values(i)
            j = i - 1
            do while (j >= 1)
                if (values(j) <= key) exit
                values(j + 1) = values(j)
                j = j - 1
            end do
            values(j + 1) = key
        end do
    end subroutine sort_values

end module fortnum_enzyme_fixture_support
