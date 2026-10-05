program bench_certified_reductions
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_certified_reductions, only: reduction_t, reduction_ok, &
        aggregate_sum, aggregate_dot, interval_sum, interval_dot
    use fortnum_reduction_fixtures, only: validate_reductions, &
        make_reduction_workload, require_encloses
    use fortnum_interval, only: interval_t, interval, exp
    implicit none

    integer, parameter :: samples = 15, warmups = 3
    integer :: n, repetitions, sample, i, repetition, ios
    integer(int64) :: start, finish, rate
    real(dp), allocatable :: x(:), y(:)
    type(interval_t), allocatable :: values(:)
    type(reduction_t) :: result
    real(dp) :: sink, ns
    real(qp) :: exact_sum, exact_dot, exact
    character(32) :: candidate, product, pattern, argument

    call get_command_argument(1, candidate)
    call validate_reductions()
    if (trim(candidate) == "validate") then
        print '(a)', "certified reduction fixtures PASS"
        stop
    end if
    call get_command_argument(2, product)
    call get_command_argument(3, pattern)
    n = 4096
    repetitions = 100
    if (command_argument_count() >= 4) then
        call get_command_argument(4, argument)
        read (argument, *, iostat=ios) n
        if (ios /= 0) error stop "invalid length"
    end if
    if (command_argument_count() >= 5) then
        call get_command_argument(5, argument)
        read (argument, *, iostat=ios) repetitions
        if (ios /= 0) error stop "invalid repetitions"
    end if
    if (n < 1 .or. n > 1048576 .or. repetitions < 1) error stop "invalid workload"
    if (trim(candidate) /= "aggregate" .and. trim(candidate) /= "interval") &
        error stop "candidate must be aggregate or interval"
    if (trim(product) /= "sum" .and. trim(product) /= "dot" .and. &
        trim(product) /= "exp") error stop "product must be sum, dot or exp"
    if (trim(product) == "exp") then
        if (trim(candidate) /= "interval" .or. trim(pattern) /= "negative") &
            error stop "exp baseline requires interval exp negative"
    else
        if (trim(pattern) /= "positive" .and. trim(pattern) /= "cancel") &
            error stop "pattern must be positive or cancel"
    end if
    allocate (x(n), y(n))
    call make_reduction_workload(x, y, trim(pattern) == "cancel", exact_sum, exact_dot)
    if (trim(product) == "exp") then
        allocate (values(n))
        ! Bounded negative arguments representative of Gaussian-weight costs.
        do i = 1, n
            x(i) = -real(mod(i, 2049), dp)/1024.0_dp
        end do
        call evaluate_exp()
        call validate_exp()
    else
        call evaluate_reduction()
        exact = exact_sum
        if (trim(product) == "dot") exact = exact_dot
        call require_encloses(result, exact)
        print '(a,es24.16)', "# center=", result%center
        print '(a,es24.16)', "# radius=", result%radius
    end if
    print '(a)', "# wall ns/element; 3 warmups, 15 raw samples; validation outside timing"
    print '(a)', "candidate,product,pattern,length,repetitions,sample,ns_per_element,sink"
    do sample = 1 - warmups, samples
        sink = 0.0_dp
        call system_clock(start, rate)
        do repetition = 1, repetitions
            if (trim(product) == "exp") then
                call evaluate_exp()
                sink = sink + values(1)%lo + values(n)%hi
            else
                call evaluate_reduction()
                sink = sink + result%center + result%radius
            end if
        end do
        call system_clock(finish)
        if (.not. ieee_is_finite(sink)) error stop "nonfinite benchmark sink"
        if (trim(product) == "exp") then
            call validate_exp()
        else
            if (result%status /= reduction_ok) error stop "timed reduction failed"
            call require_encloses(result, exact)
        end if
        ns = real(finish - start, dp)*1.0e9_dp &
            / (real(rate, dp)*real(n, dp)*real(repetitions, dp))
        if (sample > 0) then
            write (*, '(a,",",a,",",a,3(",",i0),2(",",es24.16))') &
                trim(candidate), trim(product), trim(pattern), n, repetitions, sample, &
                ns, sink
        end if
    end do

contains

    subroutine evaluate_reduction()
        if (trim(candidate) == "aggregate") then
            if (trim(product) == "sum") then
                result = aggregate_sum(x)
            else
                result = aggregate_dot(x, y)
            end if
        else
            if (trim(product) == "sum") then
                result = interval_sum(x)
            else
                result = interval_dot(x, y)
            end if
        end if
    end subroutine evaluate_reduction

    subroutine evaluate_exp()
        integer :: j

        do j = 1, n
            values(j) = exp(interval(x(j)))
        end do
    end subroutine evaluate_exp

    subroutine validate_exp()
        type(interval_t) :: at_zero

        if (.not. all(ieee_is_finite(values%lo))) error stop "nonfinite exp lower"
        if (.not. all(ieee_is_finite(values%hi))) error stop "nonfinite exp upper"
        if (any(values%lo < 0.0_dp) .or. any(values%hi > 1.01_dp)) &
            error stop "exp baseline outside expected bounded range"
        at_zero = exp(interval(0.0_dp))
        if (at_zero%lo > 1.0_dp .or. at_zero%hi < 1.0_dp) &
            error stop "exp baseline misses exp(0)=1"
        ! This row measures the existing interval primitive, not a new exp
        ! candidate. Its general enclosure oracle remains test_fortnum_interval.
    end subroutine validate_exp

end program bench_certified_reductions
