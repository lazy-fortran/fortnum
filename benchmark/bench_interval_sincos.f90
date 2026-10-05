program bench_interval_sincos
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_interval, only: interval_t, interval, sin, cos, sincos
    use fortnum_certified_reductions, only: reduction_environment_supported
    implicit none
    integer, parameter :: count = 192, samples = 15, warmups = 3, repetitions = 100
    type(interval_t) :: x(count), s(count), c(count), ref_s(count), ref_c(count)
    integer(int64) :: start, finish, rate
    integer :: i, sample, repetition
    real(dp) :: lower, sink, ns
    real(qp) :: point
    character(16) :: candidate

    call get_command_argument(1, candidate)
    if (trim(candidate) /= 'paired' .and. trim(candidate) /= 'separate' .and. &
        trim(candidate) /= 'validate') error stop 'usage: paired|separate|validate'
    if (.not. reduction_environment_supported()) error stop 'unsupported FP environment'
    ! Exactly represented intervals of width 1/16, over negative/positive
    ! arguments and all quadrants; both candidates see identical fixed inputs.
    do i = 1, count
        lower = real(i - 97, dp)/8.0_dp
        x(i) = interval(lower, lower + 0.0625_dp)
        ref_s(i) = sin(x(i))
        ref_c(i) = cos(x(i))
    end do
    call sincos(x, s, c)
    call validate()
    if (trim(candidate) == 'validate') then
        print '(a)', 'paired interval sincos benchmark validation PASS'
    else
        print '(a)', 'candidate,intervals,repetitions,sample,wall_ns_per_interval,sink'
        do sample = 1 - warmups, samples
            sink = 0.0_dp
            call system_clock(start, rate)
            do repetition = 1, repetitions
                if (trim(candidate) == 'paired') then
                    call sincos(x, s, c)
                else
                    do i = 1, count
                        s(i) = sin(x(i))
                        c(i) = cos(x(i))
                    end do
                end if
                sink = sink + s(1)%lo + c(count)%hi
            end do
            call system_clock(finish)
            call validate()
            if (.not. ieee_is_finite(sink)) error stop 'nonfinite timing sink'
            ns = real(finish - start, dp)*1.0e9_dp &
                / (real(rate, dp)*real(count, dp)*real(repetitions, dp))
            if (sample > 0) write (*, '(a,3(",",i0),2(",",es24.16))') &
                trim(candidate), count, repetitions, sample, ns, sink
        end do
    end if

contains

    subroutine validate()
        integer :: j, p

        do j = 1, count
            if (s(j)%lo < ref_s(j)%lo .or. s(j)%lo > ref_s(j)%lo .or. &
                s(j)%hi < ref_s(j)%hi .or. s(j)%hi > ref_s(j)%hi) &
                error stop 'paired/separate sine enclosure changed'
            if (c(j)%lo < ref_c(j)%lo .or. c(j)%lo > ref_c(j)%lo .or. &
                c(j)%hi < ref_c(j)%hi .or. c(j)%hi > ref_c(j)%hi) &
                error stop 'paired/separate cosine enclosure changed'
            do p = 0, 8
                point = real(x(j)%lo, qp) &
                    + (real(x(j)%hi, qp) - real(x(j)%lo, qp))*real(p, qp)/8.0_qp
                if (real(s(j)%lo, qp) > sin(point) .or. &
                    real(s(j)%hi, qp) < sin(point)) error stop 'sine sample excluded'
                if (real(c(j)%lo, qp) > cos(point) .or. &
                    real(c(j)%hi, qp) < cos(point)) error stop 'cosine sample excluded'
            end do
        end do
    end subroutine validate

end program bench_interval_sincos
