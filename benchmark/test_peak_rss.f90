program test_peak_rss
    use, intrinsic :: iso_c_binding, only: c_int8_t, c_int64_t
    use fortnum_benchmark_memory, only: peak_rss_bytes
    implicit none

    integer, parameter :: allocation_bytes = 16*1024*1024
    integer(c_int8_t), allocatable :: buffer(:)
    integer(c_int64_t) :: before, after, checksum

    before = peak_rss_bytes()
    allocate (buffer(allocation_bytes))
    buffer = 1_c_int8_t
    checksum = sum(int(buffer, c_int64_t))
    after = peak_rss_bytes()

    if (before <= 0_c_int64_t) error stop "peak RSS query failed"
    if (checksum /= int(allocation_bytes, c_int64_t)) then
        error stop "allocation was not touched"
    end if
    if (after - before < 8_c_int64_t*1024_c_int64_t*1024_c_int64_t) then
        error stop "peak RSS did not observe the known allocation"
    end if
end program test_peak_rss
