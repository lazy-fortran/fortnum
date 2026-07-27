module fortnum_benchmark_memory
    use, intrinsic :: iso_c_binding, only: c_int64_t
    implicit none
    private

    public :: peak_rss_bytes

    interface
        function fortnum_peak_rss_bytes() bind(c) result(bytes)
            import :: c_int64_t
            integer(c_int64_t) :: bytes
        end function fortnum_peak_rss_bytes
    end interface

contains

    function peak_rss_bytes() result(bytes)
        integer(c_int64_t) :: bytes
        bytes = fortnum_peak_rss_bytes()
    end function peak_rss_bytes

end module fortnum_benchmark_memory
