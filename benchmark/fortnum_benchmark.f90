module fortnum_benchmark
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, output_unit
    implicit none
    private

    ! Extensible timing harness for fortnum routines.  A benchmark is a
    ! (name, callable) pair: the callable runs one unit of work and returns a
    ! real(dp) result.  Returning the result lets the driver consume it so the
    ! optimiser cannot elide the call, which would otherwise report fake speed.
    abstract interface
        function bench_kernel() result(sink)
            import :: dp
            real(dp) :: sink
        end function bench_kernel
    end interface

    ! Reverse-mode / forward-mode derivative kernels land with #36 (M6).  The
    ! interface mirrors bench_kernel so a JVP/VJP callable plugs into the same
    ! registry; the driver then reports the primal-vs-derivative ratio.  Kept
    ! here so the slot is reserved without forcing a redesign later.
    abstract interface
        function bench_deriv_kernel() result(sink)
            import :: dp
            real(dp) :: sink
        end function bench_deriv_kernel
    end interface

    integer, parameter :: NAME_LEN = 48

    type, public :: bench_case_t
        character(NAME_LEN)                  :: name = ""
        procedure(bench_kernel), pointer, nopass :: run => null()
        ! Optional derivative counterpart for the M6 primal/derivative ratio.
        ! Null until #36 supplies a JVP/VJP kernel for this case.
        procedure(bench_deriv_kernel), pointer, nopass :: run_deriv => null()
    end type bench_case_t

    type, public :: bench_result_t
        character(NAME_LEN) :: name = ""
        integer(int64)      :: reps = 0_int64
        real(dp)            :: ns_per_call = 0.0_dp
        logical             :: has_deriv = .false.
        real(dp)            :: deriv_ns_per_call = 0.0_dp
        real(dp)            :: deriv_ratio = 0.0_dp
    end type bench_result_t

    type, public :: bench_registry_t
        integer                          :: count = 0
        type(bench_case_t), allocatable  :: cases(:)
    contains
        procedure :: add  => registry_add
        procedure :: run  => registry_run
    end type bench_registry_t

    public :: bench_kernel
    public :: bench_deriv_kernel
    public :: time_kernel
    public :: report_header
    public :: report_result

    integer, parameter :: REGISTRY_CHUNK = 8

contains

    ! Registers a benchmark.  deriv is optional; pass it once an autodiff
    ! kernel exists so the same case reports the primal/derivative overhead.
    subroutine registry_add(self, name, run, deriv)
        class(bench_registry_t), intent(inout) :: self
        character(*),            intent(in)     :: name
        procedure(bench_kernel)                 :: run
        procedure(bench_deriv_kernel), optional :: deriv
        type(bench_case_t), allocatable :: grown(:)

        if (.not. allocated(self%cases)) then
            allocate (self%cases(REGISTRY_CHUNK))
        else if (self%count == size(self%cases)) then
            allocate (grown(size(self%cases) + REGISTRY_CHUNK))
            grown(1:self%count) = self%cases(1:self%count)
            call move_alloc(grown, self%cases)
        end if

        self%count = self%count + 1
        self%cases(self%count)%name = name
        self%cases(self%count)%run => run
        if (present(deriv)) self%cases(self%count)%run_deriv => deriv
    end subroutine registry_add

    ! Runs every registered case for reps repetitions and prints a table.
    subroutine registry_run(self, reps)
        class(bench_registry_t), intent(in) :: self
        integer(int64),          intent(in) :: reps
        type(bench_result_t) :: res
        integer :: i

        call report_header()
        do i = 1, self%count
            res%name = self%cases(i)%name
            res%reps = reps
            res%ns_per_call = time_kernel(self%cases(i)%run, reps)
            res%has_deriv = associated(self%cases(i)%run_deriv)
            if (res%has_deriv) then
                res%deriv_ns_per_call = &
                    time_kernel(self%cases(i)%run_deriv, reps)
                if (res%ns_per_call > 0.0_dp) then
                    res%deriv_ratio = res%deriv_ns_per_call / res%ns_per_call
                end if
            end if
            call report_result(res)
        end do
    end subroutine registry_run

    ! Times one kernel over reps calls and returns nanoseconds per call.  A
    ! running accumulator consumes every result so the call cannot be hoisted
    ! out of the loop or constant-folded away.
    function time_kernel(run, reps) result(ns_per_call)
        procedure(bench_kernel)    :: run
        integer(int64), intent(in) :: reps
        real(dp)                   :: ns_per_call
        integer(int64) :: tick0, tick1, rate
        real(dp)       :: acc, elapsed_ns
        integer(int64) :: k

        acc = 0.0_dp
        call system_clock(count=tick0, count_rate=rate)
        do k = 1_int64, reps
            acc = acc + run()
        end do
        call system_clock(count=tick1)

        ! Defeat dead-store elimination on acc without printing noise: a NaN
        ! guard touches the value but is unreachable for finite kernels.
        if (acc /= acc) then
            write (output_unit, '(a)') "benchmark accumulator went NaN"
        end if

        elapsed_ns = real(tick1 - tick0, dp) / real(rate, dp) * 1.0e9_dp
        if (reps > 0_int64) then
            ns_per_call = elapsed_ns / real(reps, dp)
        else
            ns_per_call = 0.0_dp
        end if
    end function time_kernel

    subroutine report_header()
        write (output_unit, '(a)') "fortnum micro-benchmark"
        write (output_unit, '(a48,1x,a12,1x,a12,1x,a12,1x,a10)') &
            "name", "reps", "ns/call", "deriv ns", "ratio"
        write (output_unit, '(a)') repeat('-', 48 + 1 + 12 + 1 + 12 + 1 + 12 + 1 + 10)
    end subroutine report_header

    subroutine report_result(res)
        type(bench_result_t), intent(in) :: res
        if (res%has_deriv) then
            write (output_unit, '(a48,1x,i12,1x,f12.3,1x,f12.3,1x,f10.2)') &
                res%name, res%reps, res%ns_per_call, &
                res%deriv_ns_per_call, res%deriv_ratio
        else
            write (output_unit, '(a48,1x,i12,1x,f12.3,1x,a12,1x,a10)') &
                res%name, res%reps, res%ns_per_call, "-", "-"
        end if
    end subroutine report_result

end module fortnum_benchmark
