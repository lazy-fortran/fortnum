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

    ! Tag for how a derivative is produced once M6 lands: analytic hand-code,
    ! implicit-function differentiation, tape/trace replay, or generated
    ! (Enzyme/flang).  "primal" marks a case with no derivative kernel yet.
    integer, parameter :: BACKEND_LEN = 16

    type, public :: bench_result_t
        character(NAME_LEN)    :: name = ""
        integer(int64)         :: reps = 0_int64
        real(dp)               :: ns_per_call = 0.0_dp
        character(BACKEND_LEN) :: backend = "primal"
        logical                :: has_deriv = .false.
        real(dp)               :: deriv_ns_per_call = 0.0_dp
        ! Derivative-product overhead ratios (deriv ns / primal ns).  All zero
        ! until #36 wires JVP/VJP/grad/HVP kernels; the schema reserves a slot
        ! for each so gate.py can track them without a later format change.
        real(dp)               :: jvp_primal = 0.0_dp
        real(dp)               :: vjp_primal = 0.0_dp
        real(dp)               :: grad_primal = 0.0_dp
        real(dp)               :: hvp_primal = 0.0_dp
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
    public :: report_json

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

    ! Runs every registered case for reps repetitions and reports the table.
    ! as_json selects machine-readable JSON over the human table; gate.py
    ! consumes the JSON form.
    subroutine registry_run(self, reps, as_json)
        class(bench_registry_t), intent(in) :: self
        integer(int64),          intent(in) :: reps
        logical, optional,       intent(in) :: as_json
        type(bench_result_t) :: res
        integer :: i
        logical :: emit_json

        emit_json = .false.
        if (present(as_json)) emit_json = as_json

        if (emit_json) then
            write (output_unit, '(a)') "{"
            write (output_unit, '(a)') '  "benchmarks": ['
        else
            call report_header()
        end if

        do i = 1, self%count
            res%name = self%cases(i)%name
            res%reps = reps
            res%ns_per_call = time_kernel(self%cases(i)%run, reps)
            res%has_deriv = associated(self%cases(i)%run_deriv)
            if (res%has_deriv) then
                res%deriv_ns_per_call = &
                    time_kernel(self%cases(i)%run_deriv, reps)
                ! Until M6 distinguishes JVP/VJP/grad/HVP kernels, the single
                ! derivative kernel reports as the JVP/primal overhead.
                if (res%ns_per_call > 0.0_dp) then
                    res%jvp_primal = res%deriv_ns_per_call / res%ns_per_call
                end if
            end if
            if (emit_json) then
                call report_json(res, i == self%count)
            else
                call report_result(res)
            end if
        end do

        if (emit_json) then
            write (output_unit, '(a)') "  ]"
            write (output_unit, '(a)') "}"
        end if
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
                res%deriv_ns_per_call, res%jvp_primal
        else
            write (output_unit, '(a48,1x,i12,1x,f12.3,1x,a12,1x,a10)') &
                res%name, res%reps, res%ns_per_call, "-", "-"
        end if
    end subroutine report_result

    ! Emits one benchmark object in the machine-readable JSON consumed by
    ! gate.py.  last suppresses the trailing comma on the final element.  The
    ! derivative-product fields are present but null when no derivative kernel
    ! exists, so the schema is stable from M5 through M6.
    subroutine report_json(res, last)
        type(bench_result_t), intent(in) :: res
        logical,              intent(in) :: last
        character(:), allocatable :: tail

        if (last) then
            tail = ""
        else
            tail = ","
        end if

        write (output_unit, '(a)') "    {"
        write (output_unit, '(3a)') &
            '      "name": "', trim(res%name), '",'
        write (output_unit, '(a,i0,a)') &
            '      "reps": ', res%reps, ','
        write (output_unit, '(a,f0.6,a)') &
            '      "ns_per_call": ', res%ns_per_call, ','
        write (output_unit, '(3a)') &
            '      "backend": "', trim(res%backend), '",'
        write (output_unit, '(2a)') &
            '      "deriv_ns_per_call": ', deriv_field(res%has_deriv, res%deriv_ns_per_call)
        write (output_unit, '(3a)') &
            '      "jvp_primal": ', trim(ratio_field(res%has_deriv, res%jvp_primal)), ','
        write (output_unit, '(3a)') &
            '      "vjp_primal": ', trim(ratio_field(.false., res%vjp_primal)), ','
        write (output_unit, '(3a)') &
            '      "grad_primal": ', trim(ratio_field(.false., res%grad_primal)), ','
        write (output_unit, '(3a)') &
            '      "hvp_primal": ', trim(ratio_field(.false., res%hvp_primal))
        write (output_unit, '(2a)') "    }", tail
    end subroutine report_json

    ! Renders a derivative-ns value as a JSON number, or null when no
    ! derivative kernel ran.  Trailing comma included (this is the last
    ! comma-terminated field is handled by the caller).
    function deriv_field(has_deriv, value) result(text)
        logical,  intent(in) :: has_deriv
        real(dp), intent(in) :: value
        character(:), allocatable :: text
        character(32) :: buf

        if (has_deriv) then
            write (buf, '(f0.6)') value
            text = trim(buf)//","
        else
            text = "null,"
        end if
    end function deriv_field

    ! Renders an overhead ratio as a JSON number, or null when the
    ! corresponding derivative product has no kernel yet.
    function ratio_field(present_flag, value) result(text)
        logical,  intent(in) :: present_flag
        real(dp), intent(in) :: value
        character(:), allocatable :: text
        character(32) :: buf

        if (present_flag) then
            write (buf, '(f0.6)') value
            text = trim(buf)
        else
            text = "null"
        end if
    end function ratio_field

end module fortnum_benchmark
