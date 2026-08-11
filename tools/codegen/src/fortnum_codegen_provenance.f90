module fortnum_codegen_provenance
    use, intrinsic :: iso_fortran_env, only: int64
    use fortsym_string, only: chars
    use fortsym_arena, only: arena_t, NK_FUNC
    use fortsym_expr, only: expr_t
    use fortsym_kernel, only: operation_count_t, count_operations
    implicit none
    private

    public :: codegen_log, codegen_log_count
    public :: fortsym_revision, generated_path
    public :: cost_block_text, insert_cost_block

    character(len=1), parameter :: NL = achar(10)

contains

    subroutine codegen_log(message)
        character(*), intent(in) :: message

        if (codegen_verbose()) print "(a)", message
    end subroutine codegen_log

    subroutine codegen_log_count(label, count)
        character(*), intent(in) :: label
        integer, intent(in) :: count

        if (codegen_verbose()) print "(a,i0)", label, count
    end subroutine codegen_log_count

    function codegen_verbose() result(verbose)
        logical :: verbose
        character(16) :: value
        integer :: status

        verbose = .false.
        call get_environment_variable("FORTNUM_CODEGEN_VERBOSE", value, &
            status=status)
        if (status /= 0) return
        select case (trim(value))
        case ("1", "true", "TRUE", "on", "ON", "yes", "YES")
            verbose = .true.
        end select
    end function codegen_verbose

    function fortsym_revision() result(revision)
        character(:), allocatable :: revision
        character(128) :: line
        integer :: unit, ios

        open (newunit=unit, file="fortsym.lock", status="old", action="read", &
            iostat=ios)
        if (ios /= 0) error stop "cannot read fortsym.lock"
        read (unit, "(a)", iostat=ios) line
        close (unit)
        if (ios /= 0) error stop "cannot parse fortsym.lock"
        if (len_trim(line) /= 40) error stop "fortsym.lock must contain a full SHA"
        revision = "fortsym@"//trim(line)
    end function fortsym_revision

    function generated_path(filename) result(path)
        character(*), intent(in) :: filename
        character(:), allocatable :: path
        character(4096) :: output_directory
        integer :: length, status

        call get_environment_variable("FORTNUM_CODEGEN_OUTPUT_DIR", &
            output_directory, length=length, status=status)
        if (status /= 0 .or. length == 0) then
            path = "../../src/generated/"//filename
            return
        end if
        if (output_directory(length:length) == "/") then
            path = output_directory(:length)//filename
        else
            path = output_directory(:length)//"/"//filename
        end if
    end function generated_path

    !> The four-count cost block that #79 adds to each generated source's
    !> manifest. n_sym counts the symbolic DAG as the generator first built
    !> it; n_emit counts the emitted (simplified, post-CSE) roots. The gap is
    !> what the generator's own simplification removed, stored explicitly so a
    !> commit that changes it is reviewable. n_machine and t_meas are absent:
    !> neither was measured here, and an absent level is not a zero.
    !>
    !> The text is a JSON object comment block, emitted byte-for-byte so
    !> check_generated.sh stays a plain `cmp`.
    function cost_block_text(symbolic_roots, emitted_roots) result(text)
        type(expr_t), intent(in) :: symbolic_roots(:), emitted_roots(:)
        character(:), allocatable :: text
        type(operation_count_t) :: sym_ops, emit_ops
        character(len=32), allocatable :: func_names(:)
        integer, allocatable :: func_counts(:)
        character(len=16) :: buffer
        integer :: nfunc, k, gap

        sym_ops = count_operations(symbolic_roots)
        emit_ops = count_operations(emitted_roots)
        gap = emit_ops%total - sym_ops%total

        call transcendental_breakdown(symbolic_roots, func_names, func_counts, nfunc)

        text = "! cost: {"//NL
        text = text//'!   "n_sym": {"flops": '//integer_text(sym_ops%total)// &
            ', "adds": '//integer_text(sym_ops%additions)// &
            ', "muls": '//integer_text(sym_ops%multiplications)// &
            ', "divs": '//integer_text(sym_ops%divisions)// &
            ', "transcendental": '
        if (nfunc == 0) then
            text = text//"{}"//","//NL
        else
            text = text//"{"
            do k = 1, nfunc
                write (buffer, "(i0)") func_counts(k)
                if (k > 1) text = text//", "
                text = text//'"'//trim(func_names(k))//'": '//trim(buffer)
            end do
            text = text//"}"//","//NL
        end if
        text = text//'!   "n_emit": {"flops": '//integer_text(emit_ops%total)// &
            ', "instructions": '//integer_text(emit_ops%total)//"}"//","//NL
        text = text//'!   "gaps": {"generator": '//integer_text(gap)//"}"//NL
        text = text//"! }"//NL
    end function cost_block_text

    !> Insert a cost-block comment block into an emit_kernel banner, directly
    !> after the "! Regenerate with:" line so the manifest stays one block.
    function insert_cost_block(code, cost_text) result(new_code)
        character(*), intent(in) :: code, cost_text
        character(:), allocatable :: new_code
        integer :: pos, line_end

        pos = index(code, "! Regenerate with:")
        if (pos == 0) then
            new_code = code
            return
        end if
        line_end = index(code(pos:), NL)
        if (line_end == 0) then
            new_code = code
            return
        end if
        line_end = pos + line_end - 1
        new_code = code(:line_end)//cost_text//code(line_end + 1:)
    end function insert_cost_block

    !> Distinct transcendental function nodes across the roots, with counts,
    !> in first-seen order. Hash-consing means a shared node is counted once,
    !> matching count_operations.
    subroutine transcendental_breakdown(roots, names, counts, n)
        type(expr_t), intent(in) :: roots(:)
        character(len=32), allocatable, intent(out) :: names(:)
        integer, allocatable, intent(out) :: counts(:)
        integer, intent(out) :: n
        type(arena_t), pointer :: a
        logical, allocatable :: visited(:)
        integer :: k

        n = 0
        if (size(roots) == 0) then
            allocate (names(0), counts(0))
            return
        end if
        a => roots(1)%a
        allocate (visited(a%size()), source=.false.)
        allocate (names(a%size()), counts(a%size()))
        do k = 1, size(roots)
            call walk(a, roots(k)%id, visited, names, counts, n)
        end do
    end subroutine transcendental_breakdown

    recursive subroutine walk(a, id, visited, names, counts, n)
        type(arena_t), intent(in) :: a
        integer, intent(in) :: id
        logical, intent(inout) :: visited(:)
        character(len=32), intent(inout) :: names(:)
        integer, intent(inout) :: counts(:)
        integer, intent(inout) :: n
        integer :: k, j
        character(len=32) :: fname

        if (visited(id)) return
        visited(id) = .true.
        do k = 1, a%nargs_of(id)
            call walk(a, a%arg_of(id, k), visited, names, counts, n)
        end do
        if (a%kind_of(id) == NK_FUNC) then
            fname = trim(chars(a%name_of(id)))
            if (is_transcendental(fname)) then
                do j = 1, n
                    if (trim(names(j)) == trim(fname)) then
                        counts(j) = counts(j) + 1
                        return
                    end if
                end do
                n = n + 1
                names(n) = fname
                counts(n) = 1
            end if
        end if
    end subroutine walk

    !> True for the transcendental function calls fortsym emits. Array-indexed
    !> symbols are also NK_FUNC nodes but carry arbitrary names ("x", "v",
    !> "u") and are not arithmetic cost; sqrt and abs are algebraic helpers
    !> rather than transcendental functions.
    pure function is_transcendental(name) result(yes)
        character(*), intent(in) :: name
        logical :: yes
        select case (trim(name))
        case ("sin", "cos", "tan", "asin", "acos", "atan", "atan2", &
              "sinh", "cosh", "tanh", "asinh", "acosh", "atanh", &
              "exp", "log", "erf", "erfc", "gamma", "besselj", &
              "legendrep", "legendreq")
            yes = .true.
        case default
            yes = .false.
        end select
    end function is_transcendental

    pure function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text
        character(len=16) :: buffer

        write (buffer, "(i0)") value
        text = trim(buffer)
    end function integer_text

end module fortnum_codegen_provenance
