module fortnum_codegen_provenance
    implicit none
    private

    public :: codegen_log, codegen_log_count
    public :: fortsym_revision, generated_path

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

end module fortnum_codegen_provenance
