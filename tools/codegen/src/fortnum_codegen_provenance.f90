module fortnum_codegen_provenance
    implicit none
    private

    public :: fortsym_revision

contains

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

end module fortnum_codegen_provenance
