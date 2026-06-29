module fortnum_kinds
    use, intrinsic :: iso_fortran_env, only: dp => real64, sp => real32, &
        int32, int64
    implicit none
    private

    ! Floating-point kinds
    public :: dp ! double precision (real64) – primary working kind
    public :: sp ! single precision (real32) – for mixed-precision interfaces

    ! Integer kinds
    public :: i4 ! 32-bit integer
    public :: i8 ! 64-bit integer – for large index spaces and counts

    integer, parameter :: i4 = int32
    integer, parameter :: i8 = int64

end module fortnum_kinds
