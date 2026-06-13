program test_fortnum_kinds
    use fortnum_kinds, only: dp, sp, i4, i8
    use, intrinsic :: iso_fortran_env, only: real64, real32, int32, int64
    implicit none

    ! Verify kind parameters resolve to the expected iso_fortran_env values.
    if (dp /= real64)  stop 1
    if (sp /= real32)  stop 1
    if (i4 /= int32)   stop 1
    if (i8 /= int64)   stop 1

    ! Confirm arithmetic uses the declared kind (catches silent kind mismatches).
    block
        real(dp) :: x
        real(sp) :: y
        integer(i4) :: n4
        integer(i8) :: n8
        x  = 1.0_dp
        y  = 1.0_sp
        n4 = 1_i4
        n8 = 1_i8
        if (kind(x) /= dp) stop 1
        if (kind(y) /= sp) stop 1
        if (kind(n4) /= i4) stop 1
        if (kind(n8) /= i8) stop 1
    end block

    stop 0
end program test_fortnum_kinds
