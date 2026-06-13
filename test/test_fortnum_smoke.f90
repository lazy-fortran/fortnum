program test_fortnum_smoke
    use fortnum_version, only: fortnum_version_string
    implicit none

    if (len_trim(fortnum_version_string) == 0) then
        error stop "fortnum_version_string is empty"
    end if

    stop 0
end program test_fortnum_smoke
