program test_fortnum_special_dawson_oracle
    ! Load the scipy.special.dawsn reference table and verify dawson(x)
    ! agrees within tolerance over the full grid (series / Rybicki /
    ! asymptotic regions, positive and negative arguments).
    !
    ! Tolerances: the Rybicki method targets ~7e-18 sampling error and the
    ! Maclaurin/asymptotic series are carried to ~1e-16 relative accuracy.
    ! We ask for atol = 1e-14, rtol = 1e-13, which gives a safe margin over
    ! the last-ULP scatter between gfortran and SciPy.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_oracle, only: oracle_table_t, oracle_read, oracle_check
    use fortnum_special_dawson, only: dawson
    implicit none

    type(oracle_table_t)   :: table
    type(fortnum_status_t) :: status
    character(len=4096)    :: path
    integer                :: arglen, argstat

    call get_command_argument(1, path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") &
            "usage: test_fortnum_special_dawson_oracle <dawson.csv>"
        stop 1
    end if

    call oracle_read(path(1:arglen), table, status)
    if (.not. status_ok(status)) then
        write (error_unit, "(a)") "oracle_read failed: "//trim(status%msg)
        stop 1
    end if

    if (size(table%x) < 1) then
        write (error_unit, "(a)") "empty oracle table"
        stop 1
    end if

    call oracle_check(table, dawson_wrapper, 1.0e-14_dp, 1.0e-13_dp, status)
    if (.not. status_ok(status)) then
        write (error_unit, "(a)") "dawson failed oracle: "//trim(status%msg)
        stop 1
    end if

    write (*, "(a,i0,a)") "dawson oracle passed: ", size(table%x), " rows"
    stop 0

contains

    pure function dawson_wrapper(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = dawson(x)
    end function dawson_wrapper

end program test_fortnum_special_dawson_oracle
