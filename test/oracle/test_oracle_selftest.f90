program test_oracle_selftest
    ! End-to-end exercise of the oracle harness: load the seeded gamma table
    ! and confirm gfortran's intrinsic gamma() reproduces it within tolerance,
    ! then confirm a deliberately wrong primal is rejected. The table path is
    ! the first command-line argument (CMake passes the source-tree location).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok, FORTNUM_CONVERGENCE_ERROR
    use fortnum_oracle, only: oracle_table_t, oracle_read, oracle_check
    implicit none

    type(oracle_table_t)   :: table
    type(fortnum_status_t) :: status
    character(len=4096)    :: path
    integer                :: arglen, argstat

    call get_command_argument(1, path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") "usage: test_oracle_selftest <gamma.csv>"
        stop 1
    end if

    call oracle_read(path(1:arglen), table, status)
    if (.not. status_ok(status)) then
        write (error_unit, "(a)") "read failed: "//trim(status%msg)
        stop 1
    end if

    ! Table must be non-empty and the seeded gamma table carries a derivative.
    if (size(table%x) < 1) stop 1
    if (.not. table%has_derivative) stop 1

    ! Correct primal: intrinsic gamma must match the reference. The reference
    ! came from an independent library, so a modest relative tolerance covers
    ! the last-bit disagreement between implementations.
    call oracle_check(table, intrinsic_gamma, 1.0e-12_dp, 1.0e-12_dp, status)
    if (.not. status_ok(status)) then
        write (error_unit, "(a)") "matching primal rejected: "//trim(status%msg)
        stop 1
    end if

    ! Wrong primal must be caught: shifting the argument breaks every row.
    call oracle_check(table, wrong_gamma, 1.0e-12_dp, 1.0e-12_dp, status)
    if (status_ok(status)) then
        write (error_unit, "(a)") "wrong primal was accepted"
        stop 1
    end if
    if (status%code /= FORTNUM_CONVERGENCE_ERROR) stop 1

    write (*, "(a,i0,a)") "oracle self-test passed: ", size(table%x), &
        " gamma rows verified"
    stop 0

contains

    ! Wrapper so the intrinsic can be passed where the abstract interface
    ! expects a user procedure.
    pure function intrinsic_gamma(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = gamma(x)
    end function intrinsic_gamma

    ! Deliberately incorrect: gamma(x+1) /= gamma(x) on the grid.
    pure function wrong_gamma(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = gamma(x + 1.0_dp)
    end function wrong_gamma

end program test_oracle_selftest
