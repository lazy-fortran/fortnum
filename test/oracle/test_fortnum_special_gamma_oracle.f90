program test_fortnum_special_gamma_oracle
    ! Oracle test for fortnum_special_gamma.
    !
    ! Reads test/oracle/data/lower_incomplete_gamma.csv, which was produced by
    !   scipy.special.gammainc(a,x) * scipy.special.gamma(a)  (gamma_lower)
    !   scipy.special.gammainc(a,x)                            (gamma_reg_p)
    !
    ! Tolerances: atol=1e-13, rtol=1e-13 — well within float64 agreement
    ! between the Lentz/series algorithms and scipy's Cephes implementation.
    !
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special_gamma, only: gamma_lower, gamma_reg_p
    implicit none

    character(len=*), parameter :: csv_default = &
        "test/oracle/data/lower_incomplete_gamma.csv"

    real(dp), parameter :: atol = 1.0e-13_dp
    real(dp), parameter :: rtol = 1.0e-13_dp

    character(len=512) :: csv_path
    integer  :: narg, unit, ios, nfail, nrow
    character(len=1024) :: line, buf
    integer  :: idx
    real(dp) :: a, x, ref_gl, ref_rp, got_gl, got_rp, err, tol

    narg = command_argument_count()
    if (narg >= 1) then
        call get_command_argument(1, csv_path)
    else
        csv_path = csv_default
    end if

    open(newunit=unit, file=trim(csv_path), status="old", action="read", &
         iostat=ios)
    if (ios /= 0) then
        write(error_unit, "(a)") "FAIL: cannot open " // trim(csv_path)
        stop 1
    end if

    nfail = 0
    nrow  = 0

    do
        read(unit, "(a)", iostat=ios) line
        if (ios /= 0) exit
        ! skip comments and blank lines
        if (len_trim(line) == 0) cycle
        buf = adjustl(line)
        if (buf(1:1) == "#") cycle

        ! parse: index,a,x,gamma_lower,gamma_reg_p
        call replace_commas(line)
        read(line, *, iostat=ios) idx, a, x, ref_gl, ref_rp
        if (ios /= 0) then
            write(error_unit, "(a,i0)") "parse error at row ", nrow
            nfail = nfail + 1
            cycle
        end if
        nrow = nrow + 1

        got_gl = gamma_lower(a, x)
        got_rp = gamma_reg_p(a, x)

        ! check gamma_lower
        err = abs(got_gl - ref_gl)
        tol = atol + rtol * abs(ref_gl)
        if (.not. (err <= tol)) then
            nfail = nfail + 1
            write(error_unit, "(a,i0,a,2(es14.7,a),es10.3)") &
                "FAIL gamma_lower row ", idx, &
                ": a=", a, " x=", x, " err=", err
        end if

        ! check gamma_reg_p
        err = abs(got_rp - ref_rp)
        tol = atol + rtol * abs(ref_rp)
        if (.not. (err <= tol)) then
            nfail = nfail + 1
            write(error_unit, "(a,i0,a,2(es14.7,a),es10.3)") &
                "FAIL gamma_reg_p row ", idx, &
                ": a=", a, " x=", x, " err=", err
        end if
    end do

    close(unit)

    if (nrow == 0) then
        write(error_unit, "(a)") "FAIL: no data rows read from " // trim(csv_path)
        stop 1
    end if

    if (nfail > 0) then
        write(error_unit, "(i0,a,i0,a)") nfail, " failures in ", nrow, " rows"
        stop 1
    end if

    write(*, "(a,i0,a)") "PASS: ", nrow, " rows within tolerance"

contains

    subroutine replace_commas(s)
        character(len=*), intent(inout) :: s
        integer :: i
        do i = 1, len(s)
            if (s(i:i) == ",") s(i:i) = " "
        end do
    end subroutine replace_commas

end program test_fortnum_special_gamma_oracle
