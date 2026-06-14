program test_fortnum_special_hyperg_oracle
    ! Oracle test for fortnum_special_hypergeometric_1f1 against mpmath.hyp1f1.
    !
    ! Reads test/oracle/data/hyperg_1f1.csv produced by mpmath.hyp1f1 (dps=50),
    ! cross-checked on the real axis against scipy.special.hyp1f1.
    !   columns: index,a_re,a_im,b_re,b_im,z_re,z_im,m_re,m_im
    !
    ! Tolerance: relative error 5e-12 on the complex value.  The Taylor and
    ! Kummer branches reach near machine precision; the large-|z| asymptotic
    ! branch is the worst case (optimal-truncation tail) and sets the bound.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_special_hypergeometric_1f1, only: hyperg_1f1
    implicit none

    real(dp), parameter :: rel_tol = 5.0e-12_dp

    character(len=4096) :: csv_path
    integer  :: arglen, argstat, unit, ios, nfail, nrow, idx
    character(len=1024) :: line, buf
    real(dp) :: a_re, a_im, b_re, b_im, z_re, z_im, m_re, m_im
    real(dp) :: err, scale, worst
    complex(dp) :: a, b, z, ref, got
    type(fortnum_status_t) :: status

    call get_command_argument(1, csv_path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") &
            "usage: test_fortnum_special_hyperg_oracle <hyperg_1f1.csv>"
        stop 1
    end if

    open (newunit=unit, file=csv_path(1:arglen), status="old", &
          action="read", iostat=ios)
    if (ios /= 0) then
        write (error_unit, "(a)") "cannot open " // csv_path(1:arglen)
        stop 1
    end if

    nfail = 0
    nrow  = 0
    worst = 0.0_dp

    do
        read (unit, "(a)", iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        buf = adjustl(line)
        if (buf(1:1) == "#") cycle

        call commas_to_spaces(line)
        read (line, *, iostat=ios) idx, a_re, a_im, b_re, b_im, &
            z_re, z_im, m_re, m_im
        if (ios /= 0) then
            write (error_unit, "(a,i0)") "parse error near row ", nrow
            nfail = nfail + 1
            cycle
        end if
        nrow = nrow + 1

        a   = cmplx(a_re, a_im, dp)
        b   = cmplx(b_re, b_im, dp)
        z   = cmplx(z_re, z_im, dp)
        ref = cmplx(m_re, m_im, dp)

        call hyperg_1f1(a, b, z, got, status)
        if (.not. status_ok(status)) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,a)") "FAIL row ", idx, &
                " status: ", trim(status%msg)
            cycle
        end if

        scale = max(abs(ref), 1.0_dp)
        err   = abs(got - ref) / scale
        worst = max(worst, err)
        if (err > rel_tol) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,es12.5,a,2es24.16,a,2es24.16)") &
                "FAIL row ", idx, " relerr=", err, &
                " ref=", real(ref, dp), aimag(ref), &
                " got=", real(got, dp), aimag(got)
        end if
    end do
    close (unit)

    if (nrow == 0) then
        write (error_unit, "(a)") "FAIL: no data rows read"
        stop 1
    end if

    write (*, "(a,es13.6)") "worst rel err: ", worst
    if (nfail > 0) then
        write (error_unit, "(i0,a,i0,a)") nfail, " failures in ", nrow, " rows"
        stop 1
    end if
    write (*, "(a,i0,a)") "PASS: ", nrow, " rows within tolerance"
    stop 0

contains

    subroutine commas_to_spaces(s)
        character(len=*), intent(inout) :: s
        integer :: i
        do i = 1, len(s)
            if (s(i:i) == ",") s(i:i) = " "
        end do
    end subroutine commas_to_spaces

end program test_fortnum_special_hyperg_oracle
