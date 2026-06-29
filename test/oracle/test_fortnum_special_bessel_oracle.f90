program test_fortnum_special_bessel_oracle
    ! Oracle test for fortnum_special_bessel against scipy.special.iv / kn.
    ! CSV columns: index,func,n,x,expected
    !   func 0 -> I_n(x) via bessel_in
    !   func 1 -> K_n(x) via bessel_kn
    ! Tolerances:
    !   - values below 1e-280 on either side count as matching underflow
    !   - relative tolerance 2e-13 elsewhere (I_n n=200 near series/Miller
    !     transition reaches ~1.2e-13 vs scipy; the prior libneo oracle used 5e-13)
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special_bessel, only: bessel_in, bessel_kn

    implicit none

    ! 2e-13: I_n at n=200 near the series/Miller boundary reaches ~1.2e-13
    ! relative error vs scipy; the prior libneo oracle used 5e-13.
    real(dp), parameter :: rel_tol        = 2.0e-13_dp
    real(dp), parameter :: underflow_tol  = 1.0e-280_dp

    character(len=4096) :: csv_path
    integer             :: arglen, argstat, unit, ios, row_idx
    character(len=512)  :: line

    integer  :: func, n, nfail
    real(dp) :: x, expected, got, err, scale
    real(dp) :: worst_i, worst_k

    call get_command_argument(1, csv_path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") &
            "usage: test_fortnum_special_bessel_oracle <bessel.csv>"
        stop 1
    end if

    open (newunit=unit, file=csv_path(1:arglen), status="old", &
        action="read", iostat=ios)
    if (ios /= 0) then
        write (error_unit, "(a)") &
            "cannot open "//csv_path(1:arglen)
        stop 1
    end if

    nfail   = 0
    row_idx = 0
    worst_i = 0.0_dp
    worst_k = 0.0_dp

    do
        read (unit, "(a)", iostat=ios) line
        if (ios /= 0) exit
        ! skip comment / blank
        if (is_comment(line) .or. len_trim(line) == 0) cycle

        call parse_row(line, row_idx, func, n, x, expected, ios)
        if (ios /= 0) then
            write (error_unit, "(a,i0)") &
                "malformed row at index ", row_idx
            nfail = nfail + 1
            cycle
        end if

        select case (func)
        case (0) ! I_n
            got = bessel_in(n, x)
        case (1) ! K_n
            got = bessel_kn(n, x)
        case default
            write (error_unit, "(a,i0,a,i0)") &
                "unknown func ", func, " at row ", row_idx
            nfail = nfail + 1
            cycle
        end select

        ! Underflow both sides -> pass
        if (abs(got) < underflow_tol .and. abs(expected) < underflow_tol) cycle

        scale = max(abs(expected), underflow_tol)
        err   = abs(got - expected)/scale
        if (func == 0) worst_i = max(worst_i, err)
        if (func == 1) worst_k = max(worst_k, err)

        if (err > rel_tol) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,i0,a,i0,a,es13.6,a,es24.16,a,es24.16,a,es13.6)") &
                "FAIL row ", row_idx, " func=", func, " n=", n, " x=", x, &
                " expected=", expected, " got=", got, " relerr=", err
        end if
    end do
    close (unit)

    write (*, "(a,es13.6)") "worst I_n rel err: ", worst_i
    write (*, "(a,es13.6)") "worst K_n rel err: ", worst_k

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " failures"
        stop 1
    end if

    write (*, "(a)") "PASS"
    stop 0

contains

    pure logical function is_comment(line)
        character(*), intent(in) :: line
        integer :: p
        is_comment = .false.
        p = verify(line, " ")
        if (p > 0) is_comment = (line(p:p) == "#")
    end function is_comment

    subroutine parse_row(line, idx, func, n, x, expected, ios)
        character(*), intent(in)  :: line
        integer,      intent(out) :: idx, func, n, ios
        real(dp),     intent(out) :: x, expected

        character(len=len(line)) :: buf
        integer :: i

        buf = line
        do i = 1, len(buf)
            if (buf(i:i) == ",") buf(i:i) = " "
        end do
        read (buf, *, iostat=ios) idx, func, n, x, expected
    end subroutine parse_row

end program test_fortnum_special_bessel_oracle
