program test_fortnum_levin_oracle
    ! Oracle test for fortnum_levin: accelerates the partial sums of each test
    ! series with levin_u_accel and compares the result against the closed-form
    ! limit in the reference table (confirmed by mpmath's u-variant Levin).
    !
    ! CSV layout (levin_u.csv):
    !   H,case_id,n,expected,tol        -- one header row per case
    !   T,case_id,i,term                -- n term rows per case (i = 0..n-1)
    ! The per-case tol is the honest float64 level the transform reaches on that
    ! series: ~1e-16 for the MEPHIT/KiLCA geometric Kummer ratio (cases 3-5),
    ! float64-cancellation limited (~1e-9) for the alternating cases.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_levin, only: levin_u_accel
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    character(len=4096) :: path
    integer             :: arglen, argstat, unit, ios

    integer, parameter :: max_cases = 32
    integer, parameter :: max_terms = 200

    character(len=1)   :: tag
    character(len=512) :: line

    integer  :: case_id(max_cases), case_n(max_cases), ncases
    real(dp) :: case_expected(max_cases), case_tol(max_cases)
    real(dp) :: terms(max_terms, max_cases)

    integer  :: c, idx, n, nfail, ci
    real(dp) :: got, abserr, err
    type(fortnum_status_t) :: status

    call get_command_argument(1, path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") &
            "usage: test_fortnum_levin_oracle <levin_u.csv>"
        stop 1
    end if

    open (newunit=unit, file=path(1:arglen), status="old", action="read", &
          iostat=ios)
    if (ios /= 0) then
        write (error_unit, "(a)") "cannot open: "//path(1:arglen)
        stop 1
    end if

    ncases = 0
    do
        read (unit, "(a)", iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        if (line(verify(line, " "):verify(line, " ")) == "#") cycle
        tag = line(1:1)
        call replace_commas(line)
        if (tag == "H") then
            ncases = ncases + 1
            if (ncases > max_cases) then
                write (error_unit, "(a)") "oracle: too many cases"
                stop 1
            end if
            read (line(2:), *, iostat=ios) case_id(ncases), case_n(ncases), &
                case_expected(ncases), case_tol(ncases)
            if (ios /= 0) then
                write (error_unit, "(a)") "oracle: header parse error"
                stop 1
            end if
            if (case_n(ncases) > max_terms) then
                write (error_unit, "(a)") "oracle: too many terms"
                stop 1
            end if
        else if (tag == "T") then
            read (line(2:), *, iostat=ios) ci, idx, got
            if (ios /= 0) then
                write (error_unit, "(a)") "oracle: term parse error"
                stop 1
            end if
            c = case_index(case_id, ncases, ci)
            if (c == 0) then
                write (error_unit, "(a,i0)") "oracle: term for unknown case ", ci
                stop 1
            end if
            terms(idx + 1, c) = got
        else
            write (error_unit, "(a)") "oracle: unknown row tag: "//tag
            stop 1
        end if
    end do
    close (unit)

    nfail = 0
    do c = 1, ncases
        n = case_n(c)
        call levin_u_accel(terms(1:n, c), n, got, abserr, status)
        if (.not. status_ok(status)) then
            write (error_unit, "(a,i0,a)") "oracle FAIL case ", case_id(c), &
                ": status "//trim(status%msg)
            nfail = nfail + 1
            cycle
        end if
        err = abs(got - case_expected(c))
        if (err > case_tol(c)) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,es24.16,a,es24.16,a,es12.4,a,es12.4)") &
                "oracle FAIL case ", case_id(c), ": got=", got, &
                " want=", case_expected(c), " abserr=", err, " tol=", case_tol(c)
        end if
    end do

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " oracle case(s) failed"
        stop 1
    end if

    write (*, "(a,i0,a)") "oracle passed: ", ncases, " Levin-u cases verified"
    stop 0

contains

    pure integer function case_index(ids, ncases, want) result(c)
        integer, intent(in) :: ids(:), ncases, want
        integer :: k
        c = 0
        do k = 1, ncases
            if (ids(k) == want) then
                c = k
                return
            end if
        end do
    end function case_index

    subroutine replace_commas(buf)
        character(*), intent(inout) :: buf
        integer :: k
        do k = 1, len(buf)
            if (buf(k:k) == ",") buf(k:k) = " "
        end do
    end subroutine replace_commas

end program test_fortnum_levin_oracle
