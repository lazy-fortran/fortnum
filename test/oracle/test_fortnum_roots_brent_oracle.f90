program test_fortnum_roots_brent_oracle
    ! Oracle tests for root_brent: results must agree with high-precision
    ! mpmath references in roots_brent.csv within 2 ulp at dp.
    !
    ! Usage:
    !   test_fortnum_roots_brent_oracle <roots_brent.csv>
    !
    ! CSV columns: index,a,b,root  (comment lines start with #)
    ! Each case dispatches to a numbered function; Brent is run on [a, b].
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_roots, only: root_brent
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    character(len=4096) :: csv_path
    integer :: arglen, argstat, nfail, total

    call get_command_argument(1, csv_path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") &
            "usage: test_fortnum_roots_brent_oracle <roots_brent.csv>"
        stop 1
    end if

    nfail = 0
    total = 0

    call run_all(trim(csv_path), nfail, total)

    if (nfail > 0) then
        write (error_unit, "(i0,a,i0,a)") nfail, " oracle case(s) failed out of ", &
            total, " checked"
        stop 1
    end if
    write (*, "(a,i0,a)") "oracle passed: ", total, " cases verified"
    stop 0

contains

    subroutine run_all(path, nfail, total)
        character(len=*), intent(in) :: path
        integer,          intent(inout) :: nfail, total

        integer, parameter :: MAX_CASES = 64
        ! Brent converges to machine epsilon; allow ~9 ulp for rounding at
        ! the bracket endpoints and tolerance arithmetic.
        real(dp), parameter :: ATOL = 2.0e-15_dp

        integer  :: unit, ios, idx
        character(len=512) :: line, buf
        real(dp) :: a_arr(MAX_CASES), b_arr(MAX_CASES), ref_arr(MAX_CASES)
        integer  :: idx_arr(MAX_CASES), ncases, i

        open (newunit=unit, file=path, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            write (error_unit, "(a)") "FAIL: cannot open " // trim(path)
            nfail = nfail + 1
            return
        end if

        ncases = 0
        do
            read (unit, "(a)", iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (is_comment(line)) cycle
            ncases = ncases + 1
            if (ncases > MAX_CASES) then
                write (error_unit, "(a)") "oracle: too many rows in " // trim(path)
                stop 1
            end if
            buf = line
            call replace_commas(buf)
            read (buf, *, iostat=ios) idx, a_arr(ncases), b_arr(ncases), &
                ref_arr(ncases)
            idx_arr(ncases) = idx
            if (ios /= 0) then
                write (error_unit, "(a,i0)") "oracle: parse error at row ", ncases
                nfail = nfail + 1
                ncases = ncases - 1
            end if
        end do
        close (unit)

        do i = 1, ncases
            call check_brent(idx_arr(i), a_arr(i), b_arr(i), ref_arr(i), &
                ATOL, nfail, total)
        end do
    end subroutine run_all

    subroutine check_brent(idx, a, b, ref, atol, nfail, total)
        integer,  intent(in)    :: idx
        real(dp), intent(in)    :: a, b, ref, atol
        integer,  intent(inout) :: nfail, total

        type(fortnum_status_t) :: s
        real(dp) :: x, err

        select case (idx)
        case (0); call root_brent(fn0, a, b, x, s)
        case (1); call root_brent(fn1, a, b, x, s)
        case (2); call root_brent(fn2, a, b, x, s)
        case (3); call root_brent(fn3, a, b, x, s)
        case (4); call root_brent(fn4, a, b, x, s)
        case (5); call root_brent(fn5, a, b, x, s)
        case default
            write (error_unit, "(a,i0)") "oracle: unknown case ", idx
            nfail = nfail + 1
            return
        end select

        total = total + 1
        if (s%code /= FORTNUM_OK) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,i0,a,a)") &
                "FAIL brent case ", idx, " status=", s%code, " ", trim(s%msg)
            return
        end if
        err = abs(x - ref)
        if (err > atol) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,es12.4,a,es24.16,a,es24.16)") &
                "FAIL brent case ", idx, " abserr=", err, &
                " got=", x, " ref=", ref
        end if
    end subroutine check_brent

    ! ------------------------------------------------------------------ f defs

    ! case 0: x*exp(x) - 1 = 0,  root = W(1)
    pure function fn0(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x * exp(x) - 1.0_dp
    end function fn0

    ! case 1: atan(x) - pi/6 = 0,  root = tan(pi/6) = 1/sqrt(3)
    pure function fn1(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = atan(x) - 0.52359877559829887308_dp ! pi/6
    end function fn1

    ! case 2: x^5 - x - 1 = 0
    pure function fn2(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x*x*x*x*x - x - 1.0_dp
    end function fn2

    ! case 3: sin(x) - 4/5 = 0,  root = asin(4/5)
    pure function fn3(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = sin(x) - 0.8_dp
    end function fn3

    ! case 4: x^2 - 3 = 0,  root = sqrt(3)  (wide bracket [1, 3])
    pure function fn4(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x * x - 3.0_dp
    end function fn4

    ! case 5: x^3 + 4*x^2 - 10 = 0
    pure function fn5(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x*x*x + 4.0_dp*x*x - 10.0_dp
    end function fn5

    ! ------------------------------------------------------------------ utils

    pure logical function is_comment(line)
        character(*), intent(in) :: line
        integer :: p
        is_comment = .false.
        p = verify(line, " ")
        if (p > 0) is_comment = (line(p:p) == "#")
    end function is_comment

    pure subroutine replace_commas(buf)
        character(*), intent(inout) :: buf
        integer :: i
        do i = 1, len(buf)
            if (buf(i:i) == ",") buf(i:i) = " "
        end do
    end subroutine replace_commas

end program test_fortnum_roots_brent_oracle
