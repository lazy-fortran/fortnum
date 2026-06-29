program test_fortnum_roots_oracle
    ! Oracle tests for fortnum_roots: bisection and Newton agree with
    ! high-precision reference values (mpmath / analytic) in roots.csv.
    !
    ! Usage:
    !   test_fortnum_roots_oracle <roots.csv>
    !
    ! CSV columns: index,a,b,root  (comment lines start with #)
    ! Each case is run through both root_bisect and root_newton.
    ! Tolerance: 2 ulp at double precision ~ 4.5e-16.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_roots, only: root_bisect, root_newton
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    character(len=4096) :: csv_path
    integer :: arglen, argstat, nfail, total

    call get_command_argument(1, csv_path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") "usage: test_fortnum_roots_oracle <roots.csv>"
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
        integer, intent(inout)       :: nfail, total

        integer, parameter :: MAX_CASES = 64
        ! Bisection converges to within 4*epsilon of the bracket; the
        ! returned midpoint can sit up to ~3 ulp from the true root.
        ! Newton is quadratically convergent and typically delivers the
        ! last bit.  Both are checked at the same tolerance; 2e-15 covers
        ! bisection's inherent rounding without masking algorithmic errors.
        real(dp), parameter :: ATOL = 2.0e-15_dp ! ~9 ulp at dp

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
            call check_bisect(idx_arr(i), a_arr(i), b_arr(i), ref_arr(i), &
                ATOL, nfail, total)
            call check_newton(idx_arr(i), a_arr(i), b_arr(i), ref_arr(i), &
                ATOL, nfail, total)
        end do
    end subroutine run_all

    subroutine check_bisect(idx, a, b, ref, atol, nfail, total)
        integer,  intent(in)    :: idx
        real(dp), intent(in)    :: a, b, ref, atol
        integer,  intent(inout) :: nfail, total

        type(fortnum_status_t) :: s
        real(dp) :: x, err

        select case (idx)
        case (0); call root_bisect(fn0, a, b, x, s)
        case (1); call root_bisect(fn1, a, b, x, s)
        case (2); call root_bisect(fn2, a, b, x, s)
        case (3); call root_bisect(fn3, a, b, x, s)
        case (4); call root_bisect(fn4, a, b, x, s)
        case default
            write (error_unit, "(a,i0)") "oracle: unknown case ", idx
            nfail = nfail + 1
            return
        end select

        total = total + 1
        if (s%code /= FORTNUM_OK) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,i0,a,a)") &
                "FAIL bisect case ", idx, " status=", s%code, " ", trim(s%msg)
            return
        end if
        err = abs(x - ref)
        if (err > atol) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,es12.4,a,es24.16,a,es24.16)") &
                "FAIL bisect case ", idx, " abserr=", err, &
                " got=", x, " ref=", ref
        end if
    end subroutine check_bisect

    subroutine check_newton(idx, a, b, ref, atol, nfail, total)
        integer,  intent(in)    :: idx
        real(dp), intent(in)    :: a, b, ref, atol
        integer,  intent(inout) :: nfail, total

        type(fortnum_status_t) :: s
        real(dp) :: x, xmid, err

        xmid = 0.5_dp * (a + b)
        select case (idx)
        case (0); call root_newton(fdf0, a, b, xmid, x, s)
        case (1); call root_newton(fdf1, a, b, xmid, x, s)
        case (2); call root_newton(fdf2, a, b, xmid, x, s)
        case (3); call root_newton(fdf3, a, b, xmid, x, s)
        case (4); call root_newton(fdf4, a, b, xmid, x, s)
        case default
            write (error_unit, "(a,i0)") "oracle: unknown case ", idx
            nfail = nfail + 1
            return
        end select

        total = total + 1
        if (s%code /= FORTNUM_OK) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,i0,a,a)") &
                "FAIL newton case ", idx, " status=", s%code, " ", trim(s%msg)
            return
        end if
        err = abs(x - ref)
        if (err > atol) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,es12.4,a,es24.16,a,es24.16)") &
                "FAIL newton case ", idx, " abserr=", err, &
                " got=", x, " ref=", ref
        end if
    end subroutine check_newton

    ! ------------------------------------------------------------------ f/fdf

    ! case 0: x^2 - 2,  root = sqrt(2)
    pure function fn0(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x * x - 2.0_dp
    end function fn0

    pure subroutine fdf0(x, fx, dfx)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: fx, dfx
        fx  = x * x - 2.0_dp
        dfx = 2.0_dp * x
    end subroutine fdf0

    ! case 1: cos(x) - x,  Dottie number
    pure function fn1(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = cos(x) - x
    end function fn1

    pure subroutine fdf1(x, fx, dfx)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: fx, dfx
        fx  = cos(x) - x
        dfx = -sin(x) - 1.0_dp
    end subroutine fdf1

    ! case 2: x^3 - x - 2
    pure function fn2(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x*x*x - x - 2.0_dp
    end function fn2

    pure subroutine fdf2(x, fx, dfx)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: fx, dfx
        fx  = x*x*x - x - 2.0_dp
        dfx = 3.0_dp * x*x - 1.0_dp
    end subroutine fdf2

    ! case 3: exp(x) - 3,  root = ln(3)
    pure function fn3(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = exp(x) - 3.0_dp
    end function fn3

    pure subroutine fdf3(x, fx, dfx)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: fx, dfx
        fx  = exp(x) - 3.0_dp
        dfx = exp(x)
    end subroutine fdf3

    ! case 4: x^3 - 7x + 6,  root = 1 (exact)
    pure function fn4(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x*x*x - 7.0_dp*x + 6.0_dp
    end function fn4

    pure subroutine fdf4(x, fx, dfx)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: fx, dfx
        fx  = x*x*x - 7.0_dp*x + 6.0_dp
        dfx = 3.0_dp * x*x - 7.0_dp
    end subroutine fdf4

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

end program test_fortnum_roots_oracle
