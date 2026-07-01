program test_fortnum_cquad_oracle
    ! Verify integrate_cquad against scipy.integrate.quad (QUADPACK) reference
    ! values in test/oracle/data/cquad.csv. Both integrators converge to the
    ! true integral, so agreement to full double precision is expected.
    !
    ! CSV columns: case_id, a, b, expected. The integrand is a functional
    ! argument selected by case_id, parsed directly here.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_cquad, only: integrate_cquad
    implicit none

    character(len=4096) :: path
    integer :: arglen, argstat, nfail
    integer :: active_case

    call get_command_argument(1, path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") "usage: test_fortnum_cquad_oracle <cquad.csv>"
        stop 1
    end if
    path = path(1:arglen)

    nfail = 0
    call run_file(trim(path), 1.0e-11_dp, nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " case(s) failed"
        stop 1
    end if
    write (*, "(a)") "cquad oracle: all cases passed"
    stop 0

contains

    subroutine run_file(fpath, rtol, nfail)
        character(*), intent(in)    :: fpath
        real(dp),     intent(in)    :: rtol
        integer,      intent(inout) :: nfail
        type(fortnum_status_t) :: st
        character(len=512) :: line
        integer  :: unit, ios, cid
        real(dp) :: a, b, ref, val, err, tol

        open (newunit=unit, file=fpath, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            write (error_unit, "(a)") "cannot open: "//fpath
            nfail = nfail + 1
            return
        end if
        do
            read (unit, "(a)", iostat=ios) line
            if (ios /= 0) exit
            if (is_comment(line) .or. len_trim(line) == 0) cycle
            call parse_row(line, cid, a, b, ref, ios)
            if (ios /= 0) then
                write (error_unit, "(a)") "malformed row: "//trim(line)
                nfail = nfail + 1
                cycle
            end if
            active_case = cid
            call integrate_cquad(dispatch, a, b, val, st, epsabs=0.0_dp, &
                epsrel=1.0e-13_dp)
            err = abs(val - ref)
            tol = max(rtol*abs(ref), 1.0e-11_dp)
            if (st%code /= FORTNUM_OK .or. .not. (err <= tol)) then
                nfail = nfail + 1
                write (error_unit, "(a,i0,3(a,es20.12e3),a,i0)") &
                    "FAIL case ", cid, "  ref=", ref, "  got=", val, &
                    "  diff=", err, "  status=", st%code
            end if
        end do
        close (unit)
    end subroutine run_file

    pure logical function is_comment(ln)
        character(*), intent(in) :: ln
        integer :: p
        is_comment = .false.
        p = verify(ln, " ")
        if (p > 0) is_comment = (ln(p:p) == "#")
    end function is_comment

    subroutine parse_row(ln, case_id, av, bv, refv, stat)
        character(*), intent(in)  :: ln
        integer,      intent(out) :: case_id, stat
        real(dp),     intent(out) :: av, bv, refv
        character(len=len(ln)) :: buf
        integer :: k
        buf = ln
        do k = 1, len(buf)
            if (buf(k:k) == ",") buf(k:k) = " "
        end do
        read (buf, *, iostat=stat) case_id, av, bv, refv
    end subroutine parse_row

    ! Integrand dispatcher; case ids match gen_cquad_tables.py.
    function dispatch(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        select case (active_case)
        case (0);  fx = 3.0_dp*x**2 + 2.0_dp*x + 1.0_dp
        case (1);  fx = exp(x)
        case (2);  fx = 1.0_dp/(1.0_dp + x*x)
        case (3);  fx = exp(-x*x)
        case (4);  fx = sqrt(x)
        case (5);  fx = x*x*exp(-x*x)
        case (6);  fx = x**4*exp(-x*x)
        case (7);  fx = x*exp(-x)
        case (8);  fx = cos(5.0_dp*x)*exp(-x)
        case default; fx = 0.0_dp
        end select
    end function dispatch

end program test_fortnum_cquad_oracle
