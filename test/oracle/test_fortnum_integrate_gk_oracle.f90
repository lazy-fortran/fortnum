program test_fortnum_integrate_gk_oracle
    ! Verify integrate_gk against scipy.integrate.quad / analytic reference
    ! values loaded from test/oracle/data/integrate_gk.csv.
    !
    ! CSV columns: case_id, a, b, expected_integral
    ! (no single-argument x column; integrand is a functional argument, so the
    ! generic fortnum_oracle helper does not apply here).
    !
    ! Tolerance: atol=1e-10, rtol=1e-10. The GK21 rule with adaptive
    ! bisection (limit=200) achieves ~1e-13 on smooth integrands; the chosen
    ! bound leaves room for the sqrt integrand (slower convergence near 0).

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_integrate_gk, only: integrate_gk
    implicit none

    integer,  parameter :: NCASES = 6
    real(dp), parameter :: ATOL   = 1.0e-10_dp
    real(dp), parameter :: RTOL   = 1.0e-10_dp
    real(dp), parameter :: PI     = 3.14159265358979323846264338327950288_dp

    character(len=4096) :: csv_path
    integer             :: arglen, argstat, unit, ios, nfail, i
    character(len=512)  :: line
    integer             :: cid
    real(dp)            :: a, b, ref
    real(dp)            :: result, abserr, err, tol
    integer             :: ierr

    ! Integrand dispatch via a module-level selector; avoids heap closures.
    integer :: active_case

    call get_command_argument(1, csv_path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") &
            "usage: test_fortnum_integrate_gk_oracle <integrate_gk.csv>"
        stop 1
    end if

    open (newunit=unit, file=csv_path(1:arglen), status="old", action="read", &
        iostat=ios)
    if (ios /= 0) then
        write (error_unit, "(a)") &
            "cannot open: "//csv_path(1:arglen)
        stop 1
    end if

    nfail = 0
    do
        read (unit, "(a)", iostat=ios) line
        if (ios /= 0) exit
        ! skip header/comment lines
        if (is_comment(line) .or. len_trim(line) == 0) cycle
        call parse_row(line, cid, a, b, ref, ios)
        if (ios /= 0) then
            write (error_unit, "(a)") "malformed row: "//trim(line)
            nfail = nfail + 1
            cycle
        end if
        active_case = cid
        call integrate_gk(dispatch, a, b, 0.0_dp, 1.0e-12_dp, result, abserr, &
            ierr, key=21, limit=200)
        if (ierr /= 0) then
            write (error_unit, "(a,i0,a,i0)") &
                "case ", cid, ": ierr = ", ierr
            nfail = nfail + 1
            cycle
        end if
        err = abs(result - ref)
        tol = ATOL + RTOL*abs(ref)
        if (.not. (err <= tol)) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,4(a,es20.12e3))") &
                "FAIL case ", cid, &
                "  ref=", ref, "  got=", result, &
                "  abserr=", abserr, "  diff=", err
        end if
    end do
    close (unit)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " case(s) failed"
        stop 1
    end if
    write (*, "(a,i0,a)") &
        "integrate_gk oracle: ", NCASES, " cases passed"
    stop 0

contains

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

    ! Integrand dispatcher; active_case is set before each integrate_gk call.
    ! Not pure: reads module variable active_case.
    function dispatch(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx
        select case (active_case)
        case (0) ! exp(x) on [0,1]; analytic = e-1
            fx = exp(x)
        case (1) ! 3x^2+2x+1 on [0,2]; analytic = 14
            fx = 3.0_dp*x**2 + 2.0_dp*x + 1.0_dp
        case (2) ! sin(x) on [0,pi]; analytic = 2
            fx = sin(x)
        case (3) ! cos(x) on [0,pi/2]; analytic = 1
            fx = cos(x)
        case (4) ! sqrt(x) on [0,1]; analytic = 2/3
            fx = sqrt(x)
        case (5) ! exp(-x^2) on [0,1]; ref from scipy
            fx = exp(-x*x)
        case default
            fx = 0.0_dp
        end select
    end function dispatch

end program test_fortnum_integrate_gk_oracle
