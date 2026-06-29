program test_fortnum_integrate_singular_oracle
    ! Verify integrate_qagp / integrate_qagiu against scipy.integrate.quad
    ! (QUADPACK dqagp/dqagi) reference values in
    ! test/oracle/data/integrate_qagp_singular.csv and
    ! integrate_qagiu_infinite.csv.
    !
    ! QAGP CSV columns: case_id, a, b, p1, p2, expected (p2 = NaN means one
    ! break point). QAGIU CSV columns: case_id, bound, inf, expected. The
    ! integrand is selected by case_id through host association, so this test
    ! parses the CSVs directly rather than via the generic oracle helper.
    !
    ! Request epsrel=1e-9 and match within 1e-6 relative, 1e-7 absolute floor.
    ! The interior algebraic singularities (|x-c|^(-1/2)) converge through the
    ! epsilon table, whose error estimate is approximate, so the attainable
    ! accuracy on those cases is ~1e-8; 1e-6 is the comfortable assertion band.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t
    use fortnum_integrate, only: integrate_qagp, integrate_qagiu, &
        integrate_workspace_t, integrate_epstab_t, &
        integrate_result_t
    implicit none

    character(len=4096) :: qagp_path, qagiu_path
    integer :: arglen, argstat, nfail

    integer :: active_qagp_case, active_qagiu_case

    call get_command_argument(1, qagp_path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") "usage: "// &
            "test_fortnum_integrate_singular_oracle "// &
            "<integrate_qagp_singular.csv> <integrate_qagiu_infinite.csv>"
        stop 1
    end if
    qagp_path = qagp_path(1:arglen)
    call get_command_argument(2, qagiu_path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") "usage: "// &
            "test_fortnum_integrate_singular_oracle "// &
            "<integrate_qagp_singular.csv> <integrate_qagiu_infinite.csv>"
        stop 1
    end if
    qagiu_path = qagiu_path(1:arglen)

    nfail = 0
    call run_qagp(trim(qagp_path), nfail)
    call run_qagiu(trim(qagiu_path), nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " case(s) failed"
        stop 1
    end if
    write (*, "(a)") "integrate singular oracle: all cases passed"
    stop 0

contains

    subroutine run_qagp(path, nfail)
        character(*), intent(in)    :: path
        integer,      intent(inout) :: nfail
        real(dp), parameter :: req = 1.0e-9_dp, atol = 1.0e-6_dp
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        character(len=512) :: line
        integer  :: unit, ios, cid, npts
        real(dp) :: a, b, p1, p2, ref, err, tol
        real(dp) :: pts(2)
        open (newunit=unit, file=path, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            write (error_unit, "(a)") "cannot open: "//path
            nfail = nfail + 1
            return
        end if
        do
            read (unit, "(a)", iostat=ios) line
            if (ios /= 0) exit
            if (is_comment(line) .or. len_trim(line) == 0) cycle
            call parse_qagp(line, cid, a, b, p1, p2, ref, ios)
            if (ios /= 0) then
                write (error_unit, "(a)") "malformed row: "//trim(line)
                nfail = nfail + 1
                cycle
            end if
            ! p2 != p2 marks the NaN sentinel for a one-break-point case.
            if (p2 == p2) then
                npts = 2
                pts(1) = p1
                pts(2) = p2
            else
                npts = 1
                pts(1) = p1
            end if
            active_qagp_case = cid
            call integrate_qagp(dispatch_qagp, a, b, pts(1:npts), 0.0_dp, &
                req, ws, eps, res, st, limit=500)
            err = abs(res%value - ref)
            tol = max(atol*abs(ref), 1.0e-7_dp)
            if (.not. (err <= tol)) then
                nfail = nfail + 1
                write (error_unit, "(a,i0,3(a,es20.12e3),a,i0)") &
                    "FAIL qagp case ", cid, "  ref=", ref, "  got=", &
                    res%value, "  diff=", err, "  status=", st%code
            end if
        end do
        close (unit)
    end subroutine run_qagp

    subroutine run_qagiu(path, nfail)
        character(*), intent(in)    :: path
        integer,      intent(inout) :: nfail
        real(dp), parameter :: req = 1.0e-9_dp, atol = 1.0e-6_dp
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        character(len=512) :: line
        integer  :: unit, ios, cid, inf
        real(dp) :: bound, ref, err, tol
        open (newunit=unit, file=path, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            write (error_unit, "(a)") "cannot open: "//path
            nfail = nfail + 1
            return
        end if
        do
            read (unit, "(a)", iostat=ios) line
            if (ios /= 0) exit
            if (is_comment(line) .or. len_trim(line) == 0) cycle
            call parse_qagiu(line, cid, bound, inf, ref, ios)
            if (ios /= 0) then
                write (error_unit, "(a)") "malformed row: "//trim(line)
                nfail = nfail + 1
                cycle
            end if
            active_qagiu_case = cid
            call integrate_qagiu(dispatch_qagiu, bound, inf, 0.0_dp, req, &
                ws, eps, res, st, limit=500)
            err = abs(res%value - ref)
            tol = max(atol*abs(ref), 1.0e-7_dp)
            if (.not. (err <= tol)) then
                nfail = nfail + 1
                write (error_unit, "(a,i0,3(a,es20.12e3),a,i0)") &
                    "FAIL qagiu case ", cid, "  ref=", ref, "  got=", &
                    res%value, "  diff=", err, "  status=", st%code
            end if
        end do
        close (unit)
    end subroutine run_qagiu

    pure logical function is_comment(ln)
        character(*), intent(in) :: ln
        integer :: p
        is_comment = .false.
        p = verify(ln, " ")
        if (p > 0) is_comment = (ln(p:p) == "#")
    end function is_comment

    subroutine parse_qagp(ln, cid, av, bv, p1v, p2v, refv, stat)
        character(*), intent(in)  :: ln
        integer,      intent(out) :: cid, stat
        real(dp),     intent(out) :: av, bv, p1v, p2v, refv
        character(len=len(ln)) :: buf
        integer :: k
        buf = ln
        do k = 1, len(buf)
            if (buf(k:k) == ",") buf(k:k) = " "
        end do
        read (buf, *, iostat=stat) cid, av, bv, p1v, p2v, refv
    end subroutine parse_qagp

    subroutine parse_qagiu(ln, cid, boundv, infv, refv, stat)
        character(*), intent(in)  :: ln
        integer,      intent(out) :: cid, infv, stat
        real(dp),     intent(out) :: boundv, refv
        character(len=len(ln)) :: buf
        integer :: k
        buf = ln
        do k = 1, len(buf)
            if (buf(k:k) == ",") buf(k:k) = " "
        end do
        read (buf, *, iostat=stat) cid, boundv, infv, refv
    end subroutine parse_qagiu

    ! QAGP integrands. Case ids match the python generator.
    function dispatch_qagp(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        select case (active_qagp_case)
        case (0);  fx = abs(x - 0.5_dp)**(-0.5_dp)
        case (1);  fx = abs(x - 1.0_dp/3.0_dp)**(-0.5_dp)
        case (2);  fx = abs(x - 0.3_dp) + abs(x - 0.7_dp)
        case (3);  fx = log(abs(x - 0.5_dp))
        case default; fx = 0.0_dp
        end select
    end function dispatch_qagp

    ! QAGIU integrands. Case ids match the python generator.
    function dispatch_qagiu(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        select case (active_qagiu_case)
        case (0);  fx = exp(-x)
        case (1);  fx = 1.0_dp/(1.0_dp + x*x)
        case (2);  fx = exp(-x)
        case (3);  fx = exp(x)
        case (4);  fx = exp(-x*x)
        case (5);  fx = 1.0_dp/(1.0_dp + x*x)
        case default; fx = 0.0_dp
        end select
    end function dispatch_qagiu

end program test_fortnum_integrate_singular_oracle
