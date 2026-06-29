program test_fortnum_integrate_oracle
    ! Verify integrate_qag / integrate_qags against scipy.integrate.quad
    ! (QUADPACK) reference values in test/oracle/data/integrate_qag_smooth.csv
    ! and integrate_qags_singular.csv.
    !
    ! CSV columns: case_id, a, b, expected. The integrand is a functional
    ! argument selected by case_id, so the generic fortnum_oracle helper does
    ! not apply; this test parses the CSV directly.
    !
    ! Match within the requested relative tolerance: epsrel=1e-9 on smooth
    ! cases, 1e-7 on the singular/peaked QAGS cases (where scipy itself runs
    ! near its roundoff floor on the sharpest peak).

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t
    use fortnum_integrate, only: integrate_qag, integrate_qags, &
                                 integrate_workspace_t, integrate_epstab_t, &
                                 integrate_result_t
    implicit none

    character(len=4096) :: qag_path, qags_path
    integer :: arglen, argstat, nfail

    ! Integrand selectors read by dispatch through host association.
    integer :: active_case
    logical :: active_qags

    call get_command_argument(1, qag_path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") "usage: test_fortnum_integrate_oracle "// &
            "<integrate_qag_smooth.csv> <integrate_qags_singular.csv>"
        stop 1
    end if
    qag_path = qag_path(1:arglen)
    call get_command_argument(2, qags_path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") "usage: test_fortnum_integrate_oracle "// &
            "<integrate_qag_smooth.csv> <integrate_qags_singular.csv>"
        stop 1
    end if
    qags_path = qags_path(1:arglen)

    nfail = 0
    call run_file(trim(qag_path), .false., 1.0e-9_dp, nfail)
    call run_file(trim(qags_path), .true., 1.0e-7_dp, nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " case(s) failed"
        stop 1
    end if
    write (*, "(a)") "integrate oracle: all cases passed"
    stop 0

contains

    subroutine run_file(path, use_qags, rtol, nfail)
        character(*), intent(in)    :: path
        logical,      intent(in)    :: use_qags
        real(dp),     intent(in)    :: rtol
        integer,      intent(inout) :: nfail
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        character(len=512) :: line
        integer  :: unit, ios, cid
        real(dp) :: a, b, ref, err, tol
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
            call parse_row(line, cid, a, b, ref, ios)
            if (ios /= 0) then
                write (error_unit, "(a)") "malformed row: "//trim(line)
                nfail = nfail + 1
                cycle
            end if
            active_case = cid
            active_qags = use_qags
            if (use_qags) then
                call integrate_qags(dispatch, a, b, 0.0_dp, rtol, ws, eps, &
                                    res, st, limit=500)
            else
                call integrate_qag(dispatch, a, b, 0.0_dp, rtol, ws, res, st, &
                                   key=21, limit=500)
            end if
            err = abs(res%value - ref)
            ! Absolute floor mirrors scipy's own attainable accuracy.
            tol = max(rtol*abs(ref), 1.0e-9_dp)
            if (.not. (err <= tol)) then
                nfail = nfail + 1
                write (error_unit, "(a,l1,a,i0,3(a,es20.12e3),a,i0)") &
                    "FAIL qags=", use_qags, " case ", cid, &
                    "  ref=", ref, "  got=", res%value, "  diff=", err, &
                    "  status=", st%code
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

    ! Integrand dispatcher. active_case/active_qags select the integrand to
    ! match the python generator's case ids; both files reuse case ids 0..n.
    function dispatch(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        if (active_qags) then
            select case (active_case)
            case (0);  fx = 1.0_dp/sqrt(x)
            case (1);  fx = log(x)
            case (2);  fx = 1.0_dp/(1.0_dp + ((x - 0.5_dp)/1.0e-2_dp)**2)
            case (3);  fx = 1.0_dp/(1.0_dp + ((x - 0.3_dp)/1.0e-3_dp)**2)
            case (4);  fx = x**(-0.5_dp)*exp(-x)
            case default; fx = 0.0_dp
            end select
        else
            select case (active_case)
            case (0);  fx = 3.0_dp*x**2 + 2.0_dp*x + 1.0_dp
            case (1);  fx = exp(x)
            case (2);  fx = cos(x)
            case (3);  fx = sin(10.0_dp*x) + 2.0_dp
            case (4);  fx = 1.0_dp/(1.0_dp + x*x)
            case (5);  fx = exp(-x*x)
            case default; fx = 0.0_dp
            end select
        end if
    end function dispatch

end program test_fortnum_integrate_oracle

