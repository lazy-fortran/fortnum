program test_fortnum_bspline_oracle
    ! Oracle test for fortnum_bspline: basis functions B_{i,k}(x) and their
    ! derivatives agree with scipy.interpolate.BSpline (unit-coefficient splines)
    ! on the breakpoint sets and orders NEO-2 collop_bspline uses.
    !
    ! Usage:
    !   test_fortnum_bspline_oracle <bspline.csv>
    !
    ! File rows (comma-separated, '#' comments):
    !   brk,case,order,nbreak,b_0,...,b_{nbreak-1}
    !   val,case,order,nbreak,x,deriv,ncoef,v_0,...,v_{ncoef-1}
    ! A 'brk' row sets up the current workspace (breakpoints -> clamped knots).
    ! Each following 'val' row evaluates bspline_eval_deriv at x up to its deriv
    ! order and compares the case's basis column values against the reference.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, &
        bspline_set_knots, bspline_eval_deriv, bspline_eval_basis
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    ! scipy/de Boor and the Cox-de Boor recurrence here are the same rational
    ! arithmetic; agreement is at rounding level. 1e-12 covers the conditioning
    ! of the stretched NEO-2 knot set without masking an algorithmic error.
    real(dp), parameter :: ATOL = 1.0e-12_dp

    character(len=4096) :: csv_path
    integer :: arglen, argstat, nfail, total

    call get_command_argument(1, csv_path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") "usage: test_fortnum_bspline_oracle <bspline.csv>"
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
    write (*, "(a,i0,a)") "oracle passed: ", total, " basis/derivative rows verified"
    stop 0

contains

    subroutine run_all(path, nfail, total)
        character(len=*), intent(in) :: path
        integer, intent(inout)       :: nfail, total

        integer, parameter :: MAXTOK = 256
        type(bspline_workspace_t) :: ws
        type(fortnum_status_t)    :: s
        character(len=4096)       :: line, buf
        character(len=16)         :: tag
        real(dp)                  :: tok(MAXTOK)
        real(dp), allocatable     :: breakpts(:), dvals(:, :), ref(:)
        integer  :: unit, ios, ntok, ci, order, nbreak, ncoef, d, j
        integer  :: cur_order, cur_nbreak
        real(dp) :: x, worst
        logical  :: have_ws

        cur_order  = 0
        cur_nbreak = 0
        have_ws = .false.

        open (newunit=unit, file=path, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            write (error_unit, "(a)") "FAIL: cannot open " // trim(path)
            nfail = nfail + 1
            return
        end if

        do
            read (unit, "(a)", iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (is_comment(line)) cycle

            call split_tag(line, tag, buf)
            call parse_reals(buf, tok, ntok)

            select case (trim(tag))
            case ("brk")
                ! tok: case, order, nbreak, b_0..b_{nbreak-1}
                ci     = nint(tok(1))
                order  = nint(tok(2))
                nbreak = nint(tok(3))
                if (allocated(breakpts)) deallocate (breakpts)
                allocate (breakpts(nbreak))
                do j = 1, nbreak
                    breakpts(j) = tok(3 + j)
                end do
                call bspline_init(ws, order, nbreak, s)
                if (s%code /= FORTNUM_OK) then
                    write (error_unit, "(a,i0,a,a)") "FAIL init case ", ci, " ", trim(s%msg)
                    nfail = nfail + 1
                    cycle
                end if
                call bspline_set_knots(ws, breakpts, s)
                if (s%code /= FORTNUM_OK) then
                    write (error_unit, "(a,i0,a,a)") "FAIL set_knots case ", ci, " ", trim(s%msg)
                    nfail = nfail + 1
                    cycle
                end if
                cur_order  = order
                cur_nbreak = nbreak
                have_ws = .true.

            case ("val")
                ! tok: case, order, nbreak, x, deriv, ncoef, v_0..v_{ncoef-1}
                if (.not. have_ws) then
                    write (error_unit, "(a)") "FAIL: val row before brk row"
                    nfail = nfail + 1
                    cycle
                end if
                ci     = nint(tok(1))
                order  = nint(tok(2))
                nbreak = nint(tok(3))
                x      = tok(4)
                d      = nint(tok(5))
                ncoef  = nint(tok(6))
                if (order /= cur_order .or. nbreak /= cur_nbreak) then
                    write (error_unit, "(a,i0)") "FAIL: val/ws mismatch case ", ci
                    nfail = nfail + 1
                    cycle
                end if
                if (ncoef /= ws%ncoef) then
                    write (error_unit, "(a,i0,a,i0,a,i0)") "FAIL: ncoef case ", ci, &
                        " got ", ws%ncoef, " want ", ncoef
                    nfail = nfail + 1
                    cycle
                end if

                if (allocated(ref)) deallocate (ref)
                allocate (ref(ncoef))
                do j = 1, ncoef
                    ref(j) = tok(6 + j)
                end do

                if (allocated(dvals)) deallocate (dvals)
                allocate (dvals(0:d, ncoef))
                call bspline_eval_deriv(ws, x, d, dvals, s)
                if (s%code /= FORTNUM_OK) then
                    write (error_unit, "(a,i0,a,a)") "FAIL eval_deriv case ", ci, " ", trim(s%msg)
                    nfail = nfail + 1
                    total = total + 1
                    cycle
                end if

                worst = 0.0_dp
                do j = 1, ncoef
                    worst = max(worst, abs(dvals(d, j) - ref(j)))
                end do
                total = total + 1
                if (worst > ATOL) then
                    nfail = nfail + 1
                    write (error_unit, "(a,i0,a,i0,a,es12.4,a,f10.5)") &
                        "FAIL case ", ci, " deriv ", d, " worst_abserr=", worst, &
                        " x=", x
                end if

                ! For deriv 0 also cross-check bspline_eval_basis returns row 0.
                if (d == 0) call check_basis_matches(ws, x, ref, ATOL, ci, nfail, total)

            case default
                write (error_unit, "(a,a)") "FAIL: unknown row tag ", trim(tag)
                nfail = nfail + 1
            end select
        end do
        close (unit)
    end subroutine run_all

    subroutine check_basis_matches(ws, x, ref, atol, ci, nfail, total)
        type(bspline_workspace_t), intent(in) :: ws
        real(dp), intent(in)    :: x, ref(:), atol
        integer,  intent(in)    :: ci
        integer,  intent(inout) :: nfail, total
        type(fortnum_status_t)  :: s
        real(dp), allocatable   :: vals(:)
        real(dp) :: worst
        integer  :: j

        allocate (vals(ws%ncoef))
        call bspline_eval_basis(ws, x, vals, s)
        total = total + 1
        if (s%code /= FORTNUM_OK) then
            write (error_unit, "(a,i0)") "FAIL eval_basis case ", ci
            nfail = nfail + 1
            return
        end if
        worst = 0.0_dp
        do j = 1, ws%ncoef
            worst = max(worst, abs(vals(j) - ref(j)))
        end do
        if (worst > atol) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,es12.4)") &
                "FAIL eval_basis vs deriv-row0 case ", ci, " worst_abserr=", worst
        end if
    end subroutine check_basis_matches

    ! ------------------------------------------------------------------ utils

    pure logical function is_comment(line)
        character(*), intent(in) :: line
        integer :: p
        is_comment = .false.
        p = verify(line, " ")
        if (p > 0) is_comment = (line(p:p) == "#")
    end function is_comment

    ! Split the leading "tag," token off a row; rest is the comma list.
    subroutine split_tag(line, tag, rest)
        character(*), intent(in)  :: line
        character(*), intent(out) :: tag, rest
        integer :: p
        p = index(line, ",")
        if (p <= 0) then
            tag  = adjustl(line)
            rest = ""
        else
            tag  = adjustl(line(1:p-1))
            rest = line(p+1:)
        end if
    end subroutine split_tag

    ! Parse a comma-separated list of reals from buf into tok(1:ntok).
    ! Counts whitespace-delimited fields, then reads them in one list-directed
    ! read so the cost is linear in the row length.
    subroutine parse_reals(buf, tok, ntok)
        character(*), intent(in)  :: buf
        real(dp),     intent(out) :: tok(:)
        integer,      intent(out) :: ntok
        character(len=4096) :: work
        integer :: i, ios
        logical :: in_field
        work = buf
        do i = 1, len(work)
            if (work(i:i) == ",") work(i:i) = " "
        end do
        ! count fields
        ntok = 0
        in_field = .false.
        do i = 1, len_trim(work)
            if (work(i:i) /= " ") then
                if (.not. in_field) then
                    ntok = ntok + 1
                    in_field = .true.
                end if
            else
                in_field = .false.
            end if
        end do
        if (ntok > size(tok)) ntok = size(tok)
        if (ntok == 0) return
        read (work, *, iostat=ios) tok(1:ntok)
        if (ios /= 0) ntok = 0
    end subroutine parse_reals

end program test_fortnum_bspline_oracle
