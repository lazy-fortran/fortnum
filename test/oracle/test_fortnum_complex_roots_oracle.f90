program test_fortnum_complex_roots_oracle
    ! Oracle tests for the complex analytic-zero finder against independent
    ! references (numpy.roots polynomials / mpmath sin zeros).
    !
    ! Usage:
    !   test_fortnum_complex_roots_oracle <complex_roots.csv>
    !
    ! For each case the finder is run over the box; the test asserts the
    ! distinct-zero count, each zero (matched to a reference within newtonz),
    ! its multiplicity, the total winding number, and |f| at each zero.  One
    ! case is re-run with a reduced m_max to exercise the ICON=3 subdivision
    ! recursion, and one case re-runs the independent winding-number check.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_roots_complex, only: complex_region_roots
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    ! The active case selects which f the dispatcher evaluates.  Set before
    ! every call to complex_region_roots; read only inside the f wrappers.
    integer, save :: active_case = -1

    real(dp), parameter :: ZTOL = 5.0e-8_dp ! newtonz default
    real(dp), parameter :: FTOL = 1.0e-7_dp ! |f| acceptance at a polished zero

    character(len=4096) :: path
    integer :: alen, astat, nfail, total

    call get_command_argument(1, path, alen, astat)
    if (astat /= 0) then
        write (error_unit, "(a)") &
            "usage: test_fortnum_complex_roots_oracle <complex_roots.csv>"
        stop 1
    end if

    nfail = 0
    total = 0
    call run_table(trim(path), nfail, total)
    call run_subdivision(nfail, total)

    if (nfail > 0) then
        write (error_unit, "(i0,a,i0,a)") nfail, &
            " oracle check(s) failed out of ", total, " checked"
        stop 1
    end if
    write (*, "(a,i0,a)") "oracle passed: ", total, " checks verified"
    stop 0

contains

    ! ------------------------------------------------------------- f dispatch

    ! Case-indexed analytic functions.  Polynomial cases are evaluated from
    ! their construction roots so the Fortran f matches the reference exactly;
    ! case 4/5 are sin(z).
    subroutine f_case(kr, fk, ctx)
        complex(dp), intent(in)  :: kr
        complex(dp), intent(out) :: fk
        class(*),    intent(in), optional :: ctx

        select case (active_case)
        case (0)
            fk = (kr - cmplx(0.5_dp, 0.5_dp, dp)) &
                * (kr - cmplx(-0.6_dp, -0.3_dp, dp))
        case (1)
            fk = (kr - cmplx(0.4_dp, 0.7_dp, dp))**2 &
                * (kr - cmplx(-0.8_dp, 0.1_dp, dp))
        case (2)
            fk = (kr - cmplx(-0.3_dp, 0.6_dp, dp))**3
        case (3)
            fk = (kr - cmplx(0.20_dp, 0.10_dp, dp)) &
                * (kr - cmplx(0.50_dp, 0.10_dp, dp)) &
                * (kr - cmplx(-0.9_dp, -0.6_dp, dp))**2
        case (4, 5)
            fk = sin(kr)
        case default
            fk = (0.0_dp, 0.0_dp)
        end select
        if (present(ctx)) continue
    end subroutine f_case

    ! ------------------------------------------------------------- table run

    subroutine run_table(p, nfail, total)
        character(len=*), intent(in)    :: p
        integer,          intent(inout) :: nfail, total

        integer :: u, ios, cid, ntot, nd, k
        character(len=512)  :: line
        character(len=64)   :: kind
        real(dp) :: llr, lli, urr, uri
        real(dp), allocatable :: rre(:), rim(:)
        integer,  allocatable :: rmul(:)

        open (newunit=u, file=p, status="old", action="read", iostat=ios)
        if (ios /= 0) then
            write (error_unit, "(a,a)") "cannot open ", p
            stop 1
        end if
        do
            read (u, "(a)", iostat=ios) line
            if (ios /= 0) exit
            if (len_trim(line) == 0) cycle
            if (line(1:1) == "#") cycle
            call parse_case(line, cid, kind, llr, lli, urr, uri, ntot, nd, &
                rre, rim, rmul)
            active_case = cid
            call check_case(cid, kind, &
                cmplx(llr, lli, dp), cmplx(urr, uri, dp), ntot, nd, &
                rre, rim, rmul, 5, nfail, total)
        end do
        close (u)
    end subroutine run_table

    ! Run one case through the finder and assert against its reference triples.
    subroutine check_case(cid, kind, ll, ur, ntot, nd, rre, rim, rmul, mmax, &
            nfail, total)
        integer,          intent(in)    :: cid, ntot, nd, mmax
        character(len=*), intent(in)    :: kind
        complex(dp),      intent(in)    :: ll, ur
        real(dp),         intent(in)    :: rre(:), rim(:)
        integer,          intent(in)    :: rmul(:)
        integer,          intent(inout) :: nfail, total

        complex(dp), allocatable :: roots(:), fvals(:)
        integer,     allocatable :: mult(:)
        integer :: nfound, j, mret, sumret
        type(fortnum_status_t) :: status
        logical :: ok

        call complex_region_roots(f_case, ll, ur, roots, fvals, mult, &
            nfound, status, m_max=mmax)

        total = total + 1
        if (status%code /= FORTNUM_OK) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,i0,a,a)") "FAIL case ", cid, &
                " status=", status%code, " ", trim(status%msg)
            return
        end if

        ! Distinct-zero count.
        total = total + 1
        if (nfound /= nd) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,i0,a,i0)") "FAIL case ", cid, &
                " ndistinct got=", nfound, " expected=", nd
        end if

        ! Total winding number (sum of recovered multiplicities).
        sumret = 0
        do j = 1, nfound
            sumret = sumret + mult(j)
        end do
        total = total + 1
        if (sumret /= ntot) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,i0,a,i0)") "FAIL case ", cid, &
                " ntotal got=", sumret, " expected=", ntot
        end if

        ! Match each reference zero to a returned zero within ZTOL.
        do j = 1, nd
            total = total + 1
            ok = match_root(cmplx(rre(j), rim(j), dp), roots, fvals, mult, &
                nfound, rmul(j), mret)
            if (.not. ok) then
                nfail = nfail + 1
                write (error_unit, "(a,i0,a,a,a,es12.4,a,es12.4,a)") &
                    "FAIL case ", cid, " (", trim(kind), &
                    ") missing zero (", rre(j), ",", rim(j), ")"
            end if
        end do
    end subroutine check_case

    ! Find a returned zero within ZTOL of zref; assert its multiplicity and
    ! that |f| there is below FTOL.  Returns the matched multiplicity in mret.
    logical function match_root(zref, roots, fvals, mult, nfound, mexp, mret)
        complex(dp), intent(in) :: zref, roots(:), fvals(:)
        integer,     intent(in) :: mult(:), nfound, mexp
        integer,     intent(out):: mret
        integer :: j, jbest
        real(dp) :: d, dbest

        match_root = .false.
        mret = 0
        jbest = 0
        dbest = huge(1.0_dp)
        do j = 1, nfound
            d = abs(roots(j) - zref)
            if (d < dbest) then
                dbest = d
                jbest = j
            end if
        end do
        if (jbest == 0) return
        if (dbest > ZTOL) return
        if (abs(fvals(jbest)) > FTOL) return
        mret = mult(jbest)
        match_root = (mret == mexp)
    end function match_root

    ! ------------------------------------------------------- subdivision run

    ! Re-run sin(z) over the wide strip (3 zeros) with m_max=2, which forces
    ! the box to bisect (ICON=3 recursion) until each cell holds <= 2 zeros.
    ! The recovered zeros must still match -pi, 0, pi.
    subroutine run_subdivision(nfail, total)
        integer, intent(inout) :: nfail, total

        complex(dp), allocatable :: roots(:), fvals(:)
        integer,     allocatable :: mult(:)
        integer :: nfound, j
        type(fortnum_status_t) :: status
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        logical :: hit_m, hit_0, hit_p

        active_case = 4
        call complex_region_roots(f_case, cmplx(-4.0_dp, -1.0_dp, dp), &
            cmplx(4.0_dp, 1.0_dp, dp), roots, fvals, mult, nfound, status, &
            m_max=2)

        total = total + 1
        if (status%code /= FORTNUM_OK .or. nfound /= 3) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,i0)") &
                "FAIL subdivision status=", status%code, " nfound=", nfound
            return
        end if

        hit_m = .false.; hit_0 = .false.; hit_p = .false.
        do j = 1, nfound
            if (abs(roots(j) - cmplx(-pi, 0.0_dp, dp)) < ZTOL) hit_m = .true.
            if (abs(roots(j) - cmplx(0.0_dp, 0.0_dp, dp)) < ZTOL) hit_0 = .true.
            if (abs(roots(j) - cmplx(pi, 0.0_dp, dp)) < ZTOL) hit_p = .true.
        end do
        total = total + 1
        if (.not. (hit_m .and. hit_0 .and. hit_p)) then
            nfail = nfail + 1
            write (error_unit, "(a)") &
                "FAIL subdivision: did not recover -pi, 0, pi via bisection"
        end if
    end subroutine run_subdivision

    ! ----------------------------------------------------------- CSV parsing

    subroutine parse_case(line, cid, kind, llr, lli, urr, uri, ntot, nd, &
            rre, rim, rmul)
        character(len=*), intent(in)  :: line
        integer,          intent(out) :: cid, ntot, nd
        character(len=*), intent(out) :: kind
        real(dp),         intent(out) :: llr, lli, urr, uri
        real(dp), allocatable, intent(out) :: rre(:), rim(:)
        integer,  allocatable, intent(out) :: rmul(:)

        character(len=64) :: tok(256)
        integer :: ntok, j

        call split_commas(line, tok, ntok)
        read (tok(1), *) cid
        kind = trim(tok(2))
        read (tok(3), *) llr
        read (tok(4), *) lli
        read (tok(5), *) urr
        read (tok(6), *) uri
        read (tok(7), *) ntot
        read (tok(8), *) nd
        allocate(rre(nd), rim(nd), rmul(nd))
        do j = 1, nd
            read (tok(8 + 3*(j-1) + 1), *) rre(j)
            read (tok(8 + 3*(j-1) + 2), *) rim(j)
            read (tok(8 + 3*(j-1) + 3), *) rmul(j)
        end do
    end subroutine parse_case

    subroutine split_commas(line, tok, ntok)
        character(len=*), intent(in)  :: line
        character(len=*), intent(out) :: tok(:)
        integer,          intent(out) :: ntok
        integer :: i, start, n

        ntok = 0
        start = 1
        n = len_trim(line)
        do i = 1, n
            if (line(i:i) == ",") then
                ntok = ntok + 1
                tok(ntok) = line(start:i-1)
                start = i + 1
            end if
        end do
        ntok = ntok + 1
        tok(ntok) = line(start:n)
    end subroutine split_commas

end program test_fortnum_complex_roots_oracle
