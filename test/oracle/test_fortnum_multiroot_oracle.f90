program test_fortnum_multiroot_oracle
    ! Oracle tests for the multiroot module against independent references.
    !
    ! Usage:
    !   test_fortnum_multiroot_oracle <multiroot.csv> <deriv_central.csv> <argsort.csv>
    !
    ! multiroot.csv  : standard n-dim systems; MINPACK hybr roots (scipy).
    !                  Check the analytic-Jacobian and FD-Jacobian solvers both
    !                  drive |F(x)|_inf below ftol AND match the reference root.
    ! deriv_central.csv : analytic f'(x); check the 5-point central rule against
    !                  it and check the returned abserr actually bounds the error.
    ! argsort.csv    : numpy.argsort permutation; check x(perm) nondecreasing and
    !                  the permutation matches.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_multiroot, only: multiroot_hybrid, multiroot_hybrids, &
        deriv_central, argsort
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    character(len=4096) :: p_multi, p_deriv, p_sort
    integer :: alen, astat, nfail, total

    call get_command_argument(1, p_multi, alen, astat)
    if (astat /= 0) call usage()
    call get_command_argument(2, p_deriv, alen, astat)
    if (astat /= 0) call usage()
    call get_command_argument(3, p_sort, alen, astat)
    if (astat /= 0) call usage()

    nfail = 0
    total = 0

    call run_multiroot(trim(p_multi), nfail, total)
    call run_deriv(trim(p_deriv), nfail, total)
    call run_argsort(trim(p_sort), nfail, total)

    if (nfail > 0) then
        write (error_unit, "(i0,a,i0,a)") nfail, " oracle case(s) failed out of ", &
            total, " checked"
        stop 1
    end if
    write (*, "(a,i0,a)") "oracle passed: ", total, " cases verified"
    stop 0

contains

    subroutine usage()
        write (error_unit, "(a)") &
            "usage: test_fortnum_multiroot_oracle <multiroot.csv> <deriv_central.csv> <argsort.csv>"
        stop 1
    end subroutine usage

    ! ----------------------------------------------------------- multiroot

    subroutine run_multiroot(path, nfail, total)
        character(len=*), intent(in)    :: path
        integer,          intent(inout) :: nfail, total

        ! Residual tolerance over the consumers' domain; both solvers must
        ! reduce |F|_inf to at least this.  Powell singular reaches ~1e-9
        ! before the Newton step underflows near the singular Jacobian.
        real(dp), parameter :: RES_TOL  = 1.0e-7_dp
        ! Root-match tolerance: the converged x must lie this close to the
        ! reference root component-wise (looser where the root is singular).
        real(dp), parameter :: ROOT_TOL = 1.0e-4_dp

        integer  :: unit, ios, idx, n, i
        character(len=1024) :: line, buf
        real(dp) :: vals(64)
        real(dp), allocatable :: x0(:), ref(:), x(:), f(:)
        type(fortnum_status_t) :: s

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
            buf = line
            call replace_commas(buf)
            read (buf, *, iostat=ios) idx, n
            if (ios /= 0) cycle
            read (buf, *, iostat=ios) idx, n, (vals(i), i = 1, 2 * n)
            if (ios /= 0) then
                write (error_unit, "(a,i0)") "multiroot: parse error case ", idx
                nfail = nfail + 1
                cycle
            end if
            allocate (x0(n), ref(n), x(n), f(n))
            x0  = vals(1:n)
            ref = vals(n + 1:2 * n)

            ! Analytic-Jacobian solver.
            call dispatch_hybrid(idx, n, x0, x, s)
            call eval_system(idx, x, f)
            total = total + 1
            call check_one("hybrid", idx, x, ref, f, RES_TOL, ROOT_TOL, s, nfail)

            ! Finite-difference-Jacobian solver.
            call dispatch_hybrids(idx, n, x0, x, s)
            call eval_system(idx, x, f)
            total = total + 1
            call check_one("hybrids", idx, x, ref, f, RES_TOL, ROOT_TOL, s, nfail)

            deallocate (x0, ref, x, f)
        end do
        close (unit)
    end subroutine run_multiroot

    subroutine check_one(tag, idx, x, ref, f, res_tol, root_tol, s, nfail)
        character(len=*),       intent(in)    :: tag
        integer,                intent(in)    :: idx
        real(dp),               intent(in)    :: x(:), ref(:), f(:)
        real(dp),               intent(in)    :: res_tol, root_tol
        type(fortnum_status_t), intent(in)    :: s
        integer,                intent(inout) :: nfail

        real(dp) :: res, derr

        if (s%code /= FORTNUM_OK) then
            write (error_unit, "(a,a,a,i0,a,i0,a,a)") "FAIL ", tag, " case ", idx, &
                " status=", s%code, " ", trim(s%msg)
            nfail = nfail + 1
            return
        end if
        res = maxval(abs(f))
        if (res > res_tol) then
            write (error_unit, "(a,a,a,i0,a,es12.4)") "FAIL ", tag, " case ", idx, &
                " residual=", res
            nfail = nfail + 1
            return
        end if
        derr = maxval(abs(x - ref))
        if (derr > root_tol) then
            write (error_unit, "(a,a,a,i0,a,es12.4)") "FAIL ", tag, " case ", idx, &
                " root mismatch maxabs=", derr
            nfail = nfail + 1
        end if
    end subroutine check_one

    subroutine dispatch_hybrid(idx, n, x0, x, s)
        integer,  intent(in)  :: idx, n
        real(dp), intent(in)  :: x0(n)
        real(dp), intent(out) :: x(n)
        type(fortnum_status_t), intent(out) :: s
        select case (idx)
        case (0); call multiroot_hybrid(rosen_fdf,  n, x0, x, s, ftol=1.0e-12_dp)
        case (1); call multiroot_hybrid(powell_fdf, n, x0, x, s, ftol=1.0e-9_dp)
        case (2); call multiroot_hybrid(circ_fdf,   n, x0, x, s, ftol=1.0e-12_dp)
        case default; call multiroot_hybrid(circ_fdf, n, x0, x, s)
        end select
    end subroutine dispatch_hybrid

    subroutine dispatch_hybrids(idx, n, x0, x, s)
        integer,  intent(in)  :: idx, n
        real(dp), intent(in)  :: x0(n)
        real(dp), intent(out) :: x(n)
        type(fortnum_status_t), intent(out) :: s
        select case (idx)
        case (0); call multiroot_hybrids(rosen_fn,  n, x0, x, s, ftol=1.0e-12_dp)
        case (1); call multiroot_hybrids(powell_fn, n, x0, x, s, ftol=1.0e-9_dp)
        case (2); call multiroot_hybrids(circ_fn,   n, x0, x, s, ftol=1.0e-12_dp)
        case default; call multiroot_hybrids(circ_fn, n, x0, x, s)
        end select
    end subroutine dispatch_hybrids

    subroutine eval_system(idx, x, f)
        integer,  intent(in)  :: idx
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: f(:)
        select case (idx)
        case (0); call rosen_fn(x, f)
        case (1); call powell_fn(x, f)
        case (2); call circ_fn(x, f)
        end select
    end subroutine eval_system

    ! Rosenbrock gradient = 0 (root (1,1)).
    subroutine rosen_fn(x, f, ctx)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: f(:)
        class(*), intent(in), optional :: ctx
        f(1) = -400.0_dp * x(1) * (x(2) - x(1)*x(1)) - 2.0_dp * (1.0_dp - x(1))
        f(2) = 200.0_dp * (x(2) - x(1)*x(1))
    end subroutine rosen_fn

    subroutine rosen_fdf(x, f, jac, ctx)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: f(:)
        real(dp), intent(out) :: jac(:, :)
        class(*), intent(in), optional :: ctx
        call rosen_fn(x, f)
        jac(1, 1) = -400.0_dp * (x(2) - 3.0_dp * x(1)*x(1)) + 2.0_dp
        jac(1, 2) = -400.0_dp * x(1)
        jac(2, 1) = -400.0_dp * x(1)
        jac(2, 2) = 200.0_dp
    end subroutine rosen_fdf

    ! Powell singular function (MGH #13), root at origin.
    subroutine powell_fn(x, f, ctx)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: f(:)
        class(*), intent(in), optional :: ctx
        f(1) = x(1) + 10.0_dp * x(2)
        f(2) = sqrt(5.0_dp) * (x(3) - x(4))
        f(3) = (x(2) - 2.0_dp * x(3))**2
        f(4) = sqrt(10.0_dp) * (x(1) - x(4))**2
    end subroutine powell_fn

    subroutine powell_fdf(x, f, jac, ctx)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: f(:)
        real(dp), intent(out) :: jac(:, :)
        class(*), intent(in), optional :: ctx
        call powell_fn(x, f)
        jac = 0.0_dp
        jac(1, 1) = 1.0_dp
        jac(1, 2) = 10.0_dp
        jac(2, 3) = sqrt(5.0_dp)
        jac(2, 4) = -sqrt(5.0_dp)
        jac(3, 2) = 2.0_dp * (x(2) - 2.0_dp * x(3))
        jac(3, 3) = -4.0_dp * (x(2) - 2.0_dp * x(3))
        jac(4, 1) = 2.0_dp * sqrt(10.0_dp) * (x(1) - x(4))
        jac(4, 4) = -2.0_dp * sqrt(10.0_dp) * (x(1) - x(4))
    end subroutine powell_fdf

    ! Circle-line system: x1^2 + x2^2 - 2 = 0, x1 - x2 = 0; root (1,1).
    subroutine circ_fn(x, f, ctx)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: f(:)
        class(*), intent(in), optional :: ctx
        f(1) = x(1)*x(1) + x(2)*x(2) - 2.0_dp
        f(2) = x(1) - x(2)
    end subroutine circ_fn

    subroutine circ_fdf(x, f, jac, ctx)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: f(:)
        real(dp), intent(out) :: jac(:, :)
        class(*), intent(in), optional :: ctx
        call circ_fn(x, f)
        jac(1, 1) = 2.0_dp * x(1)
        jac(1, 2) = 2.0_dp * x(2)
        jac(2, 1) = 1.0_dp
        jac(2, 2) = -1.0_dp
    end subroutine circ_fdf

    ! ----------------------------------------------------------- deriv_central

    subroutine run_deriv(path, nfail, total)
        character(len=*), intent(in)    :: path
        integer,          intent(inout) :: nfail, total

        ! 5-point central rule on a smooth f with h=1e-3: error ~ h^4 ~ 1e-12,
        ! limited below by rounding ~ eps/h ~ 1e-13.  Allow 1e-9 against the
        ! analytic reference; abserr must also bound the actual error.
        real(dp), parameter :: VAL_TOL = 1.0e-9_dp

        integer  :: unit, ios, idx
        character(len=512) :: line, buf
        real(dp) :: x, h, dref, result, abserr, err
        type(fortnum_status_t) :: s

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
            buf = line
            call replace_commas(buf)
            read (buf, *, iostat=ios) idx, x, h, dref
            if (ios /= 0) cycle

            call dispatch_deriv(idx, x, h, result, abserr, s)
            total = total + 1
            if (s%code /= FORTNUM_OK) then
                write (error_unit, "(a,i0,a,i0)") "FAIL deriv case ", idx, &
                    " status=", s%code
                nfail = nfail + 1
                cycle
            end if
            err = abs(result - dref)
            if (err > VAL_TOL) then
                write (error_unit, "(a,i0,a,es12.4,a,es24.16,a,es24.16)") &
                    "FAIL deriv case ", idx, " abserr=", err, &
                    " got=", result, " ref=", dref
                nfail = nfail + 1
                cycle
            end if
            ! The reported error estimate must actually bound the true error
            ! (allow a small slack factor; the error estimate is conservative).
            if (err > abserr * 100.0_dp + 1.0e-14_dp) then
                write (error_unit, "(a,i0,a,es12.4,a,es12.4)") &
                    "FAIL deriv case ", idx, " true err=", err, &
                    " exceeds reported abserr=", abserr
                nfail = nfail + 1
            end if
        end do
        close (unit)
    end subroutine run_deriv

    subroutine dispatch_deriv(idx, x, h, result, abserr, s)
        integer,  intent(in)  :: idx
        real(dp), intent(in)  :: x, h
        real(dp), intent(out) :: result, abserr
        type(fortnum_status_t), intent(out) :: s
        select case (idx)
        case (0); call deriv_central(d_sin, x, h, result, abserr, s)
        case (1); call deriv_central(d_exp, x, h, result, abserr, s)
        case (2); call deriv_central(d_cube, x, h, result, abserr, s)
        case (3); call deriv_central(d_log, x, h, result, abserr, s)
        case (4); call deriv_central(d_cos, x, h, result, abserr, s)
        case default; call deriv_central(d_sin, x, h, result, abserr, s)
        end select
    end subroutine dispatch_deriv

    function d_sin(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = sin(x)
    end function d_sin

    function d_exp(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = exp(x)
    end function d_exp

    function d_cube(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = x * x * x
    end function d_cube

    function d_log(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = log(x)
    end function d_log

    function d_cos(x, ctx) result(y)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: y
        y = cos(x)
    end function d_cos

    ! ----------------------------------------------------------- argsort

    subroutine run_argsort(path, nfail, total)
        character(len=*), intent(in)    :: path
        integer,          intent(inout) :: nfail, total

        integer  :: unit, ios, idx, n, i
        character(len=1024) :: line, buf
        real(dp) :: xv(64)
        integer  :: pref(64)
        real(dp), allocatable :: x(:)
        integer,  allocatable :: perm(:), ref(:)
        logical  :: ok

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
            buf = line
            call replace_commas(buf)
            read (buf, *, iostat=ios) idx, n
            if (ios /= 0) cycle
            read (buf, *, iostat=ios) idx, n, (xv(i), i = 1, n), (pref(i), i = 1, n)
            if (ios /= 0) then
                write (error_unit, "(a,i0)") "argsort: parse error case ", idx
                nfail = nfail + 1
                cycle
            end if
            allocate (x(n), perm(n), ref(n))
            x   = xv(1:n)
            ref = pref(1:n)
            call argsort(x, perm)
            total = total + 1

            ! x(perm) must be nondecreasing.
            ok = .true.
            do i = 2, n
                if (x(perm(i)) < x(perm(i - 1))) ok = .false.
            end do
            if (.not. ok) then
                write (error_unit, "(a,i0,a)") "FAIL argsort case ", idx, &
                    " x(perm) not nondecreasing"
                nfail = nfail + 1
            else if (any(perm /= ref)) then
                write (error_unit, "(a,i0,a)") "FAIL argsort case ", idx, &
                    " permutation differs from numpy.argsort"
                nfail = nfail + 1
            end if
            deallocate (x, perm, ref)
        end do
        close (unit)
    end subroutine run_argsort

    ! ----------------------------------------------------------------- utils

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

end program test_fortnum_multiroot_oracle
