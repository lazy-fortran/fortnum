program test_fortnum_multiroot_hybrids_minpack
    ! multiroot_hybrids on problems of the MINPACK hybrd test set (More,
    ! Garbow & Hillstrom, ACM TOMS 7 (1981) 17; problem numbers of the
    ! MINPACK-1 driver hybdrv in brackets), standard starting points.
    !
    ! Oracles: the residual |F|_inf at the returned point, the known root
    ! where it is closed-form, and the function-evaluation count of MINPACK
    ! hybrd1 itself.  MINPACK_NFEV is the evaluation at which hybrd1 (unit
    ! scaling, factor 100, forward-difference Jacobian) first reaches the
    ! residual tolerance used here; hybdrv's published totals at its own
    ! stopping test (xtol = sqrt(eps)) are 22, 108, 27, 21 and 16.  A wrong
    ! dogleg, trust-region update, Broyden update or Jacobian-refresh rule
    ! changes the trajectory and exceeds the ceiling of MINPACK_NFEV plus 25%.
    ! OLD_NFEV are the counts of the previous Newton/line-search
    ! implementation with a central-difference Jacobian at every iterate
    ! (fortnum f856c80); none may regress.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_multiroot, only: multiroot_hybrids
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    integer, parameter :: ROSENBROCK = 1, POWELL_SINGULAR = 2, &
        HELICAL_VALLEY = 3, BROYDEN_TRIDIAGONAL = 4, DISCRETE_BVP = 5, &
        NO_ROOT = 6
    integer, parameter :: OLD_ROSENBROCK = 78, OLD_POWELL = 154, &
        OLD_HELICAL = 72, OLD_TRIDIAGONAL = 106, OLD_BVP = 64
    integer :: nfail, nfev

    nfail = 0
    ! hybdrv problem 1
    call check_case("rosenbrock", ROSENBROCK, [-1.2_dp, 1.0_dp], &
        [1.0_dp, 1.0_dp], 1.0e-12_dp, 22, OLD_ROSENBROCK)
    ! hybdrv problem 2; the root is singular, hence the looser ftol
    call check_case("powell singular", POWELL_SINGULAR, &
        [3.0_dp, -1.0_dp, 0.0_dp, 1.0_dp], [0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp], &
        1.0e-9_dp, 29, OLD_POWELL)
    ! hybdrv problem 5
    call check_case("helical valley", HELICAL_VALLEY, [-1.0_dp, 0.0_dp, 0.0_dp], &
        [1.0_dp, 0.0_dp, 0.0_dp], 1.0e-12_dp, 27, OLD_HELICAL)
    ! hybdrv problem 13, n = 10
    call check_case("broyden tridiagonal", BROYDEN_TRIDIAGONAL, &
        spread(-1.0_dp, 1, 10), [real(dp) ::], 1.0e-12_dp, 25, OLD_TRIDIAGONAL)
    ! hybdrv problem 9, n = 10
    call check_case("discrete boundary value", DISCRETE_BVP, &
        bvp_start(10), [real(dp) ::], 1.0e-12_dp, 16, OLD_BVP)
    call check_no_root()

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " check(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine check_case(name, problem, x0, root, ftol, minpack_nfev, old_nfev)
        character(len=*), intent(in) :: name
        integer, intent(in) :: problem, minpack_nfev, old_nfev
        real(dp), intent(in) :: x0(:), root(:), ftol

        real(dp), parameter :: ROOT_TOL = 1.0e-4_dp
        real(dp) :: x(size(x0)), f(size(x0))
        type(fortnum_status_t) :: status
        integer :: n, ceiling, used

        n = size(x0)
        ceiling = (5*minpack_nfev + 3)/4
        nfev = 0
        call multiroot_hybrids(system_fn, n, x0, x, status, ftol=ftol, &
            ctx=problem)
        used = nfev
        call system_fn(x, f, problem)
        write (*, "(a,a,i0,a,i0,a,i0)") name, ": evaluations=", used, &
            " minpack=", minpack_nfev, " old=", old_nfev
        if (status%code /= FORTNUM_OK) then
            write (error_unit, "(a,a,a)") "FAIL status: ", name, &
                " "//trim(status%msg)
            nfail = nfail + 1
        end if
        if (maxval(abs(f)) > ftol) then
            write (error_unit, "(a,a,es12.4)") "FAIL residual: ", &
                name//" ", maxval(abs(f))
            nfail = nfail + 1
        end if
        if (size(root) == n) then
            if (maxval(abs(x - root)) > ROOT_TOL) then
                write (error_unit, "(a,a,es12.4)") "FAIL root: ", &
                    name, maxval(abs(x - root))
                nfail = nfail + 1
            end if
        end if
        if (used > ceiling) then
            write (error_unit, "(a,a,a,i0,a,i0)") "FAIL evaluations vs MINPACK: ", &
                name, " ", used, " > ", ceiling
            nfail = nfail + 1
        end if
        if (used > old_nfev) then
            write (error_unit, "(a,a,a,i0,a,i0)") "FAIL evaluations vs old: ", &
                name, " ", used, " > ", old_nfev
            nfail = nfail + 1
        end if
    end subroutine check_case

    ! F(x) = x^2 + 1 has no real root: the solver must report failure from
    ! its lack-of-progress tests long before the default iteration cap.
    subroutine check_no_root()
        real(dp) :: x(1)
        type(fortnum_status_t) :: status

        nfev = 0
        call multiroot_hybrids(system_fn, 1, [1.0_dp], x, status, ctx=NO_ROOT)
        write (*, "(a,i0,a,a)") "no root: evaluations=", nfev, " ", trim(status%msg)
        if (status%code == FORTNUM_OK) then
            write (error_unit, "(a)") "FAIL no root: reported success"
            nfail = nfail + 1
        end if
        if (nfev > 100) then
            write (error_unit, "(a,i0)") "FAIL no root: evaluations ", nfev
            nfail = nfail + 1
        end if
    end subroutine check_no_root

    pure function bvp_start(n) result(x0)
        integer, intent(in) :: n
        real(dp) :: x0(n)
        integer :: i
        do i = 1, n
            x0(i) = real(i, dp)/real(n + 1, dp)*(real(i, dp)/real(n + 1, dp) &
                - 1.0_dp)
        end do
    end function bvp_start

    subroutine system_fn(x, f, ctx)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: jac(size(x), size(x))
        nfev = nfev + 1
        call evaluate(problem_of(ctx), x, f, jac)
    end subroutine system_fn

    integer function problem_of(ctx)
        class(*), intent(in), optional :: ctx
        problem_of = 0
        if (.not. present(ctx)) return
        select type (ctx)
            type is (integer)
            problem_of = ctx
        end select
    end function problem_of

    pure subroutine evaluate(problem, x, f, jac)
        integer, intent(in) :: problem
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: f(:), jac(:, :)
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: theta, r, h
        integer :: i, n

        jac = 0.0_dp
        n = size(x)
        select case (problem)
        case (ROSENBROCK)
            f(1) = 10.0_dp*(x(2) - x(1)**2)
            f(2) = 1.0_dp - x(1)
            jac(1, :) = [-20.0_dp*x(1), 10.0_dp]
            jac(2, :) = [-1.0_dp, 0.0_dp]
        case (POWELL_SINGULAR)
            f(1) = x(1) + 10.0_dp*x(2)
            f(2) = sqrt(5.0_dp)*(x(3) - x(4))
            f(3) = (x(2) - 2.0_dp*x(3))**2
            f(4) = sqrt(10.0_dp)*(x(1) - x(4))**2
            jac(1, 1:2) = [1.0_dp, 10.0_dp]
            jac(2, 3:4) = [sqrt(5.0_dp), -sqrt(5.0_dp)]
            jac(3, 2:3) = [2.0_dp, -4.0_dp]*(x(2) - 2.0_dp*x(3))
            jac(4, [1, 4]) = [2.0_dp, -2.0_dp]*sqrt(10.0_dp)*(x(1) - x(4))
        case (HELICAL_VALLEY)
            ! MINPACK branch: continuous across x2 = 0 for x1 < 0.
            if (x(1) > 0.0_dp) then
                theta = atan(x(2)/x(1))/(2.0_dp*pi)
            else if (x(1) < 0.0_dp) then
                theta = atan(x(2)/x(1))/(2.0_dp*pi) + 0.5_dp
            else
                theta = sign(0.25_dp, x(2))
            end if
            r = hypot(x(1), x(2))
            f(1) = 10.0_dp*(x(3) - 10.0_dp*theta)
            f(2) = 10.0_dp*(r - 1.0_dp)
            f(3) = x(3)
            jac(1, :) = [50.0_dp*x(2)/(pi*r**2), -50.0_dp*x(1)/(pi*r**2), &
                10.0_dp]
            jac(2, :) = [10.0_dp*x(1)/r, 10.0_dp*x(2)/r, 0.0_dp]
            jac(3, :) = [0.0_dp, 0.0_dp, 1.0_dp]
        case (BROYDEN_TRIDIAGONAL)
            do i = 1, n
                f(i) = (3.0_dp - 2.0_dp*x(i))*x(i) + 1.0_dp
                jac(i, i) = 3.0_dp - 4.0_dp*x(i)
                if (i > 1) then
                    f(i) = f(i) - x(i - 1)
                    jac(i, i - 1) = -1.0_dp
                end if
                if (i < n) then
                    f(i) = f(i) - 2.0_dp*x(i + 1)
                    jac(i, i + 1) = -2.0_dp
                end if
            end do
        case (DISCRETE_BVP)
            h = 1.0_dp/real(n + 1, dp)
            do i = 1, n
                f(i) = 2.0_dp*x(i) + 0.5_dp*h**2*(x(i) + real(i, dp)*h &
                    + 1.0_dp)**3
                jac(i, i) = 2.0_dp + 1.5_dp*h**2*(x(i) + real(i, dp)*h &
                    + 1.0_dp)**2
                if (i > 1) then
                    f(i) = f(i) - x(i - 1)
                    jac(i, i - 1) = -1.0_dp
                end if
                if (i < n) then
                    f(i) = f(i) - x(i + 1)
                    jac(i, i + 1) = -1.0_dp
                end if
            end do
        case (NO_ROOT)
            f(1) = x(1)**2 + 1.0_dp
            jac(1, 1) = 2.0_dp*x(1)
        case default
            f = huge(1.0_dp)
        end select
    end subroutine evaluate

end program test_fortnum_multiroot_hybrids_minpack
