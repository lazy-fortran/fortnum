program test_fortnum_ode_gauss_radau
    ! Behavioural tests for the adaptive 15th-order Gauss-Radau integrator.
    !
    ! Oracles are analytic or structural, never a recording of this integrator's
    ! own output:
    !   - Radau nodes are checked against the defining quadrature property
    !     (a left-Radau rule with n nodes is exact for polynomials up to degree
    !     2n-2), which is independent of how the nodes were computed.
    !   - scalar decay and a nonlinear scalar ODE have closed-form solutions
    !   - Kepler energy is an exact invariant of the true flow
    !   - accuracy is compared against DOP853, an unrelated method

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    use fortnum_ode, only: ode_problem_t, ode_solution_t
    use fortnum_ode_gauss_radau, only: ode_integrate_radau, ode_solve_radau, &
                                       radau_nodes, RADAU_STAGES
    use fortnum_ode_dop853, only: ode_solve_dop
    implicit none

    integer :: nfail

    nfail = 0
    call check_nodes(nfail)
    call check_decay(nfail)
    call check_nonlinear_scalar(nfail)
    call check_kepler_energy(nfail)
    call check_vs_dop853(nfail)
    call check_backward(nfail)
    call check_bad_input(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "fortnum_ode_gauss_radau: all tests passed"

contains

    subroutine rhs_decay(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) = -y(1)
    end subroutine rhs_decay

    ! y' = -y^2, y(0) = 1  =>  y(t) = 1/(1+t). Nonlinear, so it exercises the
    ! corrector rather than collapsing to a linear recurrence.
    subroutine rhs_nl(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt(1) = -y(1)*y(1)
    end subroutine rhs_nl

    subroutine rhs_kepler(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        real(dp) :: r3
        associate (unused_t => t); end associate
        r3 = (sqrt(y(1)**2 + y(2)**2))**3
        dydt(1) = y(3)
        dydt(2) = y(4)
        dydt(3) = -y(1)/r3
        dydt(4) = -y(2)/r3
    end subroutine rhs_kepler

    real(dp) function kepler_energy(y) result(e)
        real(dp), intent(in) :: y(:)
        e = 0.5_dp*(y(3)**2 + y(4)**2) - 1.0_dp/sqrt(y(1)**2 + y(2)**2)
    end function kepler_energy

    ! A left-Radau rule on [0,1] with n nodes (one fixed at 0) is exact for every
    ! polynomial of degree <= 2n-2. Recovering the weights by least squares from
    ! the low-degree moments and then checking the high-degree ones is an
    ! independent statement about the nodes.
    subroutine check_nodes(nfail)
        integer, intent(inout) :: nfail

        integer, parameter :: n = RADAU_STAGES
        real(dp) :: c(n), w(n), a(n, n), rhs(n)
        real(dp) :: quad, exact, worst
        integer :: i, j, d, ok_deg, nfound
        logical :: ok

        call radau_nodes(n, c, nfound)
        if (nfound /= n) then
            write (error_unit, "(a,i0,a,i0)") &
                "  Radau node search found ", nfound, " roots; expected ", n
            nfail = nfail + 1
            return
        end if

        ! Nodes must be distinct, increasing, inside [0,1], starting at 0.
        ok = .true.
        if (abs(c(1)) > 1.0e-15_dp) ok = .false.
        do i = 2, n
            if (.not. (c(i) > c(i - 1))) ok = .false.
            if (c(i) <= 0.0_dp .or. c(i) > 1.0_dp) ok = .false.
        end do
        if (.not. ok) then
            write (error_unit, "(a)") "  radau nodes not increasing within [0,1]"
            write (error_unit, "(8f12.8)") c
            nfail = nfail + 1
            return
        end if

        ! Solve for weights from exactness on degrees 0..n-1.
        do d = 0, n - 1
            do j = 1, n
                a(d + 1, j) = c(j)**d
            end do
            rhs(d + 1) = 1.0_dp/real(d + 1, dp)
        end do
        call solve_dense(a, rhs, w, ok)
        if (.not. ok) then
            write (error_unit, "(a)") "  could not solve for Radau weights"
            nfail = nfail + 1
            return
        end if

        ! Now test the degrees the rule is supposed to get for free.
        worst = 0.0_dp
        ok_deg = -1
        do d = 0, 2*n - 2
            quad = 0.0_dp
            do j = 1, n
                quad = quad + w(j)*c(j)**d
            end do
            exact = 1.0_dp/real(d + 1, dp)
            worst = max(worst, abs(quad - exact))
            if (abs(quad - exact) < 1.0e-12_dp) ok_deg = d
        end do

        write (*, "(a,i0,a,es10.3)") "  radau nodes: exact through degree ", ok_deg, &
            " (need ", 1.0e-12_dp
        if (ok_deg < 2*n - 2) then
            write (error_unit, "(a,i0,a,i0)") &
                "  Radau exactness only to degree ", ok_deg, ", expected ", 2*n - 2
            nfail = nfail + 1
        end if
    end subroutine check_nodes

    subroutine solve_dense(a_in, b_in, x, ok)
        real(dp), intent(in) :: a_in(:, :), b_in(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok

        real(dp), allocatable :: a(:, :), b(:)
        integer :: n, i, k, piv
        real(dp) :: pval, fac

        n = size(b_in)
        allocate (a(n, n), b(n))
        a = a_in
        b = b_in
        ok = .true.
        do k = 1, n
            piv = k
            pval = abs(a(k, k))
            do i = k + 1, n
                if (abs(a(i, k)) > pval) then
                    pval = abs(a(i, k))
                    piv = i
                end if
            end do
            if (pval <= 0.0_dp) then
                ok = .false.
                return
            end if
            if (piv /= k) then
                block
                    real(dp) :: rowtmp(n), btmp
                    rowtmp = a(k, :); a(k, :) = a(piv, :); a(piv, :) = rowtmp
                    btmp = b(k); b(k) = b(piv); b(piv) = btmp
                end block
            end if
            fac = a(k, k)
            a(k, :) = a(k, :)/fac
            b(k) = b(k)/fac
            do i = 1, n
                if (i == k) cycle
                fac = a(i, k)
                if (fac == 0.0_dp) cycle
                a(i, :) = a(i, :) - fac*a(k, :)
                b(i) = b(i) - fac*b(k)
            end do
        end do
        x = b
    end subroutine solve_dense

    subroutine check_decay(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: status
        real(dp) :: y0(1), err
        real(dp), allocatable :: t_out(:), y_out(:, :)

        y0(1) = 1.0_dp
        call ode_solve_radau(rhs_decay, 0.0_dp, 5.0_dp, y0, t_out, y_out, status, &
                             rtol=1.0e-13_dp)
        err = abs(y_out(1, size(t_out)) - exp(-5.0_dp))
        write (*, "(a,es10.3)") "  decay endpoint error ", err
        if (status%code /= FORTNUM_OK .or. .not. (err < 1.0e-12_dp)) then
            write (error_unit, "(a)") "  decay: did not reach exp(-5) accurately"
            nfail = nfail + 1
        end if
    end subroutine check_decay

    subroutine check_nonlinear_scalar(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: status
        real(dp) :: y0(1), err
        real(dp), allocatable :: t_out(:), y_out(:, :)

        y0(1) = 1.0_dp
        call ode_solve_radau(rhs_nl, 0.0_dp, 4.0_dp, y0, t_out, y_out, status, &
                             rtol=1.0e-13_dp)
        err = abs(y_out(1, size(t_out)) - 1.0_dp/5.0_dp)
        write (*, "(a,es10.3)") "  nonlinear scalar endpoint error ", err
        if (status%code /= FORTNUM_OK .or. .not. (err < 1.0e-12_dp)) then
            write (error_unit, "(a)") "  nonlinear scalar: endpoint inaccurate"
            nfail = nfail + 1
        end if
    end subroutine check_nonlinear_scalar

    ! The claim IAS15 is famous for: energy does not drift secularly. Over many
    ! orbits the relative energy error must stay tiny.
    subroutine check_kepler_energy(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp) :: y0(4), e0, e1, drift
        real(dp), parameter :: ecc = 0.5_dp
        real(dp), parameter :: twopi = 8.0_dp*atan(1.0_dp)
        integer :: n

        y0(1) = 1.0_dp - ecc
        y0(2) = 0.0_dp
        y0(3) = 0.0_dp
        y0(4) = sqrt((1.0_dp + ecc)/(1.0_dp - ecc))
        e0 = kepler_energy(y0)

        problem%rhs => rhs_kepler
        problem%t0 = 0.0_dp
        problem%t1 = 200.0_dp*twopi     ! 200 orbits
        problem%y0 = y0
        problem%rtol = 1.0e-12_dp
        problem%max_steps = 400000

        call ode_integrate_radau(problem, solution, status)
        n = size(solution%t)
        e1 = kepler_energy(solution%y(:, n))
        drift = abs((e1 - e0)/e0)

        write (*, "(a,es10.3,a,i0,a,i0)") "  kepler 200 orbits: rel energy drift ", &
            drift, "   steps ", solution%nsteps, "  nfev ", solution%nfev
        if (status%code /= FORTNUM_OK .or. .not. (drift < 1.0e-12_dp)) then
            write (error_unit, "(a)") "  kepler: energy drifted more than expected"
            nfail = nfail + 1
        end if
    end subroutine check_kepler_energy

    ! At matched tolerance the 15th-order method should not be worse than
    ! DOP853 on a smooth problem.
    subroutine check_vs_dop853(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: s1, s2
        real(dp) :: y0(1), er, ed
        real(dp), allocatable :: t1r(:), yr(:, :), t1d(:), yd(:, :)

        y0(1) = 1.0_dp
        call ode_solve_radau(rhs_nl, 0.0_dp, 4.0_dp, y0, t1r, yr, s1, rtol=1.0e-10_dp)
        call ode_solve_dop(rhs_nl, 0.0_dp, 4.0_dp, y0, t1d, yd, s2, rtol=1.0e-10_dp)
        er = abs(yr(1, size(t1r)) - 0.2_dp)
        ed = abs(yd(1, size(t1d)) - 0.2_dp)
        write (*, "(a,es10.3,a,es10.3)") "  vs dop853 at rtol 1e-10: radau ", er, &
            "   dop853 ", ed
        if (.not. (er <= max(ed*10.0_dp, 1.0e-11_dp))) then
            write (error_unit, "(a)") "  radau markedly less accurate than dop853"
            nfail = nfail + 1
        end if
    end subroutine check_vs_dop853

    subroutine check_backward(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_status_t) :: status
        real(dp) :: y0(1), yend, tend
        real(dp), allocatable :: t_out(:), y_out(:, :)

        y0(1) = exp(-2.0_dp)
        call ode_solve_radau(rhs_decay, 2.0_dp, 0.0_dp, y0, t_out, y_out, status, &
                             rtol=1.0e-12_dp)
        tend = t_out(size(t_out))
        yend = y_out(1, size(t_out))
        write (*, "(a,es10.3)") "  backward endpoint error ", abs(yend - 1.0_dp)
        if (status%code /= FORTNUM_OK .or. abs(tend) > 1.0e-13_dp .or. &
            .not. (abs(yend - 1.0_dp) < 1.0e-10_dp)) then
            write (error_unit, "(a)") "  backward integration failed"
            nfail = nfail + 1
        end if
    end subroutine check_backward

    subroutine check_bad_input(nfail)
        integer, intent(inout) :: nfail
        type(ode_problem_t) :: problem
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        real(dp) :: y0(1)

        y0(1) = 1.0_dp
        problem%rhs => rhs_decay
        problem%t0 = 0.0_dp
        problem%t1 = 0.0_dp        ! zero span must be rejected
        problem%y0 = y0
        call ode_integrate_radau(problem, solution, status)
        if (status%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, "(a)") "  zero-span problem was not rejected"
            nfail = nfail + 1
        end if

        problem%t1 = 1.0_dp
        problem%rtol = -1.0_dp     ! negative tolerance must be rejected
        call ode_integrate_radau(problem, solution, status)
        if (status%code /= FORTNUM_DOMAIN_ERROR) then
            write (error_unit, "(a)") "  negative rtol was not rejected"
            nfail = nfail + 1
        end if
    end subroutine check_bad_input

end program test_fortnum_ode_gauss_radau
