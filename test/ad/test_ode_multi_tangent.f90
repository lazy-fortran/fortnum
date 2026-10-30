program test_ode_multi_tangent
    ! Oracles: central differences of the frozen Cash-Karp replay, and the
    ! single-direction ode_integrate_jvp column by column.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate, ode_integrate_jvp, ode_integrate_jvp_multi
    use fortnum_ode_cash_karp, only: cash_karp_step
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none

    integer, parameter :: n = 5, np = 2, m = 3
    real(dp) :: p(np), active(np), y0(n), s0(n, m), vp(np, m)
    real(dp) :: fd(n, m), single(n, m)
    real(dp) :: plus(n), minus(n), eps, scale
    real(dp), allocatable :: s1(:, :), column(:), bad(:, :)
    type(ode_problem_t) :: problem
    type(ode_workspace_t) :: workspace
    type(ode_solution_t) :: solution
    type(fortnum_status_t) :: status
    integer :: i, k, direction

    p = [0.3_dp, 0.05_dp]
    active = p
    y0 = [(0.4_dp + 0.1_dp*real(i, dp), i = 1, n)]
    do k = 1, m
        do i = 1, n
            s0(i, k) = sin(real(i + 2*k, dp))
        end do
        vp(:, k) = [cos(real(k, dp)), 0.5_dp*real(k - 2, dp)]
    end do

    problem%rhs => rhs
    problem%t0 = 0.0_dp
    problem%t1 = 2.0_dp
    problem%y0 = y0
    problem%rtol = 1.0e-8_dp
    problem%atol = 1.0e-10_dp
    call ode_integrate(problem, workspace, solution, status)
    call require(status_ok(status), "primal")
    call require(solution%nsteps > 5, "nontrivial schedule")

    call ode_integrate_jvp_multi(problem, tangents, s0, solution, s1, status)
    call require(status_ok(status), "multi-direction tangent")
    call require(all(shape(s1) == [n, m]), "tangent shape")

    eps = 1.0e-6_dp
    do k = 1, m
        call replay(y0 + eps*s0(:, k), p + eps*vp(:, k), plus)
        call replay(y0 - eps*s0(:, k), p - eps*vp(:, k), minus)
        fd(:, k) = (plus - minus)/(2*eps)
        direction = k
        call ode_integrate_jvp(problem, tangent, s0(:, k), solution, column, &
            status)
        call require(status_ok(status), "single-direction tangent")
        single(:, k) = column
    end do
    scale = maxval(abs(fd))
    print '(a,es10.3)', "max |multi - central difference| / scale = ", &
        maxval(abs(s1 - fd))/scale
    print '(a,es10.3)', "max |multi - single| / scale = ", &
        maxval(abs(s1 - single))/scale
    call require(maxval(abs(s1 - fd)) <= 1.0e-8_dp*scale, &
        "multi-direction tangent matches central differences")
    call require(maxval(abs(s1 - single)) <= 1.0e-14_dp*scale, &
        "multi-direction tangent matches single directions")

    allocate (bad(n + 1, m))
    bad = 0
    call ode_integrate_jvp_multi(problem, tangents, bad, solution, s1, status)
    call require(.not. status_ok(status), "seed row count is checked")
    print '(a)', "ode multi-direction tangent passed"

contains

    subroutine rhs(t, y, dydt, ctx)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        call field(t, y, active, dydt)
    end subroutine rhs

    subroutine field(t, y, q, dydt)
        real(dp), intent(in) :: t, y(:), q(:)
        real(dp), intent(out) :: dydt(:)
        integer :: j
        do j = 1, size(y)
            dydt(j) = -(q(1) + 0.01_dp*real(j, dp))*y(j) + 0.03_dp*y(j)**2 + &
                q(2)*y(previous(j))*cos(t)
        end do
    end subroutine field

    ! J_f s + f_p v for one column.
    subroutine directional(t, y, s, v, ds)
        real(dp), intent(in) :: t, y(:), s(:), v(:)
        real(dp), intent(out) :: ds(:)
        integer :: j
        do j = 1, size(y)
            ds(j) = (-(p(1) + 0.01_dp*real(j, dp)) + 0.06_dp*y(j))*s(j) + &
                p(2)*s(previous(j))*cos(t) - v(1)*y(j) + &
                v(2)*y(previous(j))*cos(t)
        end do
    end subroutine directional

    subroutine tangent(t, y, s, dsdt, ctx)
        real(dp), intent(in) :: t, y(:), s(:)
        real(dp), intent(out) :: dsdt(:)
        class(*), intent(in), optional :: ctx
        call directional(t, y, s, vp(:, direction), dsdt)
    end subroutine tangent

    subroutine tangents(t, y, s, dsdt, ctx)
        real(dp), intent(in) :: t, y(:), s(:, :)
        real(dp), intent(out) :: dsdt(:, :)
        class(*), intent(in), optional :: ctx
        integer :: j
        do j = 1, size(s, 2)
            call directional(t, y, s(:, j), vp(:, j), dsdt(:, j))
        end do
    end subroutine tangents

    integer function previous(j)
        integer, intent(in) :: j
        previous = j - 1
        if (previous == 0) previous = n
    end function previous

    ! The primal map with the accepted schedule held fixed.
    subroutine replay(start, q, final)
        real(dp), intent(in) :: start(:), q(:)
        real(dp), intent(out) :: final(:)
        real(dp), dimension(n) :: y, y5, yerr, ytmp, k1, k2, k3, k4, k5, k6
        integer :: step, evaluations
        y = start
        active = q
        evaluations = 0
        do step = 1, solution%nsteps
            call cash_karp_step(rhs, solution%t(step), y, solution%h(step), &
                .false., k1, k2, k3, k4, k5, k6, ytmp, y5, yerr, evaluations)
            y = y5
        end do
        final = y
        active = p
    end subroutine replay

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (condition) return
        print '(2a)', "FAIL: ", label
        error stop 1
    end subroutine require
end program test_ode_multi_tangent
