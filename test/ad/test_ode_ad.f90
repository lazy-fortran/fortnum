module test_ode_ad_kernels
    ! Reference ODE problems and their analytic sensitivities for the #40 ode
    ! derivative sweep. Module-scope state is the test rig only (it lets the
    ! harness's plain f(x,y) / jvp(x,v,jv) callables reach a shared problem and
    ! frozen trace); fortnum_ode itself holds no global state.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
        ode_integrate, ode_integrate_jvp, ode_integrate_vjp
    implicit none
    private

    public :: set_linear, set_scalar
    public :: prim_y0, jvp_y0, vjp_y0
    public :: prim_scalar_k, jvp_scalar_k
    public :: A2, expm_Av, t1_lin

    ! 2x2 linear system y' = A y. Constant Jacobian J_f = A; the y0->y(t1) map
    ! is the matrix exponential exp(A t1), so the JVP w.r.t. y0 is exp(A t1) v.
    real(dp), parameter :: A2(2,2) = reshape( &
        [-0.3_dp, 0.7_dp, &
         -1.1_dp, 0.2_dp], shape(A2))
    real(dp), parameter :: t1_lin = 1.3_dp

    ! Scalar y' = -k y with parameter k. y(t1) = y0 exp(-k t1);
    ! dy/dk = -t1 y0 exp(-k t1). Sensitivity ODE: s' = -k s - y, s(0)=0.
    real(dp) :: kpar = 0.0_dp
    real(dp) :: y0_scalar = 0.0_dp
    real(dp) :: t1_scalar = 0.0_dp

contains

    ! ---- linear system y' = A y, sensitivity w.r.t. y0 ----

    subroutine rhs_lin(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt = matmul(A2, y)
    end subroutine rhs_lin

    ! Tangent RHS: J_f s = A s (no parameter part for y0-sensitivity).
    subroutine var_lin(t, y, s, dsdt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(in)  :: s(:)
        real(dp), intent(out) :: dsdt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t, unused_y => y); end associate
        dsdt = matmul(A2, s)
    end subroutine var_lin

    ! Adjoint tangent RHS: J_f^T lam = A^T lam.
    subroutine var_lin_adj(t, y, s, dsdt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(in)  :: s(:)
        real(dp), intent(out) :: dsdt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t, unused_y => y); end associate
        dsdt = matmul(transpose(A2), s)
    end subroutine var_lin_adj

    subroutine set_linear()
    end subroutine set_linear

    subroutine make_lin_problem(y0, problem)
        real(dp),            intent(in)  :: y0(:)
        type(ode_problem_t), intent(out) :: problem
        problem%rhs => rhs_lin
        problem%t0 = 0.0_dp
        problem%t1 = t1_lin
        problem%y0 = y0
        problem%rtol = 1.0e-10_dp
        problem%atol = 1.0e-12_dp
    end subroutine make_lin_problem

    ! Primal map y0 -> y(t1) for the linear system.
    subroutine prim_y0(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: ws
        type(ode_solution_t)   :: sol
        type(fortnum_status_t) :: st
        call make_lin_problem(x, problem)
        call ode_integrate(problem, ws, sol, st)
        y = sol%y(:, sol%nsteps + 1)
    end subroutine prim_y0

    ! Forward product J v = dy(t1)/dy0 . v via ode_integrate_jvp (seed s0 = v).
    subroutine jvp_y0(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: ws
        type(ode_solution_t)   :: sol
        type(fortnum_status_t) :: st
        real(dp), allocatable  :: s1(:)
        call make_lin_problem(x, problem)
        call ode_integrate(problem, ws, sol, st)
        call ode_integrate_jvp(problem, var_lin, v, sol, s1, st)
        jv = s1
    end subroutine jvp_y0

    ! Reverse product J^T u via ode_integrate_vjp.
    subroutine vjp_y0(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: ws
        type(ode_solution_t)   :: sol
        type(fortnum_status_t) :: st
        real(dp), allocatable  :: g(:)
        call make_lin_problem(x, problem)
        call ode_integrate(problem, ws, sol, st)
        call ode_integrate_vjp(problem, var_lin_adj, u, sol, g, st)
        jtu = g
    end subroutine vjp_y0

    ! exp(A t1) v by 2x2 closed-form matrix exponential, for the analytic JVP.
    ! exp(A t1) v in closed form. A2 has a complex-conjugate eigenvalue pair
    ! mu +/- i om, so the real matrix exponential is
    !   exp(A t) = e^{mu t} [ cos(om t) I + sin(om t)/om (A - mu I) ].
    function expm_Av(v) result(w)
        real(dp), intent(in) :: v(2)
        real(dp)             :: w(2)
        real(dp) :: tr, det, disc, mu, om, em, c0, c1
        real(dp) :: m(2,2), idn(2,2)
        tr = A2(1,1) + A2(2,2)
        det = A2(1,1)*A2(2,2) - A2(1,2)*A2(2,1)
        disc = tr*tr - 4.0_dp*det          ! < 0 for this A
        mu = 0.5_dp * tr
        om = 0.5_dp * sqrt(-disc)
        em = exp(mu * t1_lin)
        c0 = em * (cos(om*t1_lin) - mu*sin(om*t1_lin)/om)
        c1 = em * sin(om*t1_lin) / om
        idn = reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp], [2,2])
        m = c0*idn + c1*A2
        w = matmul(m, v)
    end function expm_Av

    ! ---- scalar y' = -k y, sensitivity w.r.t. parameter k ----

    subroutine set_scalar(k, y0, t1)
        real(dp), intent(in) :: k, y0, t1
        kpar = k
        y0_scalar = y0
        t1_scalar = t1
    end subroutine set_scalar

    subroutine rhs_scalar(t, y, dydt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dydt = -kpar * y
    end subroutine rhs_scalar

    ! Tangent RHS for s = dy/dk: s' = J_f s + df/dk = -k s - y.
    subroutine var_scalar_k(t, y, s, dsdt, ctx)
        real(dp), intent(in)  :: t
        real(dp), intent(in)  :: y(:)
        real(dp), intent(in)  :: s(:)
        real(dp), intent(out) :: dsdt(:)
        class(*), intent(in), optional :: ctx
        associate (unused_t => t); end associate
        dsdt = -kpar * s - y
    end subroutine var_scalar_k

    subroutine make_scalar_problem(problem)
        type(ode_problem_t), intent(out) :: problem
        problem%rhs => rhs_scalar
        problem%t0 = 0.0_dp
        problem%t1 = t1_scalar
        problem%y0 = [y0_scalar]
        problem%rtol = 1.0e-11_dp
        problem%atol = 1.0e-13_dp
    end subroutine make_scalar_problem

    ! Primal as a function of k: x(1) = k -> y(t1). Used for central FD.
    subroutine prim_scalar_k(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: ws
        type(ode_solution_t)   :: sol
        type(fortnum_status_t) :: st
        real(dp) :: ksave
        ksave = kpar
        kpar = x(1)
        call make_scalar_problem(problem)
        call ode_integrate(problem, ws, sol, st)
        y = [sol%y(1, sol%nsteps + 1)]
        kpar = ksave
    end subroutine prim_scalar_k

    ! Forward product dy(t1)/dk . v(1): seed s0 = 0 (y0 independent of k); the
    ! parameter direction enters through var_scalar_k's -y term, scaled by v(1).
    subroutine jvp_scalar_k(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        type(ode_problem_t)    :: problem
        type(ode_workspace_t)  :: ws
        type(ode_solution_t)   :: sol
        type(fortnum_status_t) :: st
        real(dp), allocatable  :: s1(:)
        real(dp) :: ksave
        ksave = kpar
        kpar = x(1)
        call make_scalar_problem(problem)
        call ode_integrate(problem, ws, sol, st)
        call ode_integrate_jvp(problem, var_scalar_k, [0.0_dp], sol, s1, st)
        jv = s1 * v(1)
        kpar = ksave
    end subroutine jvp_scalar_k

end module test_ode_ad_kernels


program test_ode_ad
    ! #40 ode derivative sweep: forward sensitivity (ode_integrate_jvp) and the
    ! discrete adjoint (ode_integrate_vjp), verified against
    !   (a) the analytic state-transition matrix exp(A t1) for y' = A y,
    !   (b) central finite difference of the primal solve, and
    !   (c) the mandatory adjoint dot-product identity u.(Jv) = v.(J^T u).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_ad_test_utils, only: check_jvp_vs_fd, dot_product_identity, &
        rel_err
    use test_ode_ad_kernels, only: set_linear, set_scalar, &
        prim_y0, jvp_y0, vjp_y0, prim_scalar_k, jvp_scalar_k, &
        A2, expm_Av, t1_lin
    implicit none

    integer :: nfail
    nfail = 0

    call test_linear_jvp_vs_analytic(nfail)
    call test_linear_jvp_vs_fd(nfail)
    call test_linear_dot_identity(nfail)
    call test_scalar_k_vs_fd(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    ! JVP for y' = A y must equal exp(A t1) v exactly (to integrator tol).
    subroutine test_linear_jvp_vs_analytic(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: y0(2), v(2), jv(2), want(2)
        real(dp) :: e
        integer  :: i
        call set_linear()
        y0 = [1.5_dp, -0.4_dp]
        v  = [0.8_dp, 1.2_dp]
        call jvp_y0(y0, v, jv)
        want = expm_Av(v)
        do i = 1, 2
            e = rel_err(jv(i), want(i))
            if (e > 1.0e-7_dp) then
                write (error_unit, '(a,i0,a,es24.16,a,es24.16,a,es12.4)') &
                    "FAIL linear jvp comp ", i, " got=", jv(i), &
                    " want=", want(i), " rel_err=", e
                nfail = nfail + 1
            end if
        end do
        write (*, '(a,2es14.6,a,2es14.6)') &
            "linear jvp     =", jv, "  exp(At1)v =", want
    end subroutine test_linear_jvp_vs_analytic

    ! Same JVP against central FD of the primal y0 -> y(t1) map.
    subroutine test_linear_jvp_vs_fd(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: y0(2), v(2)
        call set_linear()
        y0 = [1.5_dp, -0.4_dp]
        v  = [0.8_dp, 1.2_dp]
        if (.not. check_jvp_vs_fd("ode_linear_jvp_fd", prim_y0, jvp_y0, &
                y0, v, tol=1.0e-6_dp)) nfail = nfail + 1
    end subroutine test_linear_jvp_vs_fd

    ! Adjoint identity over the frozen trace: u.(Jv) = v.(J^T u).
    subroutine test_linear_dot_identity(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: y0(2), u(2), v(2)
        call set_linear()
        y0 = [1.5_dp, -0.4_dp]
        u  = [0.9_dp, -1.7_dp]
        v  = [0.8_dp, 1.2_dp]
        if (.not. dot_product_identity("ode_linear_adjoint", jvp_y0, vjp_y0, &
                y0, u, v, tol=1.0e-12_dp)) nfail = nfail + 1
    end subroutine test_linear_dot_identity

    ! Scalar y'=-k y: dy(t1)/dk = -t1 y0 exp(-k t1). JVP vs analytic and FD.
    subroutine test_scalar_k_vs_fd(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: k(1), v(1), jv(1), y0, t1, want, e
        y0 = 2.0_dp
        t1 = 0.9_dp
        k  = [1.4_dp]
        v  = [1.0_dp]
        call set_scalar(k(1), y0, t1)
        call jvp_scalar_k(k, v, jv)
        want = -t1 * y0 * exp(-k(1) * t1)
        e = rel_err(jv(1), want)
        if (e > 1.0e-7_dp) then
            write (error_unit, '(a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL scalar dk got=", jv(1), " want=", want, " rel_err=", e
            nfail = nfail + 1
        end if
        write (*, '(a,es16.8,a,es16.8)') &
            "scalar dy/dk   =", jv(1), "  analytic =", want
        if (.not. check_jvp_vs_fd("ode_scalar_k_jvp_fd", prim_scalar_k, &
                jvp_scalar_k, k, v, tol=1.0e-6_dp)) nfail = nfail + 1
    end subroutine test_scalar_k_vs_fd

end program test_ode_ad
