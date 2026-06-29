program test_dot_products
    ! Proves the mandatory adjoint identity u.(Jv)=v.(J^T u) on reference maps
    ! whose Jacobian transpose is known exactly: a rectangular linear map A x
    ! (J = A, J^T = A^T) and a nonlinear map whose JVP/VJP are hand-derived. A
    ! deliberately wrong VJP (A instead of A^T) must be caught.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_ad_test_utils, only: dot_product_identity, check_smoothness, &
        ad_status_t, AD_SMOOTH, AD_NONSMOOTH
    implicit none

    ! 3x2 linear map A: R^2 -> R^3. Module-scope so JVP/VJP share the matrix.
    real(dp), parameter :: A(3, 2) = reshape( &
        [1.0_dp, -2.0_dp, 0.5_dp, &
         3.0_dp,  4.0_dp, -1.0_dp], shape(A))

    integer :: nfail
    nfail = 0

    call test_linear_identity(nfail)
    call test_nonlinear_identity(nfail)
    call test_wrong_vjp_caught(nfail)
    call test_smoothness(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    ! Exact-to-rounding identity for the linear map; tol a few ulp.
    subroutine test_linear_identity(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x(2), u(3), v(2)
        x = [0.0_dp, 0.0_dp]            ! J constant; x irrelevant
        u = [1.0_dp, -0.5_dp, 2.0_dp]
        v = [0.75_dp, -1.25_dp]
        if (.not. dot_product_identity("linear", jvp_lin, vjp_lin, &
                x, u, v, tol=1.0e-13_dp)) nfail = nfail + 1
    end subroutine test_linear_identity

    ! Nonlinear map y = (x1^2, x1 x2, sin(x2)) : R^2 -> R^3.
    ! J = [[2 x1, 0],[x2, x1],[0, cos x2]]; JVP and VJP hand-coded below.
    subroutine test_nonlinear_identity(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x(2), u(3), v(2)
        x = [0.6_dp, -1.1_dp]
        u = [0.4_dp, 1.3_dp, -0.9_dp]
        v = [1.2_dp, 0.7_dp]
        if (.not. dot_product_identity("nonlinear", jvp_nl, vjp_nl, &
                x, u, v, tol=1.0e-13_dp)) nfail = nfail + 1
    end subroutine test_nonlinear_identity

    ! A VJP using A instead of A^T breaks the identity; the check must fail.
    ! Inverted assertion: the test passes when dot_product_identity returns
    ! .false. The expected failure is written to error_unit by the helper;
    ! that is intended diagnostic noise, not a test failure.
    subroutine test_wrong_vjp_caught(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x(2), u(3), v(2)
        logical  :: passed
        x = [0.0_dp, 0.0_dp]
        u = [1.0_dp, -0.5_dp, 2.0_dp]
        v = [0.75_dp, -1.25_dp]
        passed = dot_product_identity("wrong_vjp(expected-fail)", &
            jvp_lin, vjp_lin_wrong, x, u, v, tol=1.0e-13_dp)
        if (passed) then
            write (error_unit, '(a)') &
                "FAIL [wrong_vjp_caught] bad VJP slipped past the identity"
            nfail = nfail + 1
        end if
    end subroutine test_wrong_vjp_caught

    ! Same-trace case must report smooth; branch-change case must report
    ! non-smooth. Trace tags stand in for the branch a module recorded.
    subroutine test_smoothness(nfail)
        integer, intent(inout) :: nfail
        type(ad_status_t) :: s
        s = check_smoothness(trace_base=1, trace_pert=1, expect=AD_SMOOTH)
        if (.not. s%ok) then
            write (error_unit, '(a)') &
                "FAIL [smoothness] same trace not reported smooth"
            nfail = nfail + 1
        end if
        s = check_smoothness(trace_base=1, trace_pert=2, expect=AD_NONSMOOTH)
        if (.not. (s%ok .and. s%verdict == AD_NONSMOOTH)) then
            write (error_unit, '(a)') &
                "FAIL [smoothness] branch change not reported non-smooth"
            nfail = nfail + 1
        end if
    end subroutine test_smoothness

    ! ----------------------------------------------------- reference kernels

    subroutine jvp_lin(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        associate (unused_x => x); end associate
        jv = matmul(A, v)
    end subroutine jvp_lin

    subroutine vjp_lin(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        associate (unused_x => x); end associate
        jtu = matmul(transpose(A), u)
    end subroutine vjp_lin

    ! Wrong on purpose: drops the transpose, so J^T u is computed as A u-shaped
    ! garbage (shapes still conform only by coincidence of the 3x2 sizes, so
    ! use a hand assembly that conforms but is not the true adjoint).
    subroutine vjp_lin_wrong(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        associate (unused_x => x); end associate
        ! True adjoint is matmul(transpose(A), u); perturb one entry.
        jtu = matmul(transpose(A), u)
        jtu(1) = jtu(1) + 1.0_dp
    end subroutine vjp_lin_wrong

    subroutine jvp_nl(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        jv(1) = 2.0_dp*x(1)*v(1)
        jv(2) = x(2)*v(1) + x(1)*v(2)
        jv(3) = cos(x(2))*v(2)
    end subroutine jvp_nl

    subroutine vjp_nl(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        jtu(1) = 2.0_dp*x(1)*u(1) + x(2)*u(2)
        jtu(2) = x(1)*u(2) + cos(x(2))*u(3)
    end subroutine vjp_nl

end program test_dot_products
