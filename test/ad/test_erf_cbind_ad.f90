program test_erf_cbind_ad
    ! Derivative tests for fortnum_special_erf_cbind analytical products.
    !
    !   d/dx erf(x)  =  2/sqrt(pi) exp(-x^2)
    !   d/dx erfc(x) = -2/sqrt(pi) exp(-x^2)
    !
    ! Each analytic JVP is checked against central finite differences, and the
    ! grad (VJP) is checked against the JVP via the dot-product identity, using
    ! the fortnum_ad_test_utils harness.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64, error_unit
    use fortnum_ad_test_utils, only: check_jvp_vs_fd, dot_product_identity
    use fortnum_special_erf_cbind, only: fortnum_erf, fortnum_erfc, &
        fortnum_erf_jvp, fortnum_erfc_jvp, fortnum_erf_grad, fortnum_erfc_grad
    implicit none

    real(dp), parameter :: tol_fd  = 1.0e-7_dp ! central-FD tolerance (h ~ eps^1/3)
    real(dp), parameter :: tol_adj = 1.0e-13_dp ! adjoint identity tolerance

    abstract interface
        subroutine product_fn(x, seed, product)
            import :: dp
            real(dp), intent(in) :: x(:), seed(:)
            real(dp), intent(out) :: product(:)
        end subroutine product_fn
    end interface
    procedure(product_fn) :: manual_erf_product

    character(32) :: action, candidate, product
    integer(int64) :: iterations
    integer :: nfail

    call get_environment_variable("FORTNUM_ERF_ACTION", action)
    call get_environment_variable("FORTNUM_ERF_CANDIDATE", candidate)
    call get_environment_variable("FORTNUM_ERF_PRODUCT", product)
    call read_iterations(iterations)
    if (trim(action) == "--benchmark") then
        call benchmark_product(trim(candidate), trim(product), iterations)
        stop
    end if

    nfail = 0

    call test_erf_jvp(nfail)
    call test_erfc_jvp(nfail)
    call test_grad_adjoint(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    subroutine read_iterations(iterations)
        integer(int64), intent(out) :: iterations
        character(32) :: text
        integer :: status

        iterations = 20000000_int64
        call get_environment_variable("FORTNUM_ERF_ITERATIONS", text, &
            status=status)
        if (status == 0 .and. len_trim(text) > 0) read (text, *) iterations
    end subroutine read_iterations

    subroutine benchmark_product(candidate, product, iterations)
        character(*), intent(in) :: candidate, product
        integer(int64), intent(in) :: iterations
        real(dp) :: x(1), seed(1), result(1), sink, elapsed_ns
        integer(int64) :: first, last, rate, iteration
        procedure(product_fn), pointer :: selected

        if (candidate /= "generated" .and. candidate /= "diagnostic") then
            error stop "candidate must be generated or diagnostic"
        end if
        if (product /= "jvp" .and. product /= "vjp") then
            error stop "product must be jvp or vjp"
        end if
        if (candidate == "generated") then
            if (product == "jvp") then
                selected => fortnum_erf_jvp
            else
                selected => fortnum_erf_grad
            end if
        else
            selected => manual_erf_product
        end if
        seed(1) = 0.375_dp
        sink = 0.0_dp
        call system_clock(first, rate)
        do iteration = 1, iterations
            x(1) = 0.25_dp + 1.0e-6_dp*real(iand(iteration, 1023_int64), dp)
            call selected(x, seed, result)
            sink = sink + result(1)
        end do
        call system_clock(last)
        elapsed_ns = 1.0e9_dp*real(last - first, dp)/ &
            (real(rate, dp)*real(iterations, dp))
        write (*, "(a,',',a,',',f12.5,',',es24.16e3)") &
            candidate, product, elapsed_ns, sink
    end subroutine benchmark_product

    subroutine test_erf_jvp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: xpts(6), v(1), x(1)
        integer  :: i
        logical  :: ok

        xpts = [-3.0_dp, -0.7_dp, 0.0_dp, 0.5_dp, 1.5_dp, 4.0_dp]
        v    = [1.0_dp]
        do i = 1, size(xpts)
            x(1) = xpts(i)
            ok = check_jvp_vs_fd("fortnum_erf_jvp", f_erf, fortnum_erf_jvp, &
                x, v, tol_fd)
            if (.not. ok) nfail = nfail + 1
        end do
    end subroutine test_erf_jvp

    subroutine test_erfc_jvp(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: xpts(6), v(1), x(1)
        integer  :: i
        logical  :: ok

        xpts = [-3.0_dp, -0.7_dp, 0.0_dp, 0.5_dp, 1.5_dp, 4.0_dp]
        v    = [1.0_dp]
        do i = 1, size(xpts)
            x(1) = xpts(i)
            ok = check_jvp_vs_fd("fortnum_erfc_jvp", f_erfc, fortnum_erfc_jvp, &
                x, v, tol_fd)
            if (.not. ok) nfail = nfail + 1
        end do
    end subroutine test_erfc_jvp

    subroutine test_grad_adjoint(nfail)
        ! Scalar output: dot-product identity u*(Jv) = v*(J^T u).
        integer, intent(inout) :: nfail
        real(dp) :: x(1), u(1), v(1)
        x = [1.2_dp]
        u = [0.7_dp]
        v = [1.0_dp]
        if (.not. dot_product_identity("fortnum_erf_grad_adjoint", &
            fortnum_erf_jvp, fortnum_erf_grad, x, u, v, tol_adj)) &
            nfail = nfail + 1
        x = [-0.8_dp]
        u = [-1.3_dp]
        if (.not. dot_product_identity("fortnum_erfc_grad_adjoint", &
            fortnum_erfc_jvp, fortnum_erfc_grad, x, u, v, tol_adj)) &
            nfail = nfail + 1
    end subroutine test_grad_adjoint

    subroutine f_erf(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y(1) = fortnum_erf(x(1))
    end subroutine f_erf

    subroutine f_erfc(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        y(1) = fortnum_erfc(x(1))
    end subroutine f_erfc

end program test_erf_cbind_ad

subroutine manual_erf_product(x, seed, result)
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), parameter :: factor = 1.1283791670955126_dp
    real(dp), intent(in) :: x(:), seed(:)
    real(dp), intent(out) :: result(:)

    result(1) = seed(1)*factor*exp(-x(1)*x(1))
end subroutine manual_erf_product
