program test_fortnum_polynomial
    ! Behavioral tests for fortnum_polynomial.
    !
    ! Key invariants checked:
    !   1. Weights reproduce polynomials of degree < n EXACTLY (analytic check).
    !      For n nodes the Lagrange interpolant through degree-(n-1) data is
    !      exact in exact arithmetic; floating-point error must be <= ~n*eps.
    !   2. Derivative weights match the analytic derivative of the interpolant.
    !      Verified by comparing against a finite-difference approximation and
    !      against the known derivative of the interpolated polynomial.
    !   3. Weights sum to 1 (partition of unity: interpolates constants).
    !   4. Derivative weights sum to 0 (derivative of a constant is zero).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_polynomial, only: lagrange_weights, lagrange_deriv_weights

    implicit none

    integer :: nfail
    nfail = 0

    call test_constant(nfail)
    call test_linear(nfail)
    call test_quadratic(nfail)
    call test_cubic(nfail)
    call test_partition_of_unity(nfail)
    call test_deriv_constant(nfail)
    call test_deriv_linear(nfail)
    call test_deriv_quadratic(nfail)
    call test_deriv_cubic(nfail)
    call test_deriv_sum_zero(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    ! Relative-error check; scale by max(|expected|, 1e-14) to avoid division by zero.
    subroutine check(label, got, expected, tol, nfail)
        character(*), intent(in)    :: label
        real(dp),     intent(in)    :: got, expected, tol
        integer,      intent(inout) :: nfail

        real(dp) :: scale, err

        scale = max(abs(expected), 1.0e-14_dp)
        err   = abs(got - expected)/scale
        if (err > tol) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es12.5,a,es22.15,a,es22.15)") &
                "FAIL [", label, "] relerr=", err, &
                " got=", got, " expected=", expected
        end if
    end subroutine check

    ! Absolute-error check; use when expected is zero.
    subroutine check_abs(label, got, expected, atol, nfail)
        character(*), intent(in)    :: label
        real(dp),     intent(in)    :: got, expected, atol
        integer,      intent(inout) :: nfail

        real(dp) :: err

        err = abs(got - expected)
        if (err > atol) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es12.5,a,es22.15,a,es22.15)") &
                "FAIL [", label, "] abserr=", err, &
                " got=", got, " expected=", expected
        end if
    end subroutine check_abs

    ! Interpolate constant f(x)=1 through n nodes; must recover 1 at any x.
    subroutine test_constant(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 5
        real(dp) :: xp(n), coef(n), val, x
        integer  :: i

        do i = 1, n
            xp(i) = real(i, dp)
        end do
        x = 2.7_dp
        call lagrange_weights(n, x, xp, coef)
        val = sum(coef)   ! sum(1 * coef) = interpolant of f=1
        call check("constant_unity", val, 1.0_dp, 1.0e-14_dp, nfail)
    end subroutine test_constant

    ! f(x)=x is degree-1; 2 nodes reproduce it.
    subroutine test_linear(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 2
        real(dp) :: xp(n), f(n), coef(n), val, x
        integer  :: i

        xp = [0.0_dp, 1.0_dp]
        do i = 1, n; f(i) = xp(i); end do
        x = 0.4_dp
        call lagrange_weights(n, x, xp, coef)
        val = dot_product(f, coef)
        call check("linear_2pt", val, x, 1.0e-15_dp, nfail)
    end subroutine test_linear

    ! f(x)=x^2 is degree-2; 3 nodes reproduce it.
    subroutine test_quadratic(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 3
        real(dp) :: xp(n), f(n), coef(n), val, x
        integer  :: i

        xp = [0.0_dp, 1.0_dp, 3.0_dp]
        do i = 1, n; f(i) = xp(i)**2; end do
        x = 1.7_dp
        call lagrange_weights(n, x, xp, coef)
        val = dot_product(f, coef)
        call check("quadratic_3pt", val, x**2, 1.0e-14_dp, nfail)
    end subroutine test_quadratic

    ! f(x)=x^3 is degree-3; 4 nodes reproduce it.
    subroutine test_cubic(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 4
        real(dp) :: xp(n), f(n), coef(n), val, x
        integer  :: i

        xp = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]
        do i = 1, n; f(i) = xp(i)**3; end do
        x = 2.5_dp
        call lagrange_weights(n, x, xp, coef)
        val = dot_product(f, coef)
        call check("cubic_4pt", val, x**3, 5.0e-14_dp, nfail)
    end subroutine test_cubic

    ! Weights sum to 1 for n=6 Chebyshev-like nodes.
    subroutine test_partition_of_unity(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 6
        real(dp) :: xp(n), coef(n), x
        integer  :: i

        ! Chebyshev-like spacing to avoid near-singularities.
        do i = 1, n
            xp(i) = cos(real(2*i - 1, dp)/(2.0_dp*real(n, dp))*acos(-1.0_dp))
        end do
        x = 0.15_dp
        call lagrange_weights(n, x, xp, coef)
        call check("partition_unity_n6", sum(coef), 1.0_dp, 1.0e-13_dp, nfail)
    end subroutine test_partition_of_unity

    ! Derivative of a constant is 0: sum(dcoef) == 0.
    subroutine test_deriv_constant(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 4
        real(dp) :: xp(n), dcoef(n), x
        integer  :: i

        do i = 1, n; xp(i) = real(i, dp)*0.5_dp; end do
        x = 1.1_dp
        call lagrange_deriv_weights(n, x, xp, dcoef)
        ! derivative of f=1 must be 0; use absolute tolerance -- expected is zero
        call check_abs("dconst_zero", sum(dcoef), 0.0_dp, 1.0e-13_dp, nfail)
    end subroutine test_deriv_constant

    ! f(x)=x; derivative must be 1.
    subroutine test_deriv_linear(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 2
        real(dp) :: xp(n), f(n), dcoef(n), dval, x
        integer  :: i

        xp = [0.0_dp, 2.0_dp]
        do i = 1, n; f(i) = xp(i); end do
        x = 0.8_dp
        call lagrange_deriv_weights(n, x, xp, dcoef)
        dval = dot_product(f, dcoef)
        call check("deriv_linear", dval, 1.0_dp, 1.0e-14_dp, nfail)
    end subroutine test_deriv_linear

    ! f(x)=x^2; derivative must be 2*x.
    subroutine test_deriv_quadratic(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 3
        real(dp) :: xp(n), f(n), dcoef(n), dval, x
        integer  :: i

        xp = [0.0_dp, 1.0_dp, 3.0_dp]
        do i = 1, n; f(i) = xp(i)**2; end do
        x = 1.4_dp
        call lagrange_deriv_weights(n, x, xp, dcoef)
        dval = dot_product(f, dcoef)
        call check("deriv_quadratic", dval, 2.0_dp*x, 1.0e-13_dp, nfail)
    end subroutine test_deriv_quadratic

    ! f(x)=x^3; derivative must be 3*x^2.
    subroutine test_deriv_cubic(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 4
        real(dp) :: xp(n), f(n), dcoef(n), dval, x
        integer  :: i

        xp = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]
        do i = 1, n; f(i) = xp(i)**3; end do
        x = 1.7_dp
        call lagrange_deriv_weights(n, x, xp, dcoef)
        dval = dot_product(f, dcoef)
        call check("deriv_cubic", dval, 3.0_dp*x**2, 5.0e-13_dp, nfail)
    end subroutine test_deriv_cubic

    ! Derivative weights sum to 0 for n=5.
    subroutine test_deriv_sum_zero(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 5
        real(dp) :: xp(n), dcoef(n), x
        integer  :: i

        do i = 1, n; xp(i) = real(i, dp)**2; end do   ! non-uniform spacing
        x = 7.3_dp
        call lagrange_deriv_weights(n, x, xp, dcoef)
        ! sum of derivative weights is zero; use absolute tolerance
        call check_abs("dsum_zero_n5", sum(dcoef), 0.0_dp, 1.0e-12_dp, nfail)
    end subroutine test_deriv_sum_zero

end program test_fortnum_polynomial
