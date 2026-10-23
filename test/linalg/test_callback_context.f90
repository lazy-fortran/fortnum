module callback_context_cases
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_krylov, only: real_linear_operator_t, complex_linear_operator_t
    use fortnum_roots, only: scalar_root_function_t, scalar_root_derivative_t
    use fortnum_integrate_gk, only: gk_function_t
    implicit none
    real(dp), parameter :: matrix(3, 3) = reshape( &
        [3.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 4.0_dp, 0.0_dp, &
         0.0_dp, 1.0_dp, 5.0_dp], [3, 3])
    type, extends(real_linear_operator_t) :: real_map_t
        logical :: inverse = .false.
    contains
        procedure :: apply => real_apply
    end type
    type, extends(complex_linear_operator_t) :: complex_map_t
        logical :: inverse = .false.
    contains
        procedure :: apply => complex_apply
    end type
    type, extends(scalar_root_function_t) :: root_t
        real(dp) :: parameter = 2
    contains
        procedure :: evaluate => root_value
    end type
    type, extends(scalar_root_derivative_t) :: derivative_t
        real(dp) :: parameter = 2
    contains
        procedure :: evaluate => root_derivative
    end type
    type, extends(gk_function_t) :: polynomial_t
        real(dp) :: coefficient = 2
    contains
        procedure :: evaluate => polynomial_value
    end type
contains
    subroutine real_apply(self, input, output)
        class(real_map_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        if (self%inverse) then
            call real_inverse(input, output)
        else
            call real_product(input, output)
        end if
    end subroutine
    subroutine real_product(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        output = matmul(matrix, input)
    end subroutine
    subroutine real_inverse(input, output)
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        output = input/[3.0_dp, 4.0_dp, 5.0_dp]
    end subroutine
    subroutine complex_apply(self, input, output)
        class(complex_map_t), intent(inout) :: self
        complex(dp), intent(in) :: input(:)
        complex(dp), intent(out) :: output(:)
        if (self%inverse) then
            call complex_inverse(input, output)
        else
            call complex_product(input, output)
        end if
    end subroutine
    subroutine complex_product(input, output)
        complex(dp), intent(in) :: input(:)
        complex(dp), intent(out) :: output(:)
        output = matmul(cmplx(matrix, 0.1_dp*matrix, dp), input)
    end subroutine
    subroutine complex_inverse(input, output)
        complex(dp), intent(in) :: input(:)
        complex(dp), intent(out) :: output(:)
        output = input/cmplx([3.0_dp, 4.0_dp, 5.0_dp], &
                            [0.3_dp, 0.4_dp, 0.5_dp], dp)
    end subroutine
    pure function root_value(self, x) result(y)
        class(root_t), intent(in) :: self
        real(dp), intent(in) :: x
        real(dp) :: y
        y = x*x - self%parameter
    end function
    pure subroutine root_derivative(self, x, fx, dfx)
        class(derivative_t), intent(in) :: self
        real(dp), intent(in) :: x
        real(dp), intent(out) :: fx, dfx
        fx = x*x - self%parameter
        dfx = 2*x
    end subroutine
    pure function legacy_root(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = x*x - 2.0_dp
    end function
    pure subroutine legacy_derivative(x, fx, dfx)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: fx, dfx
        fx = x*x - 2.0_dp
        dfx = 2*x
    end subroutine
    function polynomial_value(self, x) result(fx)
        class(polynomial_t), intent(in) :: self
        real(dp), intent(in) :: x
        real(dp) :: fx
        fx = self%coefficient*x*x
    end function
    function legacy_polynomial(x) result(fx)
        real(dp), intent(in) :: x
        real(dp) :: fx
        fx = 2.0_dp*x*x
    end function
end module

program test_callback_context
    use callback_context_cases
    use, intrinsic :: iso_fortran_env, only: int64
    use fortnum_krylov, only: real_gmres, real_gmres_operator, real_fgmres, &
        real_fgmres_operator, complex_gmres, complex_gmres_operator, &
        complex_fgmres, complex_fgmres_operator
    use fortnum_roots, only: root_brent, root_bisect, root_newton
    use fortnum_integrate_gk, only: integrate_gk
    use fortnum_status, only: fortnum_status_t
    implicit none
    type(real_map_t) :: a, pre
    type(complex_map_t) :: za, zpre
    type(root_t) :: equation, second
    type(derivative_t) :: derivative
    type(polynomial_t) :: integrand
    type(fortnum_status_t) :: status, old_status
    real(dp) :: b(3), x(3), old(3), exact(3), r, old_r, root, old_root, error
    complex(dp) :: zb(3), zx(3), zold(3), zexact(3)
    integer :: info, old_info, iterations, old_iterations, mode
    exact = [0.25_dp, -0.5_dp, 1.5_dp]
    zexact = cmplx(exact, exact([3, 2, 1]), dp)
    pre%inverse = .true.; zpre%inverse = .true.
    call real_product(exact, b)
    call complex_product(zexact, zb)
    do mode = 1, 2
        x = 0; old = 0; zx = 0; zold = 0
        if (mode == 1) then
            call real_gmres(a, b, x, 1.0e-12_dp, 100, 3, info, iterations, r)
            call real_gmres_operator(real_product, b, old, 1.0e-12_dp, &
                100, 3, old_info, old_iterations, old_r)
        else
            call real_fgmres(a, pre, b, x, 1.0e-12_dp, 100, 3, info, iterations, r)
            call real_fgmres_operator(real_product, real_inverse, b, old, &
                1.0e-12_dp, 100, 3, old_info, old_iterations, old_r)
        end if
        if (info /= 0 .or. info /= old_info) error stop 'real status'
        if (iterations /= old_iterations) error stop 'real iterations'
        if (maxval(abs(x - exact)) > 1.0e-11_dp) error stop 'real oracle'
        if (any(transfer(x, [0_int64], 3) /= &
                transfer(old, [0_int64], 3))) error stop 'real adapter bits'
        if (transfer(r, 0_int64) /= transfer(old_r, 0_int64)) error stop 'residual bits'
        if (mode == 1) then
            call complex_gmres(za, zb, zx, 1.0e-12_dp, 100, 3, info, iterations, r)
            call complex_gmres_operator(complex_product, zb, zold, &
                1.0e-12_dp, 100, 3, old_info, old_iterations, old_r)
        else
            call complex_fgmres(za, zpre, zb, zx, 1.0e-12_dp, &
                100, 3, info, iterations, r)
            call complex_fgmres_operator(complex_product, complex_inverse, zb, &
                zold, 1.0e-12_dp, 100, 3, old_info, old_iterations, old_r)
        end if
        if (info /= 0 .or. info /= old_info) error stop 'complex status'
        if (iterations /= old_iterations) error stop 'complex iterations'
        if (maxval(abs(zx - zexact)) > 1.0e-11_dp) error stop 'complex oracle'
        if (any(transfer(zx, [0_int64], 6) /= &
                transfer(zold, [0_int64], 6))) error stop 'complex adapter bits'
    end do
    do mode = 1, 3
        select case (mode)
        case (1)
            call root_brent(equation, 0.0_dp, 2.0_dp, root, status)
            call root_brent(legacy_root, 0.0_dp, 2.0_dp, old_root, old_status)
        case (2)
            call root_bisect(equation, 0.0_dp, 2.0_dp, root, status)
            call root_bisect(legacy_root, 0.0_dp, 2.0_dp, old_root, old_status)
        case (3)
            call root_newton(derivative, 0.0_dp, 2.0_dp, 1.0_dp, root, status)
            call root_newton(legacy_derivative, 0.0_dp, 2.0_dp, 1.0_dp, &
                old_root, old_status)
        end select
        if (status%code /= 0 .or. old_status%code /= 0) error stop 'root status'
        if (abs(root - sqrt(2.0_dp)) > 2.0e-14_dp) error stop 'root oracle'
        if (transfer(root, 0_int64) /= transfer(old_root, 0_int64)) &
            error stop 'root adapter bits'
    end do
    second%parameter = 3
    call root_brent(second, 0.0_dp, 2.0_dp, root, status)
    if (abs(root - sqrt(3.0_dp)) > 2.0e-14_dp) error stop 'second context'
    call integrate_gk(integrand, 0.0_dp, 1.0_dp, 1.0e-13_dp, 1.0e-13_dp, &
        root, error, info)
    call integrate_gk(legacy_polynomial, 0.0_dp, 1.0_dp, &
        1.0e-13_dp, 1.0e-13_dp, old_root, old_r, old_info)
    if (info /= 0 .or. old_info /= 0) error stop 'quadrature status'
    if (abs(root - 2.0_dp/3) > 2.0e-14_dp) error stop 'quadrature oracle'
    if (transfer(root, 0_int64) /= transfer(old_root, 0_int64)) &
        error stop 'quadrature adapter bits'
    print *, 'context APIs: independent oracles and bitwise adapters passed'
end program
