module test_gpu_multi_input_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_gpu_multi_input_wrappers, only: &
        multi_input_p2_jvp_batch, multi_input_p2_vjp_batch, &
        multi_input_p4_jvp_batch, multi_input_p4_vjp_batch, &
        multi_input_p8_jvp_batch, multi_input_p8_vjp_batch, &
        multi_input_p16_jvp_batch, multi_input_p16_vjp_batch
    implicit none
    private

    public :: run_multi_input_tests

contains

    subroutine run_multi_input_tests()
        integer, parameter :: active_sizes(4) = [2, 4, 8, 16]
        integer :: i

        do i = 1, size(active_sizes)
            call test_size(active_sizes(i))
        end do
    end subroutine run_multi_input_tests

    subroutine test_size(nactive)
        integer, intent(in) :: nactive
        integer, parameter :: n = 1024
        real(dp), parameter :: tolerance = 3.0e-13_dp
        real(dp), allocatable :: x(:, :), direction(:, :), cotangent(:)
        real(dp), allocatable :: values(:), products(:), adjoints(:, :)
        real(dp) :: expected_value, expected_product, gradient
        real(dp) :: lhs, rhs, scale, sum_x
        integer :: i, j

        allocate (x(n, nactive), direction(n, nactive), cotangent(n))
        allocate (values(n), products(n), adjoints(n, nactive))
        do j = 1, nactive
            do i = 1, n
                x(i, j) = -0.4_dp + 0.8_dp* &
                    real(mod(17*i + 11*j, 257), dp)/256.0_dp
                direction(i, j) = -0.7_dp + 1.4_dp* &
                    real(mod(13*i + 19*j, 251), dp)/250.0_dp
            end do
        end do
        do i = 1, n
            cotangent(i) = -0.8_dp + 1.6_dp* &
                real(mod(23*i, 263), dp)/262.0_dp
        end do

        !$acc data copyin(x, direction, cotangent) &
        !$acc& create(values, products, adjoints)
        !$omp target data map(to: x, direction, cotangent) &
        !$omp& map(alloc: values, products, adjoints)
        select case (nactive)
        case (2)
            call multi_input_p2_jvp_batch( &
                n, x, direction, values, products)
            call multi_input_p2_vjp_batch( &
                n, x, cotangent, values, adjoints)
        case (4)
            call multi_input_p4_jvp_batch( &
                n, x, direction, values, products)
            call multi_input_p4_vjp_batch( &
                n, x, cotangent, values, adjoints)
        case (8)
            call multi_input_p8_jvp_batch( &
                n, x, direction, values, products)
            call multi_input_p8_vjp_batch( &
                n, x, cotangent, values, adjoints)
        case (16)
            call multi_input_p16_jvp_batch( &
                n, x, direction, values, products)
            call multi_input_p16_vjp_batch( &
                n, x, cotangent, values, adjoints)
        end select
        !$omp target update from(values, products, adjoints)
        !$omp end target data
        !$acc update self(values, products, adjoints)
        !$acc end data

        do i = 1, n
            sum_x = sum(x(i, :))
            expected_value = sum(sin(x(i, :))) + 0.5_dp*sum_x*sum_x
            expected_product = 0.0_dp
            do j = 1, nactive
                gradient = cos(x(i, j)) + sum_x
                expected_product = expected_product + &
                    gradient*direction(i, j)
                if (abs(adjoints(i, j) - cotangent(i)*gradient) > &
                    tolerance) then
                    error stop "multi-input VJP disagrees with CPU oracle"
                end if
            end do
            if (abs(values(i) - expected_value) > tolerance) then
                error stop "multi-input value disagrees with CPU oracle"
            end if
            if (abs(products(i) - expected_product) > tolerance) then
                error stop "multi-input JVP disagrees with CPU oracle"
            end if
        end do

        lhs = sum(cotangent*products)
        rhs = sum(direction*adjoints)
        scale = max(1.0_dp, abs(lhs), abs(rhs))
        if (abs(lhs - rhs) > tolerance*scale) then
            error stop "multi-input products violate adjoint identity"
        end if
    end subroutine test_size

end module test_gpu_multi_input_support
