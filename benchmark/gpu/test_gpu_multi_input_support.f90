module test_gpu_multi_input_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_gpu_multi_input_wrappers, only: &
        multi_input_p2_jvp_batch, multi_input_p2_vjp_batch, &
        multi_input_p4_jvp_batch, multi_input_p4_vjp_batch, &
        multi_input_p8_jvp_batch, multi_input_p8_vjp_batch, &
        multi_input_p8_jvp_active_first_batch, &
        multi_input_p8_vjp_active_first_batch, &
        multi_input_p16_jvp_batch, multi_input_p16_vjp_batch
    use fortnum_gpu_backend_selection, only: select_multi_input_gpu_backend
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
        call test_p8_layouts()
        call test_backend_selection()
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

    subroutine test_p8_layouts()
        integer, parameter :: n = 1024, p = 8
        real(dp), parameter :: tolerance = 3.0e-13_dp
        real(dp) :: x_batch(n, p), direction_batch(n, p)
        real(dp) :: x_active(p, n), direction_active(p, n)
        real(dp) :: cotangent(n), values_batch(n), values_active(n)
        real(dp) :: products_batch(n), products_active(n)
        real(dp) :: adjoints_batch(n, p), adjoints_active(p, n)
        integer :: i, j

        do j = 1, p
            do i = 1, n
                x_batch(i, j) = -0.4_dp + 0.8_dp* &
                    real(mod(17*i + 11*j, 257), dp)/256.0_dp
                direction_batch(i, j) = -0.7_dp + 1.4_dp* &
                    real(mod(13*i + 19*j, 251), dp)/250.0_dp
                x_active(j, i) = x_batch(i, j)
                direction_active(j, i) = direction_batch(i, j)
            end do
        end do
        do i = 1, n
            cotangent(i) = -0.8_dp + 1.6_dp* &
                real(mod(23*i, 263), dp)/262.0_dp
        end do

        !$acc data copyin(x_batch, direction_batch, x_active, &
        !$acc& direction_active, cotangent) create(values_batch, &
        !$acc& values_active, products_batch, products_active, &
        !$acc& adjoints_batch, adjoints_active)
        !$omp target data map(to: x_batch, direction_batch, x_active, &
        !$omp& direction_active, cotangent) map(alloc: values_batch, &
        !$omp& values_active, products_batch, products_active, &
        !$omp& adjoints_batch, adjoints_active)
        call multi_input_p8_jvp_batch( &
            n, x_batch, direction_batch, values_batch, products_batch)
        call multi_input_p8_jvp_active_first_batch( &
            n, x_active, direction_active, values_active, products_active)
        call multi_input_p8_vjp_batch( &
            n, x_batch, cotangent, values_batch, adjoints_batch)
        call multi_input_p8_vjp_active_first_batch( &
            n, x_active, cotangent, values_active, adjoints_active)
        !$omp target update from(values_batch, values_active, &
        !$omp& products_batch, products_active, adjoints_batch, &
        !$omp& adjoints_active)
        !$omp end target data
        !$acc update self(values_batch, values_active, products_batch, &
        !$acc& products_active, adjoints_batch, adjoints_active)
        !$acc end data

        if (maxval(abs(values_batch - values_active)) > tolerance) then
            error stop "multi-input layouts disagree on value"
        end if
        if (maxval(abs(products_batch - products_active)) > tolerance) then
            error stop "multi-input layouts disagree on JVP"
        end if
        if (maxval(abs(adjoints_batch - transpose(adjoints_active))) > &
            tolerance) then
            error stop "multi-input layouts disagree on VJP"
        end if
    end subroutine test_p8_layouts

    subroutine test_backend_selection()
        integer, parameter :: batch_sizes(3) = [256, 65536, 1048576]
        character(7) :: backend
        logical :: found
        integer :: i

        do i = 1, size(batch_sizes)
            call check_selection("jvp", batch_sizes(i), .true.)
            call check_selection("jvp", batch_sizes(i), .false.)
            call check_selection("vjp", batch_sizes(i), .true.)
            call check_selection("vjp", batch_sizes(i), .false.)
        end do
        call select_multi_input_gpu_backend( &
            "jvp", 1024, .true., backend, found)
        if (found .or. len_trim(backend) /= 0) then
            error stop "unmeasured GPU workload was selected"
        end if
    end subroutine test_backend_selection

    subroutine check_selection(product, batch_size, is_resident)
        character(*), intent(in) :: product
        integer, intent(in) :: batch_size
        logical, intent(in) :: is_resident
        character(7) :: backend
        logical :: found

        call select_multi_input_gpu_backend( &
            product, batch_size, is_resident, backend, found)
        if (.not. found .or. trim(backend) /= "openacc") then
            error stop "measured GPU workload selection disagrees with evidence"
        end if
    end subroutine check_selection

end module test_gpu_multi_input_support
