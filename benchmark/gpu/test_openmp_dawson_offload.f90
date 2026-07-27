program test_openmp_dawson_offload
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use omp_lib, only: omp_is_initial_device
    use fortnum_gpu_batch_wrappers, only: dawson_value_jvp_batch, &
        dawson_vjp_batch
    implicit none

    integer, parameter :: n = 4096
    real(dp), parameter :: tolerance = 2.0e-14_dp
    real(dp), allocatable :: x(:), f(:), v(:), u(:)
    real(dp), allocatable :: values(:), products(:), adjoints(:)
    logical :: executed_on_device
    integer :: i

    allocate (x(n), f(n), v(n), u(n), values(n), products(n), adjoints(n))
    do i = 1, n
        x(i) = -3.0_dp + 6.0_dp*real(i - 1, dp)/real(n - 1, dp)
        f(i) = -0.75_dp + 1.5_dp*real(mod(37*i, n), dp)/real(n - 1, dp)
        v(i) = -1.0_dp + 2.0_dp*real(mod(53*i, n), dp)/real(n - 1, dp)
        u(i) = -0.8_dp + 1.6_dp*real(mod(71*i, n), dp)/real(n - 1, dp)
    end do

    executed_on_device = .false.
    !$omp target map(from: executed_on_device)
    executed_on_device = .not. omp_is_initial_device()
    !$omp end target
    if (.not. executed_on_device) then
        error stop "OpenMP target region executed on the initial device"
    end if

    call dawson_value_jvp_batch(n, x, f, v, values, products)
    call validate_outputs()
    call dawson_vjp_batch(n, x, f, u, adjoints)
    call validate_adjoint_identity()

    !$omp target data map(to: x, f, v, u) &
    !$omp& map(alloc: values, products, adjoints)
    call dawson_value_jvp_batch(n, x, f, v, values, products)
    call dawson_vjp_batch(n, x, f, u, adjoints)
    !$omp target update from(values, products, adjoints)
    !$omp end target data
    call validate_outputs()
    call validate_adjoint_identity()

contains

    subroutine validate_outputs()
        real(dp) :: expected_value, expected_product

        do i = 1, n
            expected_value = sin(f(i)) + f(i)*f(i)
            expected_product = v(i)*(1.0_dp - 2.0_dp*x(i)*f(i))* &
                (cos(f(i)) + 2.0_dp*f(i))
            if (abs(values(i) - expected_value) > tolerance) then
                error stop &
                    "OpenMP Dawson value disagrees with analytical oracle"
            end if
            if (abs(products(i) - expected_product) > tolerance) then
                error stop &
                    "OpenMP Dawson JVP disagrees with analytical oracle"
            end if
        end do
    end subroutine validate_outputs

    subroutine validate_adjoint_identity()
        real(dp) :: lhs, rhs, scale

        lhs = sum(u*products)
        rhs = sum(v*adjoints)
        scale = max(1.0_dp, abs(lhs), abs(rhs))
        if (abs(lhs - rhs) > 2.0e-13_dp*scale) then
            error stop "OpenMP Dawson batch violates adjoint identity"
        end if
    end subroutine validate_adjoint_identity

end program test_openmp_dawson_offload
