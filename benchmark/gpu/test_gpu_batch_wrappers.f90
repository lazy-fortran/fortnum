program test_gpu_batch_wrappers
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_gpu_batch_wrappers, only: dawson_value_jvp_batch
    implicit none

    integer, parameter :: n = 17
    real(dp) :: x(n), f(n), v(n), values(n), products(n)
    real(dp) :: expected_value, expected_product
    integer :: i, failed

    failed = 0
    do i = 1, n
        x(i) = -1.2_dp + 0.15_dp*real(i - 1, dp)
        f(i) = 0.1_dp + 0.02_dp*real(i - 1, dp)
        v(i) = -0.4_dp + 0.05_dp*real(i - 1, dp)
    end do

    call dawson_value_jvp_batch(n, x, f, v, values, products)
    do i = 1, n
        expected_value = sin(f(i)) + f(i)**2
        expected_product = v(i)*(1.0_dp - 2.0_dp*x(i)*f(i))* &
            (cos(f(i)) + 2.0_dp*f(i))
        if (abs(values(i) - expected_value) > 2.0e-15_dp) failed = failed + 1
        if (abs(products(i) - expected_product) > 2.0e-15_dp) failed = failed + 1
    end do

    if (failed > 0) then
        print '(a,i0)', "gpu batch wrapper failures: ", failed
        stop 1
    end if
end program test_gpu_batch_wrappers
