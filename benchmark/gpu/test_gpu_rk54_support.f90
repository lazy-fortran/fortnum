module test_gpu_rk54_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_gpu_rk54_wrapper, only: rk54_dopri_linear_batch
    use fortnum_ode_rk54_device, only: RK54_ACCEPTED
    implicit none
    private

    public :: run_rk54_device_test

contains

    subroutine run_rk54_device_test()
        integer, parameter :: n = 8192
        real(dp), parameter :: h = 0.01_dp
        real(dp), parameter :: lambda(4) = [-1.0_dp, 0.5_dp, 1.0_dp, -2.0_dp]
        real(dp) :: y0(n, 4), yout(n, 4), expected(4), scale, error
        integer :: result(n), i, j

        do j = 1, 4
            do i = 1, n
                y0(i, j) = 0.5_dp + real(mod(17*i + 13*j, 997), dp)/997.0_dp
            end do
        end do
        !$acc data copyin(y0) create(yout, result)
        !$omp target data map(to: y0) map(alloc: yout, result)
        call rk54_dopri_linear_batch(n, y0, h, yout, result)
        !$omp target update from(yout, result)
        !$omp end target data
        !$acc update self(yout, result)
        !$acc end data

        if (any(result /= RK54_ACCEPTED)) then
            error stop "GPU DOPRI fixed step did not accept"
        end if
        do i = 1, n
            expected = y0(i, :)*exp(lambda*h)
            scale = max(1.0_dp, maxval(abs(expected)))
            error = maxval(abs(yout(i, :) - expected))
            if (error > 2.0e-13_dp*scale) then
                write (*, "(a,i0,a,es24.16,a,4es24.16)") &
                    "particle ", i, " max error ", error, " result ", yout(i, :)
                error stop "GPU DOPRI step disagrees with exponential oracle"
            end if
        end do
    end subroutine run_rk54_device_test

end module test_gpu_rk54_support
