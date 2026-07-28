module hyperg_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_special_hypergeometric_1f1, only: hyperg_1f1_a1, &
        hyperg_1f1_a1_jvp
    implicit none
    private

    complex(dp), parameter :: fixed_b = cmplx(2.5_dp, 0.0_dp, dp)

contains

    function fortnum_hyperg_kernel(x) result(value) &
            bind(c, name="fortnum_hyperg_kernel")
        real(c_double), value :: x
        real(c_double) :: value
        complex(dp) :: result
        type(fortnum_status_t) :: status

        call hyperg_1f1_a1(fixed_b, cmplx(x, 0.0_dp, dp), result, status)
        if (.not. status_ok(status)) error stop "hyperg kernel failed"
        value = real(result, dp)
    end function fortnum_hyperg_kernel

    function fortnum_hyperg_kernel_autodiff(x) result(value) &
            bind(c, name="fortnum_hyperg_kernel_autodiff")
        real(c_double), value :: x
        real(c_double) :: value

        value = real_hyperg_a1_b25(x)
    end function fortnum_hyperg_kernel_autodiff

    function fortnum_hyperg_kernel_hybrid(x) result(value) &
            bind(c, name="fortnum_hyperg_kernel_hybrid")
        real(c_double), value :: x
        real(c_double) :: value

        value = real_hyperg_a1_b25(x)
    end function fortnum_hyperg_kernel_hybrid

    function fortnum_hyperg_kernel_jvp(x, dx) result(derivative) &
            bind(c, name="fortnum_hyperg_kernel_jvp")
        real(c_double), value :: x, dx
        real(c_double) :: derivative
        real(dp) :: z(2), direction(2), product(2)

        z = [x, 0.0_dp]
        direction = [dx, 0.0_dp]
        call hyperg_1f1_a1_jvp(z, fixed_b, direction, product)
        derivative = product(1)
    end function fortnum_hyperg_kernel_jvp

    function real_hyperg_a1_b25(x) result(value)
        real(c_double), intent(in) :: x
        real(c_double) :: value

        if (x < 0.0_c_double) then
            value = exp(x)*real_hyperg_series(1.5_c_double, -x)
        else if (x <= 60.0_c_double) then
            value = real_hyperg_series(1.0_c_double, x)
        else
            value = exp(log_gamma(2.5_c_double) + x - 1.5_c_double*log(x))
        end if
    end function real_hyperg_a1_b25

    function real_hyperg_series(a, x) result(value)
        real(c_double), intent(in) :: a, x
        real(c_double) :: value, term
        integer :: k, small_run

        term = 1.0_c_double
        value = term
        small_run = 0
        do k = 0, 4999
            term = term*(a + real(k, c_double))*x &
                /((2.5_c_double + real(k, c_double))*real(k + 1, c_double))
            value = value + term
            if (abs(term) <= epsilon(1.0_c_double)*abs(value)) then
                small_run = small_run + 1
                if (small_run >= 2) return
            else
                small_run = 0
            end if
        end do
        error stop "real hyperg series did not converge"
    end function real_hyperg_series

end module hyperg_kernel
