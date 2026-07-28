module gamma_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_special_gamma, only: gamma_reg_p, gamma_reg_p_jvp
    implicit none
    private
    real(c_double), parameter :: fixed_a = 2.5_c_double

contains

    function fortnum_gamma_reg_p_kernel(x) result(value) &
            bind(c, name="fortnum_gamma_reg_p_kernel")
        real(c_double), value :: x
        real(c_double) :: value

        value = gamma_reg_p(fixed_a, x)
    end function fortnum_gamma_reg_p_kernel

    function fortnum_gamma_reg_p_kernel_autodiff(x) result(value) &
            bind(c, name="fortnum_gamma_reg_p_kernel_autodiff")
        real(c_double), value :: x
        real(c_double) :: value

        value = gamma_reg_p(fixed_a, x)
    end function fortnum_gamma_reg_p_kernel_autodiff

    function fortnum_gamma_reg_p_kernel_jvp(x, dx) result(derivative) &
            bind(c, name="fortnum_gamma_reg_p_kernel_jvp")
        real(c_double), value :: x, dx
        real(c_double) :: derivative
        real(c_double) :: inputs(2), direction(2), product(1)

        inputs = [x, fixed_a]
        direction = [dx, 0.0_c_double]
        call gamma_reg_p_jvp(inputs, direction, product)
        derivative = product(1)
    end function fortnum_gamma_reg_p_kernel_jvp

end module gamma_kernel
