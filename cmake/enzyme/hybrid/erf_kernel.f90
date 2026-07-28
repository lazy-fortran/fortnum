module erf_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_generated_erf_jvp, only: fortnum_erf_jvp
    implicit none
    private

contains

    function fortnum_erf_kernel(x) result(value) &
            bind(c, name="fortnum_erf_kernel")
        real(c_double), value :: x
        real(c_double) :: value

        value = erf(x)
    end function fortnum_erf_kernel

    function fortnum_erf_kernel_autodiff(x) result(value) &
            bind(c, name="fortnum_erf_kernel_autodiff")
        real(c_double), value :: x
        real(c_double) :: value

        value = erf(x)
    end function fortnum_erf_kernel_autodiff

    function fortnum_erf_kernel_jvp(x, dx) result(derivative) &
            bind(c, name="fortnum_erf_kernel_jvp")
        real(c_double), value :: x, dx
        real(c_double) :: derivative
        real(c_double) :: inputs(1), direction(1), product(1)

        inputs(1) = x
        direction(1) = dx
        call fortnum_erf_jvp(inputs, direction, product)
        derivative = product(1)
    end function fortnum_erf_kernel_jvp

end module erf_kernel
