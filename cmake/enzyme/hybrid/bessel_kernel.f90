module bessel_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_special_bessel, only: bessel_in
    implicit none
    private

contains

    function fortnum_bessel_i0_kernel(x) result(value) &
            bind(c, name="fortnum_bessel_i0_kernel")
        real(c_double), value :: x
        real(c_double) :: value

        value = bessel_in(0, x)
    end function fortnum_bessel_i0_kernel

    function fortnum_bessel_i1_kernel(x) result(value) &
            bind(c, name="fortnum_bessel_i1_kernel")
        real(c_double), value :: x
        real(c_double) :: value

        value = bessel_in(1, x)
    end function fortnum_bessel_i1_kernel

    function fortnum_bessel_i0_kernel_autodiff(x) result(value) &
            bind(c, name="fortnum_bessel_i0_kernel_autodiff")
        real(c_double), value :: x
        real(c_double) :: value

        value = bessel_in(0, x)
    end function fortnum_bessel_i0_kernel_autodiff

end module bessel_kernel
