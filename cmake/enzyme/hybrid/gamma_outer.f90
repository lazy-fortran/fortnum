module gamma_outer_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    private

    public :: outer

    interface
        function fortnum_gamma_reg_p_kernel(x) result(value) &
                bind(c, name="fortnum_gamma_reg_p_kernel")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: value
        end function fortnum_gamma_reg_p_kernel

        function fortnum_gamma_reg_p_kernel_autodiff(x) result(value) &
                bind(c, name="fortnum_gamma_reg_p_kernel_autodiff")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: value
        end function fortnum_gamma_reg_p_kernel_autodiff
    end interface

contains

    function outer(x) result(value) bind(c, name="fortnum_gamma_outer")
        real(c_double), value :: x
        real(c_double) :: value, inner

        inner = fortnum_gamma_reg_p_kernel(x)
        value = sin(inner) + inner*inner
    end function outer

    function outer_autodiff(x) result(value) &
            bind(c, name="fortnum_gamma_outer_autodiff")
        real(c_double), value :: x
        real(c_double) :: value, inner

        inner = fortnum_gamma_reg_p_kernel_autodiff(x)
        value = sin(inner) + inner*inner
    end function outer_autodiff

end module gamma_outer_kernel
