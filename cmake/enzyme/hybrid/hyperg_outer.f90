module hyperg_outer_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    implicit none
    private
    public :: outer

    interface
        function fortnum_hyperg_kernel_hybrid(x) result(value) &
                bind(c, name="fortnum_hyperg_kernel_hybrid")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: value
        end function fortnum_hyperg_kernel_hybrid

        function fortnum_hyperg_kernel_autodiff(x) result(value) &
                bind(c, name="fortnum_hyperg_kernel_autodiff")
            import :: c_double
            real(c_double), value :: x
            real(c_double) :: value
        end function fortnum_hyperg_kernel_autodiff
    end interface

contains

    function outer(x) result(value) bind(c, name="fortnum_hyperg_outer")
        real(c_double), value :: x
        real(c_double) :: value, inner

        inner = fortnum_hyperg_kernel_hybrid(x)
        value = log(1.0_c_double + inner*inner)
    end function outer

    function outer_autodiff(x) result(value) &
            bind(c, name="fortnum_hyperg_outer_autodiff")
        real(c_double), value :: x
        real(c_double) :: value, inner

        inner = fortnum_hyperg_kernel_autodiff(x)
        value = log(1.0_c_double + inner*inner)
    end function outer_autodiff

end module hyperg_outer_kernel
