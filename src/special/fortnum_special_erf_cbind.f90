module fortnum_special_erf_cbind
    ! Thin erf/erfc provider for the C ABI.
    !
    ! Fortran callers use the Fortran 2008 intrinsics erf/erfc directly; KAMEL
    ! C++ needs C-callable erf/erfc symbols. The
    ! module functions fortnum_erf/fortnum_erfc forward to the intrinsics with no
    ! reimplementation so the ABI layer can export them.
    !
    ! Analytical products are generated from the intrinsic expression.
    !   d/dx erf(x)  =  2/sqrt(pi) * exp(-x^2)            (DLMF 7.2.1)
    !   d/dx erfc(x) = -2/sqrt(pi) * exp(-x^2)
    ! Active argument: x (real scalar). Result is the active output.
    !   fortnum_erf_jvp / fortnum_erfc_jvp: forward product f'(x) * v.
    !   fortnum_erf_grad / fortnum_erfc_grad: scalar gradient (VJP).
    !
    ! Reference: scipy.special.erf / scipy.special.erfc; Fortran 2008 intrinsic.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_generated_erf_jvp, only: fortnum_erf_jvp
    use fortnum_generated_erf_vjp, only: fortnum_erf_grad => fortnum_erf_vjp
    use fortnum_generated_erfc_jvp, only: fortnum_erfc_jvp
    use fortnum_generated_erfc_vjp, only: fortnum_erfc_grad => fortnum_erfc_vjp
    implicit none
    private

    public :: fortnum_erf
    public :: fortnum_erfc
    public :: fortnum_erf_jvp
    public :: fortnum_erfc_jvp
    public :: fortnum_erf_grad
    public :: fortnum_erfc_grad

contains

    elemental function fortnum_erf(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = erf(x)
    end function fortnum_erf

    elemental function fortnum_erfc(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y
        y = erfc(x)
    end function fortnum_erfc

end module fortnum_special_erf_cbind
