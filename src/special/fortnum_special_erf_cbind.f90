module fortnum_special_erf_cbind
    ! Thin erf/erfc provider for the C ABI.
    !
    ! Fortran callers use the Fortran 2008 intrinsics erf/erfc directly; KAMEL
    ! C++ (the C++ caller) needs C-callable symbols. The
    ! module functions fortnum_erf/fortnum_erfc forward to the intrinsics with no
    ! reimplementation so the ABI layer can export them.
    !
    ! Derivative policy (ad.md sec 1, sec 4): transparent.
    !   d/dx erf(x)  =  2/sqrt(pi) * exp(-x^2)            (DLMF 7.2.1)
    !   d/dx erfc(x) = -2/sqrt(pi) * exp(-x^2)
    ! Active argument: x (real scalar). Result is the active output.
    !   fortnum_erf_jvp / fortnum_erfc_jvp: forward product f'(x) * v.
    !   fortnum_erf_grad / fortnum_erfc_grad: scalar gradient (VJP).
    !
    ! Reference: scipy.special.erf / scipy.special.erfc; Fortran 2008 intrinsic.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: fortnum_erf
    public :: fortnum_erfc
    public :: fortnum_erf_jvp
    public :: fortnum_erfc_jvp
    public :: fortnum_erf_grad
    public :: fortnum_erfc_grad

    ! 2/sqrt(pi), the common derivative prefactor.
    real(dp), parameter :: two_over_sqrt_pi = 1.1283791670955126_dp

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

    ! Forward product: jv = erf'(x) v, erf'(x) = 2/sqrt(pi) exp(-x^2).
    ! x and v are length-1 arrays conforming to the harness interface.
    subroutine fortnum_erf_jvp(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        jv(1) = two_over_sqrt_pi * exp(-x(1)*x(1)) * v(1)
    end subroutine fortnum_erf_jvp

    ! Forward product: jv = erfc'(x) v, erfc'(x) = -2/sqrt(pi) exp(-x^2).
    subroutine fortnum_erfc_jvp(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        jv(1) = -two_over_sqrt_pi * exp(-x(1)*x(1)) * v(1)
    end subroutine fortnum_erfc_jvp

    ! Scalar gradient (VJP): jtu = u * erf'(x). For scalar output the JVP and
    ! VJP coincide; u is the output cotangent.
    subroutine fortnum_erf_grad(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        jtu(1) = u(1) * two_over_sqrt_pi * exp(-x(1)*x(1))
    end subroutine fortnum_erf_grad

    ! Scalar gradient (VJP): jtu = u * erfc'(x).
    subroutine fortnum_erfc_grad(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        jtu(1) = -u(1) * two_over_sqrt_pi * exp(-x(1)*x(1))
    end subroutine fortnum_erfc_grad

end module fortnum_special_erf_cbind
