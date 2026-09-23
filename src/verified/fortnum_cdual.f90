!> First-order forward-mode automatic differentiation over rectangular
!> complex intervals.
!>
!> A cdual_t carries a complex-box enclosure v of a value and enclosures
!> d(k) of its derivatives with respect to n seeded complex variables. The
!> derivative is the holomorphic derivative df/dz: every chain rule here is
!> the ordinary complex chain rule, which is why the same combinators as
!> fortnum_idual (interval, real forward AD) apply verbatim with cinterval_t
!> in place of interval_t. The gradient has fixed capacity cdual_max_vars so
!> no operation allocates; entries beyond n are zero.
!>
!> Consistency with Cauchy-Riemann: for an analytic f, the enclosure d(k)
!> also contains df/dx and (1/i) df/dy along the k-th seeded direction,
!> since a holomorphic function's partial derivatives satisfy
!> df/dx = f'(z) and df/dy = i f'(z). test_fortnum_cdual checks both partials
!> by finite differences against the same enclosure.
module fortnum_cdual
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: cinterval_t, cinterval, &
        operator(+), operator(-), operator(*), operator(/), sqrt, exp, &
        sin, cos
    implicit none
    private

    integer, parameter, public :: cdual_max_vars = 8

    public :: cdual_t, cdual_var, cdual_const
    public :: operator(+), operator(-), operator(*), operator(/)
    public :: sqrt, exp, sin, cos, cdual_sincos

    type :: cdual_t
        integer :: n = 0
        type(cinterval_t) :: v
        type(cinterval_t) :: d(cdual_max_vars)
    end type cdual_t

    interface operator(+)
        module procedure add_dd, add_dc, add_cd, add_dr, add_rd, pos_d
    end interface operator(+)

    interface operator(-)
        module procedure sub_dd, sub_dc, sub_cd, sub_dr, sub_rd, neg_d
    end interface operator(-)

    interface operator(*)
        module procedure mul_dd, mul_dc, mul_cd, mul_dr, mul_rd
    end interface operator(*)

    interface operator(/)
        module procedure div_dd, div_dc, div_cd, div_dr
    end interface operator(/)

    interface sqrt
        module procedure sqrt_d
    end interface sqrt

    interface exp
        module procedure exp_d
    end interface exp

    interface sin
        module procedure sin_d
    end interface sin

    interface cos
        module procedure cos_d
    end interface cos

contains

    !> Independent complex variable k of n with value enclosure x.
    pure function cdual_var(x, k, n) result(r)
        type(cinterval_t), intent(in) :: x
        integer, intent(in) :: k, n
        type(cdual_t) :: r

        r%n = min(n, cdual_max_vars)
        r%v = x
        if (k >= 1 .and. k <= r%n) r%d(k) = cinterval(1.0_dp, 0.0_dp)
    end function cdual_var

    pure elemental function cdual_const(x, n) result(r)
        type(cinterval_t), intent(in) :: x
        integer, intent(in) :: n
        type(cdual_t) :: r

        r%n = min(n, cdual_max_vars)
        r%v = x
    end function cdual_const

    !> Result with value v and derivative fprime * a%d (chain rule).
    pure function chain(a, v, fprime) result(r)
        type(cdual_t), intent(in) :: a
        type(cinterval_t), intent(in) :: v, fprime
        type(cdual_t) :: r

        r%n = a%n
        r%v = v
        r%d(1:a%n) = fprime*a%d(1:a%n)
    end function chain

    pure elemental function pos_d(a) result(r)
        type(cdual_t), intent(in) :: a
        type(cdual_t) :: r

        r = a
    end function pos_d

    pure elemental function neg_d(a) result(r)
        type(cdual_t), intent(in) :: a
        type(cdual_t) :: r

        r%n = a%n
        r%v = -a%v
        r%d(1:a%n) = -a%d(1:a%n)
    end function neg_d

    pure elemental function add_dd(a, b) result(r)
        type(cdual_t), intent(in) :: a, b
        type(cdual_t) :: r

        r%n = max(a%n, b%n)
        r%v = a%v + b%v
        r%d(1:r%n) = a%d(1:r%n) + b%d(1:r%n)
    end function add_dd

    pure elemental function add_dc(a, b) result(r)
        type(cdual_t), intent(in) :: a
        type(cinterval_t), intent(in) :: b
        type(cdual_t) :: r

        r = a
        r%v = a%v + b
    end function add_dc

    pure elemental function add_cd(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        type(cdual_t), intent(in) :: b
        type(cdual_t) :: r

        r = add_dc(b, a)
    end function add_cd

    pure elemental function add_dr(a, b) result(r)
        type(cdual_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(cdual_t) :: r

        r = add_dc(a, cinterval(b, 0.0_dp))
    end function add_dr

    pure elemental function add_rd(a, b) result(r)
        real(dp), intent(in) :: a
        type(cdual_t), intent(in) :: b
        type(cdual_t) :: r

        r = add_dc(b, cinterval(a, 0.0_dp))
    end function add_rd

    pure elemental function sub_dd(a, b) result(r)
        type(cdual_t), intent(in) :: a, b
        type(cdual_t) :: r

        r%n = max(a%n, b%n)
        r%v = a%v - b%v
        r%d(1:r%n) = a%d(1:r%n) - b%d(1:r%n)
    end function sub_dd

    pure elemental function sub_dc(a, b) result(r)
        type(cdual_t), intent(in) :: a
        type(cinterval_t), intent(in) :: b
        type(cdual_t) :: r

        r = a
        r%v = a%v - b
    end function sub_dc

    pure elemental function sub_cd(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        type(cdual_t), intent(in) :: b
        type(cdual_t) :: r

        r = add_dc(neg_d(b), a)
    end function sub_cd

    pure elemental function sub_dr(a, b) result(r)
        type(cdual_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(cdual_t) :: r

        r = sub_dc(a, cinterval(b, 0.0_dp))
    end function sub_dr

    pure elemental function sub_rd(a, b) result(r)
        real(dp), intent(in) :: a
        type(cdual_t), intent(in) :: b
        type(cdual_t) :: r

        r = sub_cd(cinterval(a, 0.0_dp), b)
    end function sub_rd

    pure elemental function mul_dd(a, b) result(r)
        type(cdual_t), intent(in) :: a, b
        type(cdual_t) :: r

        r%n = max(a%n, b%n)
        r%v = a%v*b%v
        r%d(1:r%n) = a%d(1:r%n)*b%v + a%v*b%d(1:r%n)
    end function mul_dd

    pure elemental function mul_dc(a, b) result(r)
        type(cdual_t), intent(in) :: a
        type(cinterval_t), intent(in) :: b
        type(cdual_t) :: r

        r = chain(a, a%v*b, b)
    end function mul_dc

    pure elemental function mul_cd(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        type(cdual_t), intent(in) :: b
        type(cdual_t) :: r

        r = mul_dc(b, a)
    end function mul_cd

    pure elemental function mul_dr(a, b) result(r)
        type(cdual_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(cdual_t) :: r

        r = mul_dc(a, cinterval(b, 0.0_dp))
    end function mul_dr

    pure elemental function mul_rd(a, b) result(r)
        real(dp), intent(in) :: a
        type(cdual_t), intent(in) :: b
        type(cdual_t) :: r

        r = mul_dc(b, cinterval(a, 0.0_dp))
    end function mul_rd

    !> (a/b)' = (a' - (a/b) b')/b.
    pure elemental function div_dd(a, b) result(r)
        type(cdual_t), intent(in) :: a, b
        type(cdual_t) :: r

        r%n = max(a%n, b%n)
        r%v = a%v/b%v
        r%d(1:r%n) = (a%d(1:r%n) - r%v*b%d(1:r%n))/b%v
    end function div_dd

    pure elemental function div_dc(a, b) result(r)
        type(cdual_t), intent(in) :: a
        type(cinterval_t), intent(in) :: b
        type(cdual_t) :: r

        r%n = a%n
        r%v = a%v/b
        r%d(1:a%n) = a%d(1:a%n)/b
    end function div_dc

    pure elemental function div_cd(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        type(cdual_t), intent(in) :: b
        type(cdual_t) :: r
        type(cinterval_t) :: q

        q = a/b%v
        r = chain(b, q, -(q/b%v))
    end function div_cd

    pure elemental function div_dr(a, b) result(r)
        type(cdual_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(cdual_t) :: r

        r = div_dc(a, cinterval(b, 0.0_dp))
    end function div_dr

    pure elemental function sqrt_d(a) result(r)
        type(cdual_t), intent(in) :: a
        type(cdual_t) :: r
        type(cinterval_t) :: s

        s = sqrt(a%v)
        r = chain(a, s, cinterval(0.5_dp, 0.0_dp)/s)
    end function sqrt_d

    pure elemental function exp_d(a) result(r)
        type(cdual_t), intent(in) :: a
        type(cdual_t) :: r
        type(cinterval_t) :: e

        e = exp(a%v)
        r = chain(a, e, e)
    end function exp_d

    pure elemental function sin_d(a) result(r)
        type(cdual_t), intent(in) :: a
        type(cdual_t) :: r

        r = chain(a, sin(a%v), cos(a%v))
    end function sin_d

    pure elemental function cos_d(a) result(r)
        type(cdual_t), intent(in) :: a
        type(cdual_t) :: r

        r = chain(a, cos(a%v), -sin(a%v))
    end function cos_d

    !> Simultaneous sin and cos, sharing the one evaluation of a%v.
    pure elemental subroutine cdual_sincos(a, s, c)
        type(cdual_t), intent(in) :: a
        type(cdual_t), intent(out) :: s, c
        type(cinterval_t) :: sv, cv

        sv = sin(a%v)
        cv = cos(a%v)
        s = chain(a, sv, cv)
        c = chain(a, cv, -sv)
    end subroutine cdual_sincos

end module fortnum_cdual
