!> First-order forward-mode automatic differentiation over intervals.
!>
!> An idual_t carries an interval enclosure v of a value and enclosures d(k)
!> of its partial derivatives with respect to n seeded variables. Evaluated
!> on a box of seed values, the derivative parts enclose the gradient at
!> every point of the box, which is what mean-value and Lohner enclosures
!> need. The gradient has fixed capacity idual_max_vars so no operation
!> allocates; entries beyond n are zero.
module fortnum_idual
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, operator(+), &
        operator(-), operator(*), operator(/), operator(**), sqrt, exp, log, &
        sin, cos, sqr
    implicit none
    private

    integer, parameter, public :: idual_max_vars = 8

    public :: idual_t, idual_var, idual_const
    public :: operator(+), operator(-), operator(*), operator(/), operator(**)
    public :: sqrt, exp, log, sin, cos, sqr

    type :: idual_t
        integer :: n = 0
        type(interval_t) :: v
        type(interval_t) :: d(idual_max_vars)
    end type idual_t

    interface operator(+)
        module procedure add_dd, add_di, add_id, add_dr, add_rd, pos_d
    end interface operator(+)

    interface operator(-)
        module procedure sub_dd, sub_di, sub_id, sub_dr, sub_rd, neg_d
    end interface operator(-)

    interface operator(*)
        module procedure mul_dd, mul_di, mul_id, mul_dr, mul_rd
    end interface operator(*)

    interface operator(/)
        module procedure div_dd, div_di, div_id, div_dr, div_rd
    end interface operator(/)

    interface operator(**)
        module procedure pow_dn
    end interface operator(**)

    interface sqrt
        module procedure sqrt_d
    end interface sqrt

    interface exp
        module procedure exp_d
    end interface exp

    interface log
        module procedure log_d
    end interface log

    interface sin
        module procedure sin_d
    end interface sin

    interface cos
        module procedure cos_d
    end interface cos

    interface sqr
        module procedure sqr_d
    end interface sqr

contains

    !> Independent variable k of n with value enclosure x.
    pure function idual_var(x, k, n) result(r)
        type(interval_t), intent(in) :: x
        integer, intent(in) :: k, n
        type(idual_t) :: r

        r%n = min(n, idual_max_vars)
        r%v = x
        if (k >= 1 .and. k <= r%n) r%d(k) = interval(1.0_dp)
    end function idual_var

    pure elemental function idual_const(x, n) result(r)
        type(interval_t), intent(in) :: x
        integer, intent(in) :: n
        type(idual_t) :: r

        r%n = min(n, idual_max_vars)
        r%v = x
    end function idual_const

    !> Result with value v and derivative f' * a%d (chain rule).
    pure function chain(a, v, fprime) result(r)
        type(idual_t), intent(in) :: a
        type(interval_t), intent(in) :: v, fprime
        type(idual_t) :: r

        r%n = a%n
        r%v = v
        r%d(1:a%n) = fprime*a%d(1:a%n)
    end function chain

    pure elemental function pos_d(a) result(r)
        type(idual_t), intent(in) :: a
        type(idual_t) :: r

        r = a
    end function pos_d

    pure elemental function neg_d(a) result(r)
        type(idual_t), intent(in) :: a
        type(idual_t) :: r

        r%n = a%n
        r%v = -a%v
        r%d(1:a%n) = -a%d(1:a%n)
    end function neg_d

    pure elemental function add_dd(a, b) result(r)
        type(idual_t), intent(in) :: a, b
        type(idual_t) :: r

        r%n = max(a%n, b%n)
        r%v = a%v + b%v
        r%d(1:r%n) = a%d(1:r%n) + b%d(1:r%n)
    end function add_dd

    pure elemental function add_di(a, b) result(r)
        type(idual_t), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(idual_t) :: r

        r = a
        r%v = a%v + b
    end function add_di

    pure elemental function add_id(a, b) result(r)
        type(interval_t), intent(in) :: a
        type(idual_t), intent(in) :: b
        type(idual_t) :: r

        r = add_di(b, a)
    end function add_id

    pure elemental function add_dr(a, b) result(r)
        type(idual_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(idual_t) :: r

        r = add_di(a, interval(b))
    end function add_dr

    pure elemental function add_rd(a, b) result(r)
        real(dp), intent(in) :: a
        type(idual_t), intent(in) :: b
        type(idual_t) :: r

        r = add_di(b, interval(a))
    end function add_rd

    pure elemental function sub_dd(a, b) result(r)
        type(idual_t), intent(in) :: a, b
        type(idual_t) :: r

        r%n = max(a%n, b%n)
        r%v = a%v - b%v
        r%d(1:r%n) = a%d(1:r%n) - b%d(1:r%n)
    end function sub_dd

    pure elemental function sub_di(a, b) result(r)
        type(idual_t), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(idual_t) :: r

        r = a
        r%v = a%v - b
    end function sub_di

    pure elemental function sub_id(a, b) result(r)
        type(interval_t), intent(in) :: a
        type(idual_t), intent(in) :: b
        type(idual_t) :: r

        r = add_di(neg_d(b), a)
    end function sub_id

    pure elemental function sub_dr(a, b) result(r)
        type(idual_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(idual_t) :: r

        r = sub_di(a, interval(b))
    end function sub_dr

    pure elemental function sub_rd(a, b) result(r)
        real(dp), intent(in) :: a
        type(idual_t), intent(in) :: b
        type(idual_t) :: r

        r = sub_id(interval(a), b)
    end function sub_rd

    pure elemental function mul_dd(a, b) result(r)
        type(idual_t), intent(in) :: a, b
        type(idual_t) :: r

        r%n = max(a%n, b%n)
        r%v = a%v*b%v
        r%d(1:r%n) = a%d(1:r%n)*b%v + a%v*b%d(1:r%n)
    end function mul_dd

    pure elemental function mul_di(a, b) result(r)
        type(idual_t), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(idual_t) :: r

        r = chain(a, a%v*b, b)
    end function mul_di

    pure elemental function mul_id(a, b) result(r)
        type(interval_t), intent(in) :: a
        type(idual_t), intent(in) :: b
        type(idual_t) :: r

        r = mul_di(b, a)
    end function mul_id

    pure elemental function mul_dr(a, b) result(r)
        type(idual_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(idual_t) :: r

        r = mul_di(a, interval(b))
    end function mul_dr

    pure elemental function mul_rd(a, b) result(r)
        real(dp), intent(in) :: a
        type(idual_t), intent(in) :: b
        type(idual_t) :: r

        r = mul_di(b, interval(a))
    end function mul_rd

    !> (a/b)' = (a' - (a/b) b')/b.
    pure elemental function div_dd(a, b) result(r)
        type(idual_t), intent(in) :: a, b
        type(idual_t) :: r

        r%n = max(a%n, b%n)
        r%v = a%v/b%v
        r%d(1:r%n) = (a%d(1:r%n) - r%v*b%d(1:r%n))/b%v
    end function div_dd

    pure elemental function div_di(a, b) result(r)
        type(idual_t), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(idual_t) :: r

        r%n = a%n
        r%v = a%v/b
        r%d(1:a%n) = a%d(1:a%n)/b
    end function div_di

    pure elemental function div_id(a, b) result(r)
        type(interval_t), intent(in) :: a
        type(idual_t), intent(in) :: b
        type(idual_t) :: r
        type(interval_t) :: q

        q = a/b%v
        r = chain(b, q, -(q/b%v))
    end function div_id

    pure elemental function div_dr(a, b) result(r)
        type(idual_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(idual_t) :: r

        r = div_di(a, interval(b))
    end function div_dr

    pure elemental function div_rd(a, b) result(r)
        real(dp), intent(in) :: a
        type(idual_t), intent(in) :: b
        type(idual_t) :: r

        r = div_id(interval(a), b)
    end function div_rd

    pure elemental function pow_dn(a, n) result(r)
        type(idual_t), intent(in) :: a
        integer, intent(in) :: n
        type(idual_t) :: r

        if (n == 0) then
            r = idual_const(interval(1.0_dp), a%n)
        else
            r = chain(a, a%v**n, interval(real(n, dp))*a%v**(n - 1))
        end if
    end function pow_dn

    pure elemental function sqr_d(a) result(r)
        type(idual_t), intent(in) :: a
        type(idual_t) :: r

        r = chain(a, sqr(a%v), interval(2.0_dp)*a%v)
    end function sqr_d

    pure elemental function sqrt_d(a) result(r)
        type(idual_t), intent(in) :: a
        type(idual_t) :: r
        type(interval_t) :: s

        s = sqrt(a%v)
        r = chain(a, s, interval(0.5_dp)/s)
    end function sqrt_d

    pure elemental function exp_d(a) result(r)
        type(idual_t), intent(in) :: a
        type(idual_t) :: r
        type(interval_t) :: e

        e = exp(a%v)
        r = chain(a, e, e)
    end function exp_d

    pure elemental function log_d(a) result(r)
        type(idual_t), intent(in) :: a
        type(idual_t) :: r

        r = chain(a, log(a%v), interval(1.0_dp)/a%v)
    end function log_d

    pure elemental function sin_d(a) result(r)
        type(idual_t), intent(in) :: a
        type(idual_t) :: r

        r = chain(a, sin(a%v), cos(a%v))
    end function sin_d

    pure elemental function cos_d(a) result(r)
        type(idual_t), intent(in) :: a
        type(idual_t) :: r

        r = chain(a, cos(a%v), -sin(a%v))
    end function cos_d
end module fortnum_idual
