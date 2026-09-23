!> Real and rectangular complex interval arithmetic with outward rounding.
!>
!> Every operation returns an interval that contains the exact result for all
!> points of the arguments. Endpoints are computed in round-to-nearest and
!> moved outward with fortnum_rounding (successor formula, no rounding-mode
!> change). Only IEEE +, -, *, / and sqrt are used; exp, log, sin, cos, sinh,
!> cosh and pi are enclosed by Taylor polynomials with Lagrange remainders
!> after argument reduction against two-word enclosures of ln 2 and pi/2, so
!> no libm accuracy claim enters. log uses the libm value only as a starting
!> guess and verifies the enclosure through the rigorous exp.
!>
!> Conventions: intervals have lo <= hi. A divisor containing zero yields the
!> entire line [-inf, inf]. sqrt and log act on the part of the argument inside
!> their domain (sqrt on [0, inf), log on (0, inf)); callers that need a
!> domain check test the argument first.
module fortnum_interval
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, &
        ieee_is_nan
    use fortnum_rounding, only: round_up, round_down, sqrt_up, sqrt_down
    implicit none
    private

    public :: interval_t, cinterval_t, interval, cinterval
    public :: operator(+), operator(-), operator(*), operator(/), operator(**)
    public :: sqrt, exp, log, sin, cos, sinh, cosh, abs
    public :: sqr, interval_pi, hull, intersect, mid, rad, width, mag, mig
    public :: contains, contains_zero, subset, disjoint, entire, is_empty
    public :: real_part, imag_part, conjg, cabs_up, cabs_down
    ! Runtime interface called by fortsym-emitted rigorous kernels.
    public :: ipoint, ienclose, iadd, isub, imul, idiv, ineg, iinv, isqrt
    public :: ipowi, iscale

    type :: interval_t
        real(dp) :: lo = 0.0_dp
        real(dp) :: hi = 0.0_dp
    end type interval_t

    type :: cinterval_t
        type(interval_t) :: re
        type(interval_t) :: im
    end type cinterval_t

    interface interval
        module procedure iv_from_real, iv_from_bounds, iv_from_int
    end interface interval

    interface cinterval
        module procedure ci_from_complex, ci_from_parts, ci_from_real_parts
    end interface cinterval

    interface operator(+)
        module procedure add_ii, add_ir, add_ri, add_in, add_ni, pos_i
        module procedure add_cc, add_ci, add_ic, add_cr, add_rc
    end interface operator(+)

    interface operator(-)
        module procedure sub_ii, sub_ir, sub_ri, sub_in, sub_ni, neg_i
        module procedure sub_cc, sub_ci, sub_ic, sub_cr, sub_rc, neg_c
    end interface operator(-)

    interface operator(*)
        module procedure mul_ii, mul_ir, mul_ri, mul_in, mul_ni
        module procedure mul_cc, mul_ci, mul_ic, mul_cr, mul_rc
    end interface operator(*)

    interface operator(/)
        module procedure div_ii, div_ir, div_ri, div_in, div_ni
        module procedure div_cc, div_ci, div_ic, div_cr
    end interface operator(/)

    interface operator(**)
        module procedure pow_in, pow_cn
    end interface operator(**)

    interface sqrt
        module procedure sqrt_i, sqrt_c
    end interface sqrt

    interface exp
        module procedure exp_i, exp_c
    end interface exp

    interface log
        module procedure log_i
    end interface log

    interface sin
        module procedure sin_i, sin_c
    end interface sin

    interface cos
        module procedure cos_i, cos_c
    end interface cos

    interface sinh
        module procedure sinh_i
    end interface sinh

    interface cosh
        module procedure cosh_i
    end interface cosh

    interface abs
        module procedure abs_i
    end interface abs

    interface sqr
        module procedure sqr_i, sqr_c
    end interface sqr

    interface conjg
        module procedure conjg_c
    end interface conjg

    interface contains
        module procedure contains_ir, contains_ii
    end interface contains

    ! Two-word constants: c = c_hi + c_lo + r with |r| < ulp(c_lo), so
    ! [c_hi + round_down(c_lo), c_hi + round_up(c_lo)] encloses c. Each c_hi
    ! is the binary64 value nearest c and lies below c.
    real(dp), parameter :: pi_hi = 3.141592653589793_dp
    real(dp), parameter :: pi_lo = 1.2246467991473532e-16_dp
    real(dp), parameter :: pio2_hi = 1.5707963267948966_dp
    real(dp), parameter :: pio2_lo = 6.123233995736766e-17_dp
    real(dp), parameter :: ln2_hi = 0.6931471805599453_dp
    real(dp), parameter :: ln2_lo = 2.3190468138462996e-17_dp
    ! Cody-Waite split c_hi = c1 + c2: c1 keeps the leading 30 significand
    ! bits, so k c1 and k c2 are exact for abs(k) < 2^23.
    integer(int64), parameter :: cw_mask = not(int(z'7FFFFF', int64))
    real(dp), parameter :: pio2_c1 = transfer(iand(transfer(pio2_hi, 0_int64), &
        cw_mask), 1.0_dp)
    real(dp), parameter :: pio2_c2 = pio2_hi - pio2_c1
    real(dp), parameter :: ln2_c1 = transfer(iand(transfer(ln2_hi, 0_int64), &
        cw_mask), 1.0_dp)
    real(dp), parameter :: ln2_c2 = ln2_hi - ln2_c1
    integer, parameter :: k_exact = 2**23
    integer, parameter :: exp_degree = 18
    integer, parameter :: trig_degree = 26

contains

    pure elemental function iv_from_real(x) result(r)
        real(dp), intent(in) :: x
        type(interval_t) :: r

        r%lo = x
        r%hi = x
    end function iv_from_real

    !> The interval [min(lo, hi), max(lo, hi)].
    pure elemental function iv_from_bounds(lo, hi) result(r)
        real(dp), intent(in) :: lo, hi
        type(interval_t) :: r

        r%lo = min(lo, hi)
        r%hi = max(lo, hi)
    end function iv_from_bounds

    !> Point interval of an integer; outward rounded when abs(n) > 2^53.
    pure elemental function iv_from_int(n) result(r)
        integer, intent(in) :: n
        type(interval_t) :: r

        r%lo = real(n, dp)
        r%hi = r%lo
        if (abs(r%lo) > 2.0_dp**53) then
            r%lo = round_down(r%lo)
            r%hi = round_up(r%hi)
        end if
    end function iv_from_int

    pure function entire() result(r)
        type(interval_t) :: r

        r%hi = ieee_value(1.0_dp, ieee_positive_inf)
        r%lo = -r%hi
    end function entire

    !> True for an interval with a NaN endpoint or lo > hi.
    pure elemental function is_empty(a) result(e)
        type(interval_t), intent(in) :: a
        logical :: e

        e = .not. (a%lo <= a%hi)
    end function is_empty

    pure elemental function pos_i(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r

        r = a
    end function pos_i

    pure elemental function add_ii(a, b) result(r)
        type(interval_t), intent(in) :: a, b
        type(interval_t) :: r

        r%lo = round_down(a%lo + b%lo)
        r%hi = round_up(a%hi + b%hi)
    end function add_ii

    pure elemental function add_ir(a, b) result(r)
        type(interval_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(interval_t) :: r

        r = add_ii(a, iv_from_real(b))
    end function add_ir

    pure elemental function add_ri(a, b) result(r)
        real(dp), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(interval_t) :: r

        r = add_ii(iv_from_real(a), b)
    end function add_ri

    pure elemental function add_in(a, n) result(r)
        type(interval_t), intent(in) :: a
        integer, intent(in) :: n
        type(interval_t) :: r

        r = add_ii(a, iv_from_int(n))
    end function add_in

    pure elemental function add_ni(n, b) result(r)
        integer, intent(in) :: n
        type(interval_t), intent(in) :: b
        type(interval_t) :: r

        r = add_ii(iv_from_int(n), b)
    end function add_ni

    pure elemental function neg_i(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r

        r%lo = -a%hi
        r%hi = -a%lo
    end function neg_i

    pure elemental function sub_ii(a, b) result(r)
        type(interval_t), intent(in) :: a, b
        type(interval_t) :: r

        r%lo = round_down(a%lo - b%hi)
        r%hi = round_up(a%hi - b%lo)
    end function sub_ii

    pure elemental function sub_ir(a, b) result(r)
        type(interval_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(interval_t) :: r

        r = sub_ii(a, iv_from_real(b))
    end function sub_ir

    pure elemental function sub_ri(a, b) result(r)
        real(dp), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(interval_t) :: r

        r = sub_ii(iv_from_real(a), b)
    end function sub_ri

    pure elemental function sub_in(a, n) result(r)
        type(interval_t), intent(in) :: a
        integer, intent(in) :: n
        type(interval_t) :: r

        r = sub_ii(a, iv_from_int(n))
    end function sub_in

    pure elemental function sub_ni(n, b) result(r)
        integer, intent(in) :: n
        type(interval_t), intent(in) :: b
        type(interval_t) :: r

        r = sub_ii(iv_from_int(n), b)
    end function sub_ni

    !> Product from the four endpoint products; 0 * inf gives the entire line.
    pure elemental function mul_ii(a, b) result(r)
        type(interval_t), intent(in) :: a, b
        type(interval_t) :: r
        real(dp) :: p1, p2, p3, p4

        p1 = a%lo*b%lo
        p2 = a%lo*b%hi
        p3 = a%hi*b%lo
        p4 = a%hi*b%hi
        if (ieee_is_nan(p1) .or. ieee_is_nan(p2) .or. ieee_is_nan(p3) &
            .or. ieee_is_nan(p4)) then
            r = entire()
            return
        end if
        r%lo = round_down(min(p1, p2, p3, p4))
        r%hi = round_up(max(p1, p2, p3, p4))
    end function mul_ii

    pure elemental function mul_ir(a, b) result(r)
        type(interval_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(interval_t) :: r

        r = mul_ii(a, iv_from_real(b))
    end function mul_ir

    pure elemental function mul_ri(a, b) result(r)
        real(dp), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(interval_t) :: r

        r = mul_ii(iv_from_real(a), b)
    end function mul_ri

    pure elemental function mul_in(a, n) result(r)
        type(interval_t), intent(in) :: a
        integer, intent(in) :: n
        type(interval_t) :: r

        r = mul_ii(a, iv_from_int(n))
    end function mul_in

    pure elemental function mul_ni(n, b) result(r)
        integer, intent(in) :: n
        type(interval_t), intent(in) :: b
        type(interval_t) :: r

        r = mul_ii(iv_from_int(n), b)
    end function mul_ni

    !> Quotient; a divisor that contains zero gives the entire line.
    pure elemental function div_ii(a, b) result(r)
        type(interval_t), intent(in) :: a, b
        type(interval_t) :: r
        real(dp) :: q1, q2, q3, q4

        if (.not. (b%lo > 0.0_dp .or. b%hi < 0.0_dp)) then
            r = entire()
            return
        end if
        q1 = a%lo/b%lo
        q2 = a%lo/b%hi
        q3 = a%hi/b%lo
        q4 = a%hi/b%hi
        if (ieee_is_nan(q1) .or. ieee_is_nan(q2) .or. ieee_is_nan(q3) &
            .or. ieee_is_nan(q4)) then
            r = entire()
            return
        end if
        r%lo = round_down(min(q1, q2, q3, q4))
        r%hi = round_up(max(q1, q2, q3, q4))
    end function div_ii

    pure elemental function div_ir(a, b) result(r)
        type(interval_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(interval_t) :: r

        r = div_ii(a, iv_from_real(b))
    end function div_ir

    pure elemental function div_ri(a, b) result(r)
        real(dp), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(interval_t) :: r

        r = div_ii(iv_from_real(a), b)
    end function div_ri

    pure elemental function div_in(a, n) result(r)
        type(interval_t), intent(in) :: a
        integer, intent(in) :: n
        type(interval_t) :: r

        r = div_ii(a, iv_from_int(n))
    end function div_in

    pure elemental function div_ni(n, b) result(r)
        integer, intent(in) :: n
        type(interval_t), intent(in) :: b
        type(interval_t) :: r

        r = div_ii(iv_from_int(n), b)
    end function div_ni

    !> Square without the dependency problem of a*a.
    pure elemental function sqr_i(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r
        real(dp) :: m

        if (a%lo >= 0.0_dp) then
            r%lo = round_down(a%lo*a%lo)
            r%hi = round_up(a%hi*a%hi)
        else if (a%hi <= 0.0_dp) then
            r%lo = round_down(a%hi*a%hi)
            r%hi = round_up(a%lo*a%lo)
        else
            m = max(-a%lo, a%hi)
            r%lo = 0.0_dp
            r%hi = round_up(m*m)
        end if
        r%lo = max(r%lo, 0.0_dp)
    end function sqr_i

    !> x^n for x >= 0 rounded toward direction s (+1 up, -1 down).
    pure elemental function pow_nonneg(x, n, s) result(y)
        real(dp), intent(in) :: x
        integer, intent(in) :: n, s
        real(dp) :: y
        integer :: k

        y = 1.0_dp
        do k = 1, n
            if (s > 0) then
                y = round_up(y*x)
            else
                y = max(round_down(y*x), 0.0_dp)
            end if
        end do
    end function pow_nonneg

    !> Integer power, exact in the monotonicity of x^n (no dependency loss).
    pure elemental function pow_in(a, n) result(r)
        type(interval_t), intent(in) :: a
        integer, intent(in) :: n
        type(interval_t) :: r

        if (n < 0) then
            r = div_ii(iv_from_real(1.0_dp), pow_nonneg_exponent(a, -n))
        else
            r = pow_nonneg_exponent(a, n)
        end if
    end function pow_in

    pure elemental function pow_nonneg_exponent(a, n) result(r)
        type(interval_t), intent(in) :: a
        integer, intent(in) :: n
        type(interval_t) :: r

        if (n == 0) then
            r = iv_from_real(1.0_dp)
        else if (a%lo >= 0.0_dp) then
            r%lo = pow_nonneg(a%lo, n, -1)
            r%hi = pow_nonneg(a%hi, n, 1)
        else if (a%hi <= 0.0_dp) then
            if (modulo(n, 2) == 0) then
                r%lo = pow_nonneg(-a%hi, n, -1)
                r%hi = pow_nonneg(-a%lo, n, 1)
            else
                r%lo = -pow_nonneg(-a%lo, n, 1)
                r%hi = -pow_nonneg(-a%hi, n, -1)
            end if
        else
            if (modulo(n, 2) == 0) then
                r%lo = 0.0_dp
                r%hi = pow_nonneg(max(-a%lo, a%hi), n, 1)
            else
                r%lo = -pow_nonneg(-a%lo, n, 1)
                r%hi = pow_nonneg(a%hi, n, 1)
            end if
        end if
    end function pow_nonneg_exponent

    !> Square root of the part of a inside [0, inf).
    pure elemental function sqrt_i(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r

        r%lo = sqrt_down(max(a%lo, 0.0_dp))
        r%hi = sqrt_up(max(a%hi, 0.0_dp))
    end function sqrt_i

    pure elemental function abs_i(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r

        if (a%lo >= 0.0_dp) then
            r = a
        else if (a%hi <= 0.0_dp) then
            r = neg_i(a)
        else
            r%lo = 0.0_dp
            r%hi = max(-a%lo, a%hi)
        end if
    end function abs_i

    !> Enclosure of pi.
    pure function interval_pi() result(r)
        type(interval_t) :: r

        r%lo = round_down(pi_hi + round_down(pi_lo))
        r%hi = round_up(pi_hi + round_up(pi_lo))
    end function interval_pi

    pure elemental function hull(a, b) result(r)
        type(interval_t), intent(in) :: a, b
        type(interval_t) :: r

        r%lo = min(a%lo, b%lo)
        r%hi = max(a%hi, b%hi)
    end function hull

    !> Intersection; empty (lo > hi) when a and b are disjoint.
    pure elemental function intersect(a, b) result(r)
        type(interval_t), intent(in) :: a, b
        type(interval_t) :: r

        r%lo = max(a%lo, b%lo)
        r%hi = min(a%hi, b%hi)
    end function intersect

    !> A point of the interval near its centre.
    pure elemental function mid(a) result(m)
        type(interval_t), intent(in) :: a
        real(dp) :: m

        m = min(max(0.5_dp*a%lo + 0.5_dp*a%hi, a%lo), a%hi)
    end function mid

    !> Upper bound on the distance from mid(a) to either endpoint.
    pure elemental function rad(a) result(r)
        type(interval_t), intent(in) :: a
        real(dp) :: r, m

        m = mid(a)
        r = max(round_up(m - a%lo), round_up(a%hi - m))
    end function rad

    !> Upper bound on hi - lo.
    pure elemental function width(a) result(w)
        type(interval_t), intent(in) :: a
        real(dp) :: w

        w = round_up(a%hi - a%lo)
    end function width

    !> Magnitude max(abs(x)) over the interval.
    pure elemental function mag(a) result(m)
        type(interval_t), intent(in) :: a
        real(dp) :: m

        m = max(abs(a%lo), abs(a%hi))
    end function mag

    !> Mignitude min(abs(x)) over the interval.
    pure elemental function mig(a) result(m)
        type(interval_t), intent(in) :: a
        real(dp) :: m

        if (a%lo > 0.0_dp) then
            m = a%lo
        else if (a%hi < 0.0_dp) then
            m = -a%hi
        else
            m = 0.0_dp
        end if
    end function mig

    pure elemental function contains_ir(a, x) result(c)
        type(interval_t), intent(in) :: a
        real(dp), intent(in) :: x
        logical :: c

        c = a%lo <= x .and. x <= a%hi
    end function contains_ir

    !> True when b is a subset of a.
    pure elemental function contains_ii(a, b) result(c)
        type(interval_t), intent(in) :: a, b
        logical :: c

        c = a%lo <= b%lo .and. b%hi <= a%hi
    end function contains_ii

    pure elemental function contains_zero(a) result(c)
        type(interval_t), intent(in) :: a
        logical :: c

        c = a%lo <= 0.0_dp .and. 0.0_dp <= a%hi
    end function contains_zero

    !> True when a is a subset of b.
    pure elemental function subset(a, b) result(c)
        type(interval_t), intent(in) :: a, b
        logical :: c

        c = b%lo <= a%lo .and. a%hi <= b%hi
    end function subset

    pure elemental function disjoint(a, b) result(d)
        type(interval_t), intent(in) :: a, b
        logical :: d

        d = a%hi < b%lo .or. b%hi < a%lo
    end function disjoint

    !> Upper bound on abs(r)^n/n! * c for c >= 0.
    pure function term_bound(c, r, n) result(t)
        real(dp), intent(in) :: c, r
        integer, intent(in) :: n
        real(dp) :: t
        integer :: j

        t = c
        do j = 1, n
            t = round_up(round_up(t*r)/real(j, dp))
        end do
    end function term_bound

    !> Enclosure of exp(r) for abs(r) <= 0.4: degree-exp_degree Taylor
    !> polynomial in Horner form plus the Lagrange remainder
    !> exp(rho) rho^(N+1)/(N+1)! <= 1.5 rho^(N+1)/(N+1)!.
    pure function exp_reduced(r) result(p)
        type(interval_t), intent(in) :: r
        type(interval_t) :: p
        real(dp) :: e
        integer :: j

        p = iv_from_real(1.0_dp)
        do j = exp_degree, 1, -1
            p = add_ri(1.0_dp, div_in(mul_ii(r, p), j))
        end do
        e = term_bound(1.5_dp, mag(r), exp_degree + 1)
        p%lo = round_down(p%lo - e)
        p%hi = round_up(p%hi + e)
    end function exp_reduced

    !> Enclosure of exp(x) for a point x: x = k ln2 + r, exp(x) = 2^k exp(r).
    pure function exp_point(x) result(r)
        real(dp), intent(in) :: x
        type(interval_t) :: r, red, lnlo
        real(dp) :: kr

        if (ieee_is_nan(x)) then
            r = entire()
        else if (x > 710.0_dp) then
            r%lo = huge(1.0_dp)
            r%hi = ieee_value(1.0_dp, ieee_positive_inf)
        else if (x < -746.0_dp) then
            r%lo = 0.0_dp
            r%hi = tiny(1.0_dp)*epsilon(1.0_dp)
        else
            kr = anint(x/ln2_hi)
            lnlo = interval(round_down(ln2_lo), round_up(ln2_lo))
            red = sub_ii(iv_from_real(x), iv_from_real(kr*ln2_c1))
            red = sub_ir(red, kr*ln2_c2)
            red = sub_ii(red, mul_ri(kr, lnlo))
            r = exp_reduced(red)
            r = scale_pow2(r, nint(kr))
        end if
    end function exp_point

    !> Multiply an interval with nonnegative endpoints by 2^k.
    pure function scale_pow2(a, k) result(r)
        type(interval_t), intent(in) :: a
        integer, intent(in) :: k
        type(interval_t) :: r

        r%lo = scale(a%lo, k)
        r%hi = scale(a%hi, k)
        if (r%lo < tiny(1.0_dp)) r%lo = max(round_down(r%lo), 0.0_dp)
        if (r%hi < tiny(1.0_dp)) r%hi = round_up(r%hi)
        if (r%lo > huge(1.0_dp)) r%lo = huge(1.0_dp)
    end function scale_pow2

    pure elemental function exp_i(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r, t

        t = exp_point(a%lo)
        r%lo = t%lo
        t = exp_point(a%hi)
        r%hi = t%hi
    end function exp_i

    !> Lower (s < 0) or upper (s > 0) bound on log(x), x > 0, verified by
    !> exp: a y with exp(y) <= x is a lower bound, exp(y) >= x an upper one.
    !> The libm log only supplies the starting guess.
    pure function log_bound(x, s) result(y)
        real(dp), intent(in) :: x
        integer, intent(in) :: s
        real(dp) :: y, g, d
        type(interval_t) :: e
        integer :: it

        g = log(x)
        d = 2.0_dp*spacing(max(abs(g), 0.5_dp))
        do it = 1, 64
            if (s < 0) then
                y = g - d
                e = exp_point(y)
                if (e%hi <= x) return
            else
                y = g + d
                e = exp_point(y)
                if (e%lo >= x) return
            end if
            d = 2.0_dp*d
        end do
        y = ieee_value(1.0_dp, ieee_positive_inf)
        if (s < 0) y = -y
    end function log_bound

    !> Logarithm of the part of a inside (0, inf).
    pure elemental function log_i(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r, t

        if (.not. (a%hi > 0.0_dp)) then
            r%lo = ieee_value(1.0_dp, ieee_positive_inf)
            r%hi = -r%lo
            return
        end if
        if (a%lo > 0.0_dp) then
            r%lo = log_bound(a%lo, -1)
            if (a%lo >= 0.5_dp .and. a%lo <= 2.0_dp) then
                t = log_near_one(a%lo)
                r%lo = t%lo
            end if
        else
            r%lo = -ieee_value(1.0_dp, ieee_positive_inf)
        end if
        if (a%hi > huge(1.0_dp)) then
            r%hi = a%hi
        else if (a%hi >= 0.5_dp .and. a%hi <= 2.0_dp) then
            t = log_near_one(a%hi)
            r%hi = t%hi
        else
            r%hi = log_bound(a%hi, 1)
        end if
    end function log_i

    !> log(x) = 2 atanh(s), s = (x - 1)/(x + 1), for x in [0.5, 2]
    !> (abs(s) <= 1/3): 2 sum_(k<=K) s^(2k+1)/(2k+1) plus the remainder
    !> 2 abs(s)^(2K+3)/((2K+3)(1 - s^2)) <= 3 abs(s)^(2K+3)/(2K+3). Accurate
    !> relative to log(x) itself near x = 1.
    pure function log_near_one(x) result(r)
        real(dp), intent(in) :: x
        type(interval_t) :: r, s, s2
        real(dp) :: e, sm
        integer :: k
        integer, parameter :: kmax = 30

        s = div_ii(sub_ir(iv_from_real(x), 1.0_dp), add_ir(iv_from_real(x), 1.0_dp))
        s2 = sqr_i(s)
        r = div_in(iv_from_real(1.0_dp), 2*kmax + 1)
        do k = kmax - 1, 0, -1
            r = add_ii(div_in(iv_from_real(1.0_dp), 2*k + 1), mul_ii(s2, r))
        end do
        r = mul_in(mul_ii(s, r), 2)
        sm = mag(s)
        e = 3.0_dp
        do k = 1, 2*kmax + 3
            e = round_up(e*sm)
        end do
        e = round_up(e/real(2*kmax + 3, dp))
        r%lo = round_down(r%lo - e)
        r%hi = round_up(r%hi + e)
    end function log_near_one

    !> Enclosures of sin(r) and cos(r) by Taylor polynomials of degree
    !> 2 M + 1 and 2 M with Lagrange remainders abs(r)^(2M+3)/(2M+3)! and
    !> abs(r)^(2M+2)/(2M+2)! (all derivatives are bounded by 1).
    pure subroutine sincos_reduced(r, s, c)
        type(interval_t), intent(in) :: r
        type(interval_t), intent(out) :: s, c
        type(interval_t) :: r2
        real(dp) :: es, ec, rm
        integer :: j, m

        m = trig_degree/2
        r2 = sqr_i(r)
        rm = mag(r)
        s = iv_from_real(1.0_dp)
        c = iv_from_real(1.0_dp)
        do j = m, 1, -1
            s = sub_ri(1.0_dp, div_in(mul_ii(r2, s), (2*j)*(2*j + 1)))
            c = sub_ri(1.0_dp, div_in(mul_ii(r2, c), (2*j - 1)*(2*j)))
        end do
        s = mul_ii(r, s)
        es = term_bound(1.0_dp, rm, 2*m + 3)
        ec = term_bound(1.0_dp, rm, 2*m + 2)
        s%lo = round_down(s%lo - es)
        s%hi = round_up(s%hi + es)
        c%lo = round_down(c%lo - ec)
        c%hi = round_up(c%hi + ec)
    end subroutine sincos_reduced

    !> Enclosures of sin(x) and cos(x) for a point x: x = k pi/2 + r with a
    !> Cody-Waite reduction that is exact in its leading terms for
    !> abs(k) < 2^23 and rigorous (but wider) beyond.
    pure subroutine sincos_point(x, s, c)
        real(dp), intent(in) :: x
        type(interval_t), intent(out) :: s, c
        type(interval_t) :: r, plo, sr, cr
        real(dp) :: kr
        integer :: q

        s = interval(-1.0_dp, 1.0_dp)
        c = s
        if (ieee_is_nan(x) .or. abs(x) > 1.0e300_dp) return
        kr = anint(x/pio2_hi)
        plo = interval(round_down(pio2_lo), round_up(pio2_lo))
        if (abs(kr) < real(k_exact, dp)) then
            r = sub_ii(iv_from_real(x), iv_from_real(kr*pio2_c1))
            r = sub_ir(r, kr*pio2_c2)
            r = sub_ii(r, mul_ri(kr, plo))
        else
            r = sub_ri(x, mul_ri(kr, add_ri(pio2_hi, plo)))
        end if
        if (mag(r) > 1.0_dp) return
        q = nint(modulo(kr, 4.0_dp))
        call sincos_reduced(r, sr, cr)
        select case (q)
        case (0)
            s = sr
            c = cr
        case (1)
            s = cr
            c = neg_i(sr)
        case (2)
            s = neg_i(sr)
            c = neg_i(cr)
        case default
            s = neg_i(cr)
            c = sr
        end select
        s = clip_unit(s)
        c = clip_unit(c)
    end subroutine sincos_point

    pure elemental function clip_unit(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r

        r%lo = max(a%lo, -1.0_dp)
        r%hi = min(a%hi, 1.0_dp)
    end function clip_unit

    !> Range of sin (shift = 0.5) or cos (shift = 0) over a: hull of the
    !> endpoint enclosures plus every extremum (j + shift) pi that may lie in
    !> a (j even: maximum 1, j odd: minimum -1).
    pure function trig_range(a, shift) result(r)
        type(interval_t), intent(in) :: a
        real(dp), intent(in) :: shift
        type(interval_t) :: r, s1, c1, s2, c2, pi_iv, pj
        integer(int64) :: j, jlo, jhi

        r = interval(-1.0_dp, 1.0_dp)
        if (is_empty(a) .or. max(abs(a%lo), abs(a%hi)) > 1.0e15_dp) return
        if (a%hi - a%lo >= 6.5_dp) return
        call sincos_point(a%lo, s1, c1)
        call sincos_point(a%hi, s2, c2)
        if (shift > 0.0_dp) then
            r = hull(s1, s2)
        else
            r = hull(c1, c2)
        end if
        pi_iv = interval_pi()
        jlo = floor(a%lo/pi_hi - shift, int64) - 1
        jhi = ceiling(a%hi/pi_hi - shift, int64) + 1
        do j = jlo, jhi
            pj = mul_ri(real(j, dp) + shift, pi_iv)
            if (disjoint(pj, a)) cycle
            if (modulo(j, 2_int64) == 0) then
                r%hi = 1.0_dp
            else
                r%lo = -1.0_dp
            end if
        end do
    end function trig_range

    pure elemental function sin_i(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r

        r = trig_range(a, 0.5_dp)
    end function sin_i

    pure elemental function cos_i(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r

        r = trig_range(a, 0.0_dp)
    end function cos_i

    !> sinh at a point: Taylor for abs(x) < 1 (remainder bounded with
    !> cosh(1) < 1.6), (e - 1/e)/2 otherwise.
    pure function sinh_point(x) result(r)
        real(dp), intent(in) :: x
        type(interval_t) :: r, x2, e
        real(dp) :: t
        integer :: j, m

        if (abs(x) < 1.0_dp) then
            m = trig_degree/2
            x2 = sqr_i(iv_from_real(x))
            r = iv_from_real(1.0_dp)
            do j = m, 1, -1
                r = add_ri(1.0_dp, div_in(mul_ii(x2, r), (2*j)*(2*j + 1)))
            end do
            r = mul_ri(x, r)
            t = term_bound(1.6_dp, abs(x), 2*m + 3)
            r%lo = round_down(r%lo - t)
            r%hi = round_up(r%hi + t)
        else
            e = exp_point(x)
            r = div_in(sub_ii(e, div_ri(1.0_dp, e)), 2)
        end if
    end function sinh_point

    pure function cosh_point(x) result(r)
        real(dp), intent(in) :: x
        type(interval_t) :: r, e

        e = exp_point(abs(x))
        r = div_in(add_ii(e, div_ri(1.0_dp, e)), 2)
        r%lo = max(r%lo, 1.0_dp)
    end function cosh_point

    pure elemental function sinh_i(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r, t

        t = sinh_point(a%lo)
        r%lo = t%lo
        t = sinh_point(a%hi)
        r%hi = t%hi
    end function sinh_i

    pure elemental function cosh_i(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r, t1, t2

        t1 = cosh_point(a%lo)
        t2 = cosh_point(a%hi)
        if (contains_zero(a)) then
            r%lo = 1.0_dp
            r%hi = max(t1%hi, t2%hi)
        else if (a%lo > 0.0_dp) then
            r%lo = t1%lo
            r%hi = t2%hi
        else
            r%lo = t2%lo
            r%hi = t1%hi
        end if
    end function cosh_i

    pure elemental function ci_from_complex(z) result(r)
        complex(dp), intent(in) :: z
        type(cinterval_t) :: r

        r%re = iv_from_real(real(z, dp))
        r%im = iv_from_real(aimag(z))
    end function ci_from_complex

    pure elemental function ci_from_parts(re, im) result(r)
        type(interval_t), intent(in) :: re, im
        type(cinterval_t) :: r

        r%re = re
        r%im = im
    end function ci_from_parts

    pure elemental function ci_from_real_parts(re, im) result(r)
        real(dp), intent(in) :: re, im
        type(cinterval_t) :: r

        r%re = iv_from_real(re)
        r%im = iv_from_real(im)
    end function ci_from_real_parts

    pure elemental function real_part(a) result(r)
        type(cinterval_t), intent(in) :: a
        type(interval_t) :: r

        r = a%re
    end function real_part

    pure elemental function imag_part(a) result(r)
        type(cinterval_t), intent(in) :: a
        type(interval_t) :: r

        r = a%im
    end function imag_part

    pure elemental function conjg_c(a) result(r)
        type(cinterval_t), intent(in) :: a
        type(cinterval_t) :: r

        r%re = a%re
        r%im = neg_i(a%im)
    end function conjg_c

    pure elemental function add_cc(a, b) result(r)
        type(cinterval_t), intent(in) :: a, b
        type(cinterval_t) :: r

        r%re = add_ii(a%re, b%re)
        r%im = add_ii(a%im, b%im)
    end function add_cc

    pure elemental function add_ci(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(cinterval_t) :: r

        r%re = add_ii(a%re, b)
        r%im = a%im
    end function add_ci

    pure elemental function add_ic(a, b) result(r)
        type(interval_t), intent(in) :: a
        type(cinterval_t), intent(in) :: b
        type(cinterval_t) :: r

        r = add_ci(b, a)
    end function add_ic

    pure elemental function add_cr(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(cinterval_t) :: r

        r = add_ci(a, iv_from_real(b))
    end function add_cr

    pure elemental function add_rc(a, b) result(r)
        real(dp), intent(in) :: a
        type(cinterval_t), intent(in) :: b
        type(cinterval_t) :: r

        r = add_ci(b, iv_from_real(a))
    end function add_rc

    pure elemental function neg_c(a) result(r)
        type(cinterval_t), intent(in) :: a
        type(cinterval_t) :: r

        r%re = neg_i(a%re)
        r%im = neg_i(a%im)
    end function neg_c

    pure elemental function sub_cc(a, b) result(r)
        type(cinterval_t), intent(in) :: a, b
        type(cinterval_t) :: r

        r%re = sub_ii(a%re, b%re)
        r%im = sub_ii(a%im, b%im)
    end function sub_cc

    pure elemental function sub_ci(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(cinterval_t) :: r

        r%re = sub_ii(a%re, b)
        r%im = a%im
    end function sub_ci

    pure elemental function sub_ic(a, b) result(r)
        type(interval_t), intent(in) :: a
        type(cinterval_t), intent(in) :: b
        type(cinterval_t) :: r

        r%re = sub_ii(a, b%re)
        r%im = neg_i(b%im)
    end function sub_ic

    pure elemental function sub_cr(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(cinterval_t) :: r

        r = sub_ci(a, iv_from_real(b))
    end function sub_cr

    pure elemental function sub_rc(a, b) result(r)
        real(dp), intent(in) :: a
        type(cinterval_t), intent(in) :: b
        type(cinterval_t) :: r

        r = sub_ic(iv_from_real(a), b)
    end function sub_rc

    pure elemental function mul_cc(a, b) result(r)
        type(cinterval_t), intent(in) :: a, b
        type(cinterval_t) :: r

        r%re = sub_ii(mul_ii(a%re, b%re), mul_ii(a%im, b%im))
        r%im = add_ii(mul_ii(a%re, b%im), mul_ii(a%im, b%re))
    end function mul_cc

    pure elemental function mul_ci(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(cinterval_t) :: r

        r%re = mul_ii(a%re, b)
        r%im = mul_ii(a%im, b)
    end function mul_ci

    pure elemental function mul_ic(a, b) result(r)
        type(interval_t), intent(in) :: a
        type(cinterval_t), intent(in) :: b
        type(cinterval_t) :: r

        r = mul_ci(b, a)
    end function mul_ic

    pure elemental function mul_cr(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(cinterval_t) :: r

        r = mul_ci(a, iv_from_real(b))
    end function mul_cr

    pure elemental function mul_rc(a, b) result(r)
        real(dp), intent(in) :: a
        type(cinterval_t), intent(in) :: b
        type(cinterval_t) :: r

        r = mul_ci(b, iv_from_real(a))
    end function mul_rc

    !> a/b = a conj(b)/abs(b)^2; a box b that contains 0 gives entire parts.
    pure elemental function div_cc(a, b) result(r)
        type(cinterval_t), intent(in) :: a, b
        type(cinterval_t) :: r
        type(interval_t) :: den

        den = add_ii(sqr_i(b%re), sqr_i(b%im))
        r%re = div_ii(add_ii(mul_ii(a%re, b%re), mul_ii(a%im, b%im)), den)
        r%im = div_ii(sub_ii(mul_ii(a%im, b%re), mul_ii(a%re, b%im)), den)
    end function div_cc

    pure elemental function div_ci(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        type(interval_t), intent(in) :: b
        type(cinterval_t) :: r

        r%re = div_ii(a%re, b)
        r%im = div_ii(a%im, b)
    end function div_ci

    pure elemental function div_ic(a, b) result(r)
        type(interval_t), intent(in) :: a
        type(cinterval_t), intent(in) :: b
        type(cinterval_t) :: r

        r = div_cc(ci_from_parts(a, iv_from_real(0.0_dp)), b)
    end function div_ic

    pure elemental function div_cr(a, b) result(r)
        type(cinterval_t), intent(in) :: a
        real(dp), intent(in) :: b
        type(cinterval_t) :: r

        r = div_ci(a, iv_from_real(b))
    end function div_cr

    pure elemental function sqr_c(a) result(r)
        type(cinterval_t), intent(in) :: a
        type(cinterval_t) :: r

        r%re = sub_ii(sqr_i(a%re), sqr_i(a%im))
        r%im = mul_in(mul_ii(a%re, a%im), 2)
    end function sqr_c

    !> Integer power by binary exponentiation.
    pure elemental function pow_cn(a, n) result(r)
        type(cinterval_t), intent(in) :: a
        integer, intent(in) :: n
        type(cinterval_t) :: r, base
        integer :: m

        r = ci_from_real_parts(1.0_dp, 0.0_dp)
        base = a
        m = abs(n)
        do while (m > 0)
            if (modulo(m, 2) == 1) r = mul_cc(r, base)
            m = m/2
            if (m > 0) base = sqr_c(base)
        end do
        if (n < 0) r = div_cc(ci_from_real_parts(1.0_dp, 0.0_dp), r)
    end function pow_cn

    !> Principal square root on boxes with re > 0:
    !> sqrt(z) = t + i y/(2 t), t = sqrt((abs(z) + x)/2). Other boxes give
    !> entire parts.
    pure elemental function sqrt_c(a) result(r)
        type(cinterval_t), intent(in) :: a
        type(cinterval_t) :: r
        type(interval_t) :: m, t

        if (.not. (a%re%lo > 0.0_dp)) then
            r%re = entire()
            r%im = entire()
            return
        end if
        m = sqrt_i(add_ii(sqr_i(a%re), sqr_i(a%im)))
        t = sqrt_i(div_in(add_ii(m, a%re), 2))
        r%re = t
        r%im = div_ii(a%im, mul_in(t, 2))
    end function sqrt_c

    !> exp(x + i y) = exp(x) (cos y + i sin y).
    pure elemental function exp_c(a) result(r)
        type(cinterval_t), intent(in) :: a
        type(cinterval_t) :: r
        type(interval_t) :: e

        e = exp_i(a%re)
        r%re = mul_ii(e, cos_i(a%im))
        r%im = mul_ii(e, sin_i(a%im))
    end function exp_c

    !> sin(x + i y) = sin x cosh y + i cos x sinh y.
    pure elemental function sin_c(a) result(r)
        type(cinterval_t), intent(in) :: a
        type(cinterval_t) :: r

        r%re = mul_ii(sin_i(a%re), cosh_i(a%im))
        r%im = mul_ii(cos_i(a%re), sinh_i(a%im))
    end function sin_c

    !> cos(x + i y) = cos x cosh y - i sin x sinh y.
    pure elemental function cos_c(a) result(r)
        type(cinterval_t), intent(in) :: a
        type(cinterval_t) :: r

        r%re = mul_ii(cos_i(a%re), cosh_i(a%im))
        r%im = neg_i(mul_ii(sin_i(a%re), sinh_i(a%im)))
    end function cos_c

    !> Upper bound on abs(z) over the box.
    pure elemental function cabs_up(a) result(m)
        type(cinterval_t), intent(in) :: a
        real(dp) :: m
        type(interval_t) :: t

        t = sqrt_i(add_ii(sqr_i(a%re), sqr_i(a%im)))
        m = t%hi
    end function cabs_up

    !> Lower bound on abs(z) over the box.
    pure elemental function cabs_down(a) result(m)
        type(cinterval_t), intent(in) :: a
        real(dp) :: m
        type(interval_t) :: t

        t = sqrt_i(add_ii(sqr_i(a%re), sqr_i(a%im)))
        m = t%lo
    end function cabs_down

    ! Runtime interface for fortsym-emitted kernels. Semantics per the
    ! fortsym rigorous-runtime contract: every result encloses the exact
    ! result for all operand values; an operation outside its domain returns
    ! the entire line (isqrt of an interval with a negative part, iinv and
    ! idiv of an interval containing zero).

    pure elemental function ipoint(x) result(r)
        real(dp), intent(in) :: x
        type(interval_t) :: r

        r = iv_from_real(x)
    end function ipoint

    !> The interval [m - r, m + r] rounded outward.
    pure elemental function ienclose(m, r) result(a)
        real(dp), intent(in) :: m, r
        type(interval_t) :: a

        a%lo = round_down(m - r)
        a%hi = round_up(m + r)
    end function ienclose

    pure elemental function iadd(a, b) result(r)
        type(interval_t), intent(in) :: a, b
        type(interval_t) :: r

        r = add_ii(a, b)
    end function iadd

    pure elemental function isub(a, b) result(r)
        type(interval_t), intent(in) :: a, b
        type(interval_t) :: r

        r = sub_ii(a, b)
    end function isub

    pure elemental function imul(a, b) result(r)
        type(interval_t), intent(in) :: a, b
        type(interval_t) :: r

        r = mul_ii(a, b)
    end function imul

    pure elemental function idiv(a, b) result(r)
        type(interval_t), intent(in) :: a, b
        type(interval_t) :: r

        r = div_ii(a, b)
    end function idiv

    pure elemental function ineg(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r

        r = neg_i(a)
    end function ineg

    pure elemental function iinv(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r

        r = div_ii(iv_from_real(1.0_dp), a)
    end function iinv

    pure elemental function isqrt(a) result(r)
        type(interval_t), intent(in) :: a
        type(interval_t) :: r

        if (a%lo < 0.0_dp) then
            r = entire()
        else
            r = sqrt_i(a)
        end if
    end function isqrt

    pure elemental function ipowi(a, n) result(r)
        type(interval_t), intent(in) :: a
        integer, intent(in) :: n
        type(interval_t) :: r

        r = pow_in(a, n)
    end function ipowi

    pure elemental function iscale(a, x) result(r)
        type(interval_t), intent(in) :: a
        real(dp), intent(in) :: x
        type(interval_t) :: r

        r = mul_ir(a, x)
    end function iscale
end module fortnum_interval
