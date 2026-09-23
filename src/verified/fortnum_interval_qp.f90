!> Outward-rounded binary128 (real128) interval arithmetic.
!>
!> Used where binary64 enclosures are too wide: closed-form Gaussian integrals
!> with heavy cancellation, and certified FFT twiddle factors. Endpoints are
!> computed with correctly rounded IEEE binary128 +, -, *, / and moved one ulp
!> outward with the exact nearest intrinsic. sqrt uses the library value only
!> as a guess and verifies each endpoint by squaring, so no library accuracy
!> claim enters.
module fortnum_interval_qp
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf
    use fortnum_rounding, only: round_up, round_down
    use fortnum_interval, only: interval_t
    implicit none
    private

    public :: qp, qinterval_t, qinterval, qrat, qinterval_pi, qsqrt, qsqr
    public :: to_interval, qmid, qrad
    public :: operator(+), operator(-), operator(*), operator(/)

    type :: qinterval_t
        real(qp) :: lo = 0.0_qp
        real(qp) :: hi = 0.0_qp
    end type qinterval_t

    interface qinterval
        module procedure q_from_real, q_from_bounds, q_from_int, q_from_dp
    end interface qinterval

    interface operator(+)
        module procedure qadd, qadd_n
    end interface operator(+)

    interface operator(-)
        module procedure qsub, qsub_n, qneg
    end interface operator(-)

    interface operator(*)
        module procedure qmul, qmul_n
    end interface operator(*)

    interface operator(/)
        module procedure qdiv, qdiv_n
    end interface operator(/)

    ! pi = pi_q + r with abs(r) < ulp(pi_q): pi_q is the binary128 value
    ! nearest pi.
    real(qp), parameter :: pi_q = 3.14159265358979323846264338327950288_qp

contains

    pure elemental function qdn(x) result(y)
        real(qp), intent(in) :: x
        real(qp) :: y

        y = x
        if (abs(x) <= huge(x)) y = nearest(x, -1.0_qp)
    end function qdn

    pure elemental function qup(x) result(y)
        real(qp), intent(in) :: x
        real(qp) :: y

        y = x
        if (abs(x) <= huge(x)) y = nearest(x, 1.0_qp)
    end function qup

    pure elemental function q_from_real(x) result(r)
        real(qp), intent(in) :: x
        type(qinterval_t) :: r

        r%lo = x
        r%hi = x
    end function q_from_real

    pure elemental function q_from_dp(x) result(r)
        real(dp), intent(in) :: x
        type(qinterval_t) :: r

        r%lo = real(x, qp)
        r%hi = r%lo
    end function q_from_dp

    pure elemental function q_from_bounds(lo, hi) result(r)
        real(qp), intent(in) :: lo, hi
        type(qinterval_t) :: r

        r%lo = min(lo, hi)
        r%hi = max(lo, hi)
    end function q_from_bounds

    pure elemental function q_from_int(n) result(r)
        integer, intent(in) :: n
        type(qinterval_t) :: r

        r%lo = real(n, qp)
        r%hi = r%lo
    end function q_from_int

    !> Enclosure of the rational n/m, m /= 0.
    pure elemental function qrat(n, m) result(r)
        integer, intent(in) :: n, m
        type(qinterval_t) :: r

        r = qdiv(q_from_int(n), q_from_int(m))
    end function qrat

    pure function qinterval_pi() result(r)
        type(qinterval_t) :: r

        r%lo = qdn(pi_q)
        r%hi = qup(pi_q)
    end function qinterval_pi

    pure elemental function qadd(a, b) result(r)
        type(qinterval_t), intent(in) :: a, b
        type(qinterval_t) :: r

        r%lo = qdn(a%lo + b%lo)
        r%hi = qup(a%hi + b%hi)
    end function qadd

    pure elemental function qadd_n(a, n) result(r)
        type(qinterval_t), intent(in) :: a
        integer, intent(in) :: n
        type(qinterval_t) :: r

        r = qadd(a, q_from_int(n))
    end function qadd_n

    pure elemental function qneg(a) result(r)
        type(qinterval_t), intent(in) :: a
        type(qinterval_t) :: r

        r%lo = -a%hi
        r%hi = -a%lo
    end function qneg

    pure elemental function qsub(a, b) result(r)
        type(qinterval_t), intent(in) :: a, b
        type(qinterval_t) :: r

        r%lo = qdn(a%lo - b%hi)
        r%hi = qup(a%hi - b%lo)
    end function qsub

    pure elemental function qsub_n(a, n) result(r)
        type(qinterval_t), intent(in) :: a
        integer, intent(in) :: n
        type(qinterval_t) :: r

        r = qsub(a, q_from_int(n))
    end function qsub_n

    pure elemental function qmul(a, b) result(r)
        type(qinterval_t), intent(in) :: a, b
        type(qinterval_t) :: r
        real(qp) :: p1, p2, p3, p4

        p1 = a%lo*b%lo
        p2 = a%lo*b%hi
        p3 = a%hi*b%lo
        p4 = a%hi*b%hi
        r%lo = qdn(min(p1, p2, p3, p4))
        r%hi = qup(max(p1, p2, p3, p4))
    end function qmul

    pure elemental function qmul_n(a, n) result(r)
        type(qinterval_t), intent(in) :: a
        integer, intent(in) :: n
        type(qinterval_t) :: r

        r = qmul(a, q_from_int(n))
    end function qmul_n

    !> Quotient; a divisor containing zero gives the entire line.
    pure elemental function qdiv(a, b) result(r)
        type(qinterval_t), intent(in) :: a, b
        type(qinterval_t) :: r
        real(qp) :: q1, q2, q3, q4

        if (.not. (b%lo > 0.0_qp .or. b%hi < 0.0_qp)) then
            r%hi = ieee_value(1.0_qp, ieee_positive_inf)
            r%lo = -r%hi
            return
        end if
        q1 = a%lo/b%lo
        q2 = a%lo/b%hi
        q3 = a%hi/b%lo
        q4 = a%hi/b%hi
        r%lo = qdn(min(q1, q2, q3, q4))
        r%hi = qup(max(q1, q2, q3, q4))
    end function qdiv

    pure elemental function qdiv_n(a, n) result(r)
        type(qinterval_t), intent(in) :: a
        integer, intent(in) :: n
        type(qinterval_t) :: r

        r = qdiv(a, q_from_int(n))
    end function qdiv_n

    pure elemental function qsqr(a) result(r)
        type(qinterval_t), intent(in) :: a
        type(qinterval_t) :: r
        real(qp) :: m

        if (a%lo >= 0.0_qp) then
            r%lo = qdn(a%lo*a%lo)
            r%hi = qup(a%hi*a%hi)
        else if (a%hi <= 0.0_qp) then
            r%lo = qdn(a%hi*a%hi)
            r%hi = qup(a%lo*a%lo)
        else
            m = max(-a%lo, a%hi)
            r%lo = 0.0_qp
            r%hi = qup(m*m)
        end if
        r%lo = max(r%lo, 0.0_qp)
    end function qsqr

    !> Square root of the part of a inside [0, inf). Each endpoint is the
    !> library value moved outward until its square, rounded the other way,
    !> brackets the argument.
    pure elemental function qsqrt(a) result(r)
        type(qinterval_t), intent(in) :: a
        type(qinterval_t) :: r
        real(qp) :: x, s
        integer :: it

        x = max(a%lo, 0.0_qp)
        s = qdn(sqrt(x))
        do it = 1, 8
            if (s <= 0.0_qp) exit
            if (qup(s*s) <= x) exit
            s = qdn(s)
        end do
        r%lo = max(s, 0.0_qp)
        x = max(a%hi, 0.0_qp)
        s = qup(sqrt(x))
        do it = 1, 8
            if (qdn(s*s) >= x) exit
            s = qup(s)
        end do
        r%hi = s
    end function qsqrt

    !> Binary64 enclosure of a binary128 interval.
    pure elemental function to_interval(a) result(r)
        type(qinterval_t), intent(in) :: a
        type(interval_t) :: r

        r%lo = real(a%lo, dp)
        r%hi = real(a%hi, dp)
        if (real(r%lo, qp) > a%lo) r%lo = round_down(r%lo)
        if (real(r%hi, qp) < a%hi) r%hi = round_up(r%hi)
    end function to_interval

    pure elemental function qmid(a) result(m)
        type(qinterval_t), intent(in) :: a
        real(qp) :: m

        m = min(max(0.5_qp*a%lo + 0.5_qp*a%hi, a%lo), a%hi)
    end function qmid

    !> Upper bound on the distance from qmid(a) to either endpoint.
    pure elemental function qrad(a) result(r)
        type(qinterval_t), intent(in) :: a
        real(qp) :: r, m

        m = qmid(a)
        r = max(qup(m - a%lo), qup(a%hi - m))
    end function qrad
end module fortnum_interval_qp
