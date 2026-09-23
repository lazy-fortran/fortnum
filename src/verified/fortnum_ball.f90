!> Complex ball arithmetic: a ball (c, r) is the closed disc abs(z - c) <= r.
!>
!> Every operation returns a ball that contains the exact result for all
!> inputs in the argument balls. Centres are computed in round-to-nearest.
!> The rounding error of each IEEE operation is bounded with the standard
!> model abs(fl(x op y) - x op y) <= u abs(x op y), u = 2^-53; complex
!> products by sqrt(5) u abs(a) abs(b) (Brent, Percival, Zimmermann, Math.
!> Comp. 76 (2007) 1469-1481), covered here by 3u; every radius is accumulated
!> with the upward successor step of fortnum_rounding. An absolute floor
!> eta = 2 tiny added to every radius covers underflow. Moduli are bounded by
!> a power-of-two scaled sqrt(x^2 + y^2), so no libm hypot enters.
!>
!> The procedure names badd, bsub, bmul, bdiv, bneg, binv, bsqrt, bpowi,
!> bscale, bpoint, bcpoint and benclose form the ball runtime interface that
!> fortsym-emitted rigorous kernels call.
module fortnum_ball
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_rounding, only: round_up, round_down, unit_roundoff
    use fortnum_interval, only: interval_t, cinterval_t
    implicit none
    private

    public :: ball_t, ball, badd, bsub, bmul, bdiv, bneg, binv, bsqrt, bpowi
    public :: bscale, bmulr, bpoint, bcpoint, benclose, bconjg
    public :: babs_hi, babs_lo, bre_lo, bre_hi, bim_lo, bim_hi
    public :: bcontains, bcontains_zero, cabs_hi, cabs_lo
    public :: ball_from_cinterval, cinterval_from_ball, ball_from_interval
    public :: operator(+), operator(-), operator(*), operator(/), operator(**)

    type :: ball_t
        complex(dp) :: c = (0.0_dp, 0.0_dp)
        real(dp) :: r = 0.0_dp
    end type ball_t

    interface bpoint
        module procedure bpoint_r, bpoint_c
    end interface bpoint

    interface benclose
        module procedure benclose_r, benclose_c
    end interface benclose

    interface operator(+)
        module procedure badd, add_br, add_rb, add_bz, add_zb, bpos
    end interface operator(+)

    interface operator(-)
        module procedure bsub, sub_br, sub_rb, sub_bz, sub_zb, bneg
    end interface operator(-)

    interface operator(*)
        module procedure bmul, bscale, mul_rb, mul_bz, mul_zb
    end interface operator(*)

    interface operator(/)
        module procedure bdiv, div_br
    end interface operator(/)

    interface operator(**)
        module procedure bpowi
    end interface operator(**)

    real(dp), parameter :: u = unit_roundoff
    real(dp), parameter :: eta = 2.0_dp*tiny(1.0_dp)

contains

    !> Upper bound on abs(z). The power-of-two scaling is exact except for
    !> subnormal loss in the smaller component (below 2^-1074 relative); the
    !> scaled sqrt(x^2 + y^2) has relative error below 2u + u^2 and four
    !> upward steps (each at least u relative) cover both.
    pure elemental function cabs_hi(z) result(a)
        complex(dp), intent(in) :: z
        real(dp) :: a, x, y, m
        integer :: e

        x = abs(real(z, dp))
        y = abs(aimag(z))
        m = max(x, y)
        if (.not. (m > 0.0_dp)) then
            a = m
            return
        end if
        if (m > huge(1.0_dp)) then
            a = m
            return
        end if
        e = exponent(m)
        x = scale(x, -e)
        y = scale(y, -e)
        a = round_up(round_up(round_up(round_up(sqrt(x*x + y*y)))))
        a = scale(a, e)
        if (a < tiny(1.0_dp)) a = round_up(a)
    end function cabs_hi

    !> Lower bound on abs(z), never negative.
    pure elemental function cabs_lo(z) result(a)
        complex(dp), intent(in) :: z
        real(dp) :: a, x, y, m
        integer :: e

        x = abs(real(z, dp))
        y = abs(aimag(z))
        m = max(x, y)
        if (.not. (m > 0.0_dp) .or. m > huge(1.0_dp)) then
            a = m
            return
        end if
        e = exponent(m)
        x = scale(x, -e)
        y = scale(y, -e)
        a = round_down(round_down(round_down(round_down(sqrt(x*x + y*y)))))
        a = scale(max(a, 0.0_dp), e)
        if (a < tiny(1.0_dp)) a = max(round_down(a), 0.0_dp)
    end function cabs_lo

    pure elemental function ball(c, r) result(b)
        complex(dp), intent(in) :: c
        real(dp), intent(in) :: r
        type(ball_t) :: b

        b%c = c
        b%r = r
    end function ball

    pure elemental function bpoint_r(x) result(b)
        real(dp), intent(in) :: x
        type(ball_t) :: b

        b%c = cmplx(x, 0.0_dp, dp)
        b%r = 0.0_dp
    end function bpoint_r

    pure elemental function bpoint_c(c) result(b)
        complex(dp), intent(in) :: c
        type(ball_t) :: b

        b%c = c
        b%r = 0.0_dp
    end function bpoint_c

    pure elemental function bcpoint(x, y) result(b)
        real(dp), intent(in) :: x, y
        type(ball_t) :: b

        b%c = cmplx(x, y, dp)
        b%r = 0.0_dp
    end function bcpoint

    pure elemental function benclose_r(m, r) result(b)
        real(dp), intent(in) :: m, r
        type(ball_t) :: b

        b%c = cmplx(m, 0.0_dp, dp)
        b%r = r
    end function benclose_r

    pure elemental function benclose_c(c, r) result(b)
        complex(dp), intent(in) :: c
        real(dp), intent(in) :: r
        type(ball_t) :: b

        b%c = c
        b%r = r
    end function benclose_c

    pure elemental function bpos(a) result(b)
        type(ball_t), intent(in) :: a
        type(ball_t) :: b

        b = a
    end function bpos

    pure elemental function bneg(a) result(b)
        type(ball_t), intent(in) :: a
        type(ball_t) :: b

        b%c = -a%c
        b%r = a%r
    end function bneg

    pure elemental function bconjg(a) result(b)
        type(ball_t), intent(in) :: a
        type(ball_t) :: b

        b%c = conjg(a%c)
        b%r = a%r
    end function bconjg

    pure elemental function badd(a, b) result(s)
        type(ball_t), intent(in) :: a, b
        type(ball_t) :: s

        s%c = a%c + b%c
        s%r = round_up(round_up(round_up(a%r + b%r) &
            + round_up(2.0_dp*u*cabs_hi(s%c))) + eta)
    end function badd

    pure elemental function bsub(a, b) result(s)
        type(ball_t), intent(in) :: a, b
        type(ball_t) :: s

        s%c = a%c - b%c
        s%r = round_up(round_up(round_up(a%r + b%r) &
            + round_up(2.0_dp*u*cabs_hi(s%c))) + eta)
    end function bsub

    !> (A + a)(B + b) - A B = A b + a B + a b.
    pure elemental function bmul(a, b) result(p)
        type(ball_t), intent(in) :: a, b
        type(ball_t) :: p
        real(dp) :: ma, mb, prop, rnd

        p%c = a%c*b%c
        ma = cabs_hi(a%c)
        mb = cabs_hi(b%c)
        prop = round_up(round_up(round_up(ma*b%r) + round_up(a%r*mb)) &
            + round_up(a%r*b%r))
        rnd = round_up(3.0_dp*u*round_up(ma*mb))
        p%r = round_up(round_up(prop + rnd) + eta)
    end function bmul

    !> Product with an exact real x.
    pure elemental function bscale(a, x) result(p)
        type(ball_t), intent(in) :: a
        real(dp), intent(in) :: x
        type(ball_t) :: p

        p%c = a%c*x
        p%r = round_up(round_up(round_up(a%r*abs(x)) &
            + round_up(2.0_dp*u*cabs_hi(p%c))) + eta)
    end function bscale

    pure elemental function bmulr(a, x) result(p)
        type(ball_t), intent(in) :: a
        real(dp), intent(in) :: x
        type(ball_t) :: p

        p = bscale(a, x)
    end function bmulr

    !> The image of the disc under z -> 1/z is the disc with centre
    !> conj(c)/D and radius r/D, D = abs(c)^2 - r^2 > 0. D is enclosed in
    !> [dlo, dhi], 1/D in [ilo, ihi]; the computed centre conj(c) ihi differs
    !> from conj(c)/D by at most abs(c) (ihi - ilo) plus its own rounding
    !> u abs(c) ihi. A ball that may contain 0 gives an infinite radius.
    pure elemental function binv(a) result(q)
        type(ball_t), intent(in) :: a
        type(ball_t) :: q
        real(dp) :: x, y, dlo, dhi, ilo, ihi, ac, cr

        x = real(a%c, dp)
        y = aimag(a%c)
        dlo = round_down(round_down(round_down(x*x) + round_down(y*y)) &
            - round_up(a%r*a%r))
        dhi = round_up(round_up(round_up(x*x) + round_up(y*y)) &
            - round_down(a%r*a%r))
        if (.not. (dlo > 0.0_dp)) then
            q%c = (0.0_dp, 0.0_dp)
            q%r = huge(1.0_dp)
            return
        end if
        ilo = round_down(1.0_dp/dhi)
        ihi = round_up(1.0_dp/dlo)
        q%c = cmplx(x*ihi, -(y*ihi), dp)
        ac = cabs_hi(a%c)
        cr = round_up(ac*round_up(round_up(ihi - ilo) + round_up(u*ihi)))
        q%r = round_up(round_up(round_up(a%r*ihi) + cr) + eta)
    end function binv

    pure elemental function bdiv(a, b) result(q)
        type(ball_t), intent(in) :: a, b
        type(ball_t) :: q

        q = bmul(a, binv(b))
    end function bdiv

    !> Principal square root. For any z, abs(sqrt z - sqrt c) =
    !> abs(z - c)/abs(sqrt z + sqrt c) <= r/Re(sqrt c) because Re sqrt z >= 0.
    !> The computed centre w is certified a posteriori:
    !> abs(w - sqrt c) = abs(w^2 - c)/abs(w + sqrt c) <= abs(w^2 - c)/Re w.
    !> Centres on (-inf, 0] give an infinite radius.
    pure elemental function bsqrt(a) result(q)
        type(ball_t), intent(in) :: a
        type(ball_t) :: q
        type(ball_t) :: w2
        real(dp) :: x, y, m, s, t, rew, e, relo

        x = real(a%c, dp)
        y = aimag(a%c)
        if (x <= 0.0_dp .and. .not. (abs(y) > 0.0_dp)) then
            q%c = (0.0_dp, 0.0_dp)
            q%r = huge(1.0_dp)
            return
        end if
        m = abs(a%c)
        if (x >= 0.0_dp) then
            t = sqrt(0.5_dp*(m + x))
            s = y/(2.0_dp*t)
        else
            s = sign(sqrt(0.5_dp*(m - x)), y)
            t = abs(y)/(2.0_dp*abs(s))
        end if
        q%c = cmplx(t, s, dp)
        rew = t
        w2 = bsub(bmul(bpoint_c(q%c), bpoint_c(q%c)), bpoint_c(a%c))
        e = round_up(babs_hi(w2)/round_down(rew))
        relo = round_down(rew - e)
        if (.not. (rew > 0.0_dp) .or. .not. (relo > 0.0_dp)) then
            q%r = huge(1.0_dp)
            return
        end if
        q%r = round_up(round_up(round_up(a%r/relo) + e) + eta)
    end function bsqrt

    !> Integer power by binary exponentiation of balls; n < 0 inverts.
    pure elemental function bpowi(a, n) result(p)
        type(ball_t), intent(in) :: a
        integer, intent(in) :: n
        type(ball_t) :: p, base
        integer :: m

        p = bpoint_r(1.0_dp)
        base = a
        m = abs(n)
        do while (m > 0)
            if (modulo(m, 2) == 1) p = bmul(p, base)
            m = m/2
            if (m > 0) base = bmul(base, base)
        end do
        if (n < 0) p = binv(p)
    end function bpowi

    pure elemental function add_br(a, x) result(s)
        type(ball_t), intent(in) :: a
        real(dp), intent(in) :: x
        type(ball_t) :: s

        s = badd(a, bpoint_r(x))
    end function add_br

    pure elemental function add_rb(x, a) result(s)
        real(dp), intent(in) :: x
        type(ball_t), intent(in) :: a
        type(ball_t) :: s

        s = badd(bpoint_r(x), a)
    end function add_rb

    pure elemental function add_bz(a, z) result(s)
        type(ball_t), intent(in) :: a
        complex(dp), intent(in) :: z
        type(ball_t) :: s

        s = badd(a, bpoint_c(z))
    end function add_bz

    pure elemental function add_zb(z, a) result(s)
        complex(dp), intent(in) :: z
        type(ball_t), intent(in) :: a
        type(ball_t) :: s

        s = badd(bpoint_c(z), a)
    end function add_zb

    pure elemental function sub_br(a, x) result(s)
        type(ball_t), intent(in) :: a
        real(dp), intent(in) :: x
        type(ball_t) :: s

        s = bsub(a, bpoint_r(x))
    end function sub_br

    pure elemental function sub_rb(x, a) result(s)
        real(dp), intent(in) :: x
        type(ball_t), intent(in) :: a
        type(ball_t) :: s

        s = bsub(bpoint_r(x), a)
    end function sub_rb

    pure elemental function sub_bz(a, z) result(s)
        type(ball_t), intent(in) :: a
        complex(dp), intent(in) :: z
        type(ball_t) :: s

        s = bsub(a, bpoint_c(z))
    end function sub_bz

    pure elemental function sub_zb(z, a) result(s)
        complex(dp), intent(in) :: z
        type(ball_t), intent(in) :: a
        type(ball_t) :: s

        s = bsub(bpoint_c(z), a)
    end function sub_zb

    pure elemental function mul_rb(x, a) result(p)
        real(dp), intent(in) :: x
        type(ball_t), intent(in) :: a
        type(ball_t) :: p

        p = bscale(a, x)
    end function mul_rb

    pure elemental function mul_bz(a, z) result(p)
        type(ball_t), intent(in) :: a
        complex(dp), intent(in) :: z
        type(ball_t) :: p

        p = bmul(a, bpoint_c(z))
    end function mul_bz

    pure elemental function mul_zb(z, a) result(p)
        complex(dp), intent(in) :: z
        type(ball_t), intent(in) :: a
        type(ball_t) :: p

        p = bmul(bpoint_c(z), a)
    end function mul_zb

    pure elemental function div_br(a, x) result(q)
        type(ball_t), intent(in) :: a
        real(dp), intent(in) :: x
        type(ball_t) :: q

        q = bdiv(a, bpoint_r(x))
    end function div_br

    pure elemental function babs_hi(a) result(x)
        type(ball_t), intent(in) :: a
        real(dp) :: x

        x = round_up(cabs_hi(a%c) + a%r)
    end function babs_hi

    pure elemental function babs_lo(a) result(x)
        type(ball_t), intent(in) :: a
        real(dp) :: x

        x = max(round_down(cabs_lo(a%c) - a%r), 0.0_dp)
    end function babs_lo

    pure elemental function bre_lo(a) result(x)
        type(ball_t), intent(in) :: a
        real(dp) :: x

        x = round_down(real(a%c, dp) - a%r)
    end function bre_lo

    pure elemental function bre_hi(a) result(x)
        type(ball_t), intent(in) :: a
        real(dp) :: x

        x = round_up(real(a%c, dp) + a%r)
    end function bre_hi

    pure elemental function bim_lo(a) result(x)
        type(ball_t), intent(in) :: a
        real(dp) :: x

        x = round_down(aimag(a%c) - a%r)
    end function bim_lo

    pure elemental function bim_hi(a) result(x)
        type(ball_t), intent(in) :: a
        real(dp) :: x

        x = round_up(aimag(a%c) + a%r)
    end function bim_hi

    !> True when z certainly lies in the ball.
    pure elemental function bcontains(a, z) result(inside)
        type(ball_t), intent(in) :: a
        complex(dp), intent(in) :: z
        logical :: inside
        real(dp) :: d

        ! fl(z - c) = (z - c)(1 + delta) componentwise, abs(delta) <= u, and a
        ! subnormal difference is exact, so abs(z - c) <= d (1 + 4u) rounded.
        d = cabs_hi(z - a%c)
        inside = d*(1.0_dp + 4.0_dp*u) <= a%r
    end function bcontains

    !> True when 0 may lie in the ball.
    pure elemental function bcontains_zero(a) result(z)
        type(ball_t), intent(in) :: a
        logical :: z

        z = cabs_lo(a%c) <= a%r
    end function bcontains_zero

    !> A ball around the box centre that contains the box.
    pure elemental function ball_from_cinterval(z) result(b)
        type(cinterval_t), intent(in) :: z
        type(ball_t) :: b
        real(dp) :: x, y, hx, hy

        x = 0.5_dp*z%re%lo + 0.5_dp*z%re%hi
        y = 0.5_dp*z%im%lo + 0.5_dp*z%im%hi
        hx = max(round_up(x - z%re%lo), round_up(z%re%hi - x))
        hy = max(round_up(y - z%im%lo), round_up(z%im%hi - y))
        b%c = cmplx(x, y, dp)
        ! The relative margin 16u lets bcontains certify the box corners.
        b%r = cabs_hi(cmplx(hx, hy, dp))
        b%r = round_up(round_up(b%r*(1.0_dp + 16.0_dp*u)) + eta)
    end function ball_from_cinterval

    pure elemental function ball_from_interval(a) result(b)
        type(interval_t), intent(in) :: a
        type(ball_t) :: b
        real(dp) :: x

        x = 0.5_dp*a%lo + 0.5_dp*a%hi
        b%c = cmplx(x, 0.0_dp, dp)
        b%r = round_up(max(round_up(x - a%lo), round_up(a%hi - x)) + eta)
    end function ball_from_interval

    !> Bounding box of the disc.
    pure elemental function cinterval_from_ball(b) result(z)
        type(ball_t), intent(in) :: b
        type(cinterval_t) :: z

        z%re%lo = bre_lo(b)
        z%re%hi = bre_hi(b)
        z%im%lo = bim_lo(b)
        z%im%hi = bim_hi(b)
    end function cinterval_from_ball
end module fortnum_ball
