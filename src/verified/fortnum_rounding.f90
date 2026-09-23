!> Directed rounding in round-to-nearest without changing the rounding mode.
!>
!> round_up(x) = fl(x + fl(phi |x| + eta)), phi = u (1 + 2u), eta = 2^-1074,
!> is succ(x) or succ(succ(x)) for every finite x (Rump, Zimmermann, Boldo,
!> Melquiond, "Computing predecessor and successor in rounding to nearest",
!> BIT 49 (2009) 419-431, Algorithm 2), hence round_up(fl(e)) >= e for the
!> exact value e of any correctly rounded IEEE operation. The proof needs only
!> fl(phi |x| + eta) > ulp(x)/2, which also holds when the compiler fuses the
!> inner multiply-add. An overflowed +inf is a valid upper bound and stays;
!> as a lower bound it is replaced by huge, the largest value the exact
!> result is known to exceed (and symmetrically for -inf). The formula needs
!> gradual underflow (no flush-to-zero) and no library call.
module fortnum_rounding
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    public :: unit_roundoff, round_up, round_down, add_up, add_down, sub_up
    public :: sub_down, mul_up, mul_down, div_up, div_down, sqrt_up, sqrt_down
    public :: sum_up, sum_down, gamma_up

    !> Unit roundoff u = 2^-53 of binary64 round-to-nearest.
    real(dp), parameter :: unit_roundoff = epsilon(1.0_dp)/2.0_dp
    real(dp), parameter :: phi = unit_roundoff*(1.0_dp + 2.0_dp*unit_roundoff)
    real(dp), parameter :: eta = tiny(1.0_dp)*epsilon(1.0_dp)

contains

    !> Upper bound: round_up(x) >= nearest(x, +1) for every finite x.
    pure elemental function round_up(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        y = max(x + (phi*min(abs(x), huge(x)) + eta), -huge(x))
    end function round_up

    !> Lower bound: round_down(x) <= nearest(x, -1) for every finite x.
    pure elemental function round_down(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        y = min(x - (phi*min(abs(x), huge(x)) + eta), huge(x))
    end function round_down

    pure elemental function add_up(a, b) result(y)
        real(dp), intent(in) :: a, b
        real(dp) :: y

        y = round_up(a + b)
    end function add_up

    pure elemental function add_down(a, b) result(y)
        real(dp), intent(in) :: a, b
        real(dp) :: y

        y = round_down(a + b)
    end function add_down

    pure elemental function sub_up(a, b) result(y)
        real(dp), intent(in) :: a, b
        real(dp) :: y

        y = round_up(a - b)
    end function sub_up

    pure elemental function sub_down(a, b) result(y)
        real(dp), intent(in) :: a, b
        real(dp) :: y

        y = round_down(a - b)
    end function sub_down

    pure elemental function mul_up(a, b) result(y)
        real(dp), intent(in) :: a, b
        real(dp) :: y

        y = round_up(a*b)
    end function mul_up

    pure elemental function mul_down(a, b) result(y)
        real(dp), intent(in) :: a, b
        real(dp) :: y

        y = round_down(a*b)
    end function mul_down

    pure elemental function div_up(a, b) result(y)
        real(dp), intent(in) :: a, b
        real(dp) :: y

        y = round_up(a/b)
    end function div_up

    pure elemental function div_down(a, b) result(y)
        real(dp), intent(in) :: a, b
        real(dp) :: y

        y = round_down(a/b)
    end function div_down

    !> Upper bound on sqrt(x) for x >= 0 (IEEE sqrt is correctly rounded).
    pure elemental function sqrt_up(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        y = round_up(sqrt(x))
    end function sqrt_up

    !> Lower bound on sqrt(x) for x >= 0, never negative.
    pure elemental function sqrt_down(x) result(y)
        real(dp), intent(in) :: x
        real(dp) :: y

        y = max(round_down(sqrt(x)), 0.0_dp)
    end function sqrt_down

    !> Upper bound on the exact sum of x, for any signs: every partial sum is
    !> rounded upward, so by induction it bounds the exact partial sum.
    pure function sum_up(x) result(s)
        real(dp), intent(in) :: x(:)
        real(dp) :: s
        integer :: i

        s = 0.0_dp
        do i = 1, size(x)
            s = round_up(s + x(i))
        end do
    end function sum_up

    !> Lower bound on the exact sum of x, for any signs.
    pure function sum_down(x) result(s)
        real(dp), intent(in) :: x(:)
        real(dp) :: s
        integer :: i

        s = 0.0_dp
        do i = 1, size(x)
            s = round_down(s + x(i))
        end do
    end function sum_down

    !> Upper bound on Higham's gamma_k = k u/(1 - k u); huge when k u >= 1/2.
    pure elemental function gamma_up(k) result(g)
        integer, intent(in) :: k
        real(dp) :: g, num, den

        num = mul_up(real(k, dp), unit_roundoff)
        if (num >= 0.5_dp) then
            g = huge(1.0_dp)
            return
        end if
        den = sub_down(1.0_dp, num)
        g = div_up(num, den)
    end function gamma_up
end module fortnum_rounding
