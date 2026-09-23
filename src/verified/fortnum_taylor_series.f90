!> Order-by-order Taylor-series arithmetic.
!>
!> A series is an array a(0:q) of coefficient enclosures; each routine
!> computes only coefficient k of its result from coefficients 0..k of its
!> arguments (or 0..k-1 for the transcendental recurrences), so evaluating an
!> ODE right-hand side as a sequence of these calls, for k = 0, 1, ..., q,
!> yields rigorous Taylor coefficients of the composed expression with O(q^2)
!> total work per elementary operation instead of symbolic differentiation.
!> This is the standard automatic Taylor-coefficient recurrence (Moore 1966;
!> Corliss and Chang 1982), generalized here over three coefficient rings so
!> the same driver code serves plain interval enclosures, first-order real
!> forward AD (fortnum_idual, for a real Jacobian alongside the flow), and
!> first-order complex forward AD (fortnum_cdual, for a holomorphic
!> variational equation). `fortnum_validated_ode` composes these to build the
!> high-order Taylor-Lohner step.
!>
!> Recurrences, coefficient index k, ring operations + - * / and the same
!> ring's sin/cos/exp/sqrt at index 0:
!>   mul:   c_k = sum_{j=0}^k a_j b_{k-j}
!>   div:   c_k = (a_k - sum_{j=1}^k b_j c_{k-j}) / b_0
!>   sincos: s_0 = sin a_0, c_0 = cos a_0;
!>           k s_k =  sum_{j=1}^k j a_j c_{k-j},
!>           k c_k = -sum_{j=1}^k j a_j s_{k-j}
!>   exp:   e_0 = exp a_0;  k e_k = sum_{j=1}^k j a_j e_{k-j}
!>   sqrt:  s_0 = sqrt a_0; 2 s_0 s_k = a_k - sum_{j=1}^{k-1} s_j s_{k-j}
!>   lin:   c_k = alpha a_k + beta b_k  (alpha, beta real)
module fortnum_taylor_series
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, interval, cinterval, operator(+), &
        operator(-), operator(*), operator(/), sin, cos, exp, sqrt
    use fortnum_idual, only: idual_t, idual_const, operator(+), &
        operator(-), operator(*), operator(/), sin, cos, exp, sqrt
    use fortnum_cdual, only: cdual_t, cdual_const, operator(+), &
        operator(-), operator(*), operator(/), sin, cos, exp, sqrt
    implicit none
    private

    public :: ts_zero, ts_mul, ts_div, ts_sincos, ts_exp, ts_sqrt, ts_lin

    interface ts_zero
        module procedure ts_zero_i, ts_zero_d, ts_zero_c
    end interface ts_zero

    interface ts_mul
        module procedure ts_mul_i, ts_mul_d, ts_mul_c
    end interface ts_mul

    interface ts_div
        module procedure ts_div_i, ts_div_d, ts_div_c
    end interface ts_div

    interface ts_sincos
        module procedure ts_sincos_i, ts_sincos_d, ts_sincos_c
    end interface ts_sincos

    interface ts_exp
        module procedure ts_exp_i, ts_exp_d, ts_exp_c
    end interface ts_exp

    interface ts_sqrt
        module procedure ts_sqrt_i, ts_sqrt_d, ts_sqrt_c
    end interface ts_sqrt

    interface ts_lin
        module procedure ts_lin_i, ts_lin_d, ts_lin_c
    end interface ts_lin

contains

    ! ---- interval_t coefficients ----

    pure subroutine ts_zero_i(k, c)
        integer, intent(in) :: k
        type(interval_t), intent(inout) :: c(0:)
        c(k) = interval(0.0_dp)
    end subroutine ts_zero_i

    pure subroutine ts_mul_i(k, a, b, c)
        integer, intent(in) :: k
        type(interval_t), intent(in) :: a(0:), b(0:)
        type(interval_t), intent(inout) :: c(0:)
        integer :: j
        c(k) = a(0)*b(k)
        do j = 1, k
            c(k) = c(k) + a(j)*b(k - j)
        end do
    end subroutine ts_mul_i

    pure subroutine ts_div_i(k, a, b, c)
        integer, intent(in) :: k
        type(interval_t), intent(in) :: a(0:), b(0:)
        type(interval_t), intent(inout) :: c(0:)
        type(interval_t) :: acc
        integer :: j
        acc = a(k)
        do j = 1, k
            acc = acc - b(j)*c(k - j)
        end do
        c(k) = acc/b(0)
    end subroutine ts_div_i

    pure subroutine ts_sincos_i(k, a, s, c)
        integer, intent(in) :: k
        type(interval_t), intent(in) :: a(0:)
        type(interval_t), intent(inout) :: s(0:), c(0:)
        type(interval_t) :: as_, ac
        integer :: j
        if (k == 0) then
            s(0) = sin(a(0))
            c(0) = cos(a(0))
            return
        end if
        as_ = real(1, dp)*(a(1)*c(k - 1))
        ac = real(1, dp)*(a(1)*s(k - 1))
        do j = 2, k
            as_ = as_ + real(j, dp)*(a(j)*c(k - j))
            ac = ac + real(j, dp)*(a(j)*s(k - j))
        end do
        s(k) = as_/real(k, dp)
        c(k) = -(ac/real(k, dp))
    end subroutine ts_sincos_i

    pure subroutine ts_exp_i(k, a, e)
        integer, intent(in) :: k
        type(interval_t), intent(in) :: a(0:)
        type(interval_t), intent(inout) :: e(0:)
        type(interval_t) :: acc
        integer :: j
        if (k == 0) then
            e(0) = exp(a(0))
            return
        end if
        acc = real(1, dp)*(a(1)*e(k - 1))
        do j = 2, k
            acc = acc + real(j, dp)*(a(j)*e(k - j))
        end do
        e(k) = acc/real(k, dp)
    end subroutine ts_exp_i

    pure subroutine ts_sqrt_i(k, a, s)
        integer, intent(in) :: k
        type(interval_t), intent(in) :: a(0:)
        type(interval_t), intent(inout) :: s(0:)
        type(interval_t) :: acc
        integer :: j
        if (k == 0) then
            s(0) = sqrt(a(0))
            return
        end if
        acc = a(k)
        do j = 1, k - 1
            acc = acc - s(j)*s(k - j)
        end do
        s(k) = acc/(real(2, dp)*s(0))
    end subroutine ts_sqrt_i

    pure subroutine ts_lin_i(k, alpha, a, beta, b, c)
        integer, intent(in) :: k
        real(dp), intent(in) :: alpha, beta
        type(interval_t), intent(in) :: a(0:), b(0:)
        type(interval_t), intent(inout) :: c(0:)
        c(k) = alpha*a(k) + beta*b(k)
    end subroutine ts_lin_i

    ! ---- idual_t coefficients (interval value plus real gradient) ----

    pure subroutine ts_zero_d(k, n, c)
        integer, intent(in) :: k, n
        type(idual_t), intent(inout) :: c(0:)
        c(k) = idual_const(interval(0.0_dp), n)
    end subroutine ts_zero_d

    pure subroutine ts_mul_d(k, a, b, c)
        integer, intent(in) :: k
        type(idual_t), intent(in) :: a(0:), b(0:)
        type(idual_t), intent(inout) :: c(0:)
        integer :: j
        c(k) = a(0)*b(k)
        do j = 1, k
            c(k) = c(k) + a(j)*b(k - j)
        end do
    end subroutine ts_mul_d

    pure subroutine ts_div_d(k, a, b, c)
        integer, intent(in) :: k
        type(idual_t), intent(in) :: a(0:), b(0:)
        type(idual_t), intent(inout) :: c(0:)
        type(idual_t) :: acc
        integer :: j
        acc = a(k)
        do j = 1, k
            acc = acc - b(j)*c(k - j)
        end do
        c(k) = acc/b(0)
    end subroutine ts_div_d

    pure subroutine ts_sincos_d(k, a, s, c)
        integer, intent(in) :: k
        type(idual_t), intent(in) :: a(0:)
        type(idual_t), intent(inout) :: s(0:), c(0:)
        type(idual_t) :: as_, ac
        integer :: j
        if (k == 0) then
            s(0) = sin(a(0))
            c(0) = cos(a(0))
            return
        end if
        as_ = real(1, dp)*(a(1)*c(k - 1))
        ac = real(1, dp)*(a(1)*s(k - 1))
        do j = 2, k
            as_ = as_ + real(j, dp)*(a(j)*c(k - j))
            ac = ac + real(j, dp)*(a(j)*s(k - j))
        end do
        s(k) = as_/real(k, dp)
        c(k) = -(ac/real(k, dp))
    end subroutine ts_sincos_d

    pure subroutine ts_exp_d(k, a, e)
        integer, intent(in) :: k
        type(idual_t), intent(in) :: a(0:)
        type(idual_t), intent(inout) :: e(0:)
        type(idual_t) :: acc
        integer :: j
        if (k == 0) then
            e(0) = exp(a(0))
            return
        end if
        acc = real(1, dp)*(a(1)*e(k - 1))
        do j = 2, k
            acc = acc + real(j, dp)*(a(j)*e(k - j))
        end do
        e(k) = acc/real(k, dp)
    end subroutine ts_exp_d

    pure subroutine ts_sqrt_d(k, a, s)
        integer, intent(in) :: k
        type(idual_t), intent(in) :: a(0:)
        type(idual_t), intent(inout) :: s(0:)
        type(idual_t) :: acc
        integer :: j
        if (k == 0) then
            s(0) = sqrt(a(0))
            return
        end if
        acc = a(k)
        do j = 1, k - 1
            acc = acc - s(j)*s(k - j)
        end do
        s(k) = acc/(real(2, dp)*s(0))
    end subroutine ts_sqrt_d

    pure subroutine ts_lin_d(k, alpha, a, beta, b, c)
        integer, intent(in) :: k
        real(dp), intent(in) :: alpha, beta
        type(idual_t), intent(in) :: a(0:), b(0:)
        type(idual_t), intent(inout) :: c(0:)
        c(k) = alpha*a(k) + beta*b(k)
    end subroutine ts_lin_d

    ! ---- cdual_t coefficients (complex-box value, holomorphic derivative) --

    pure subroutine ts_zero_c(k, n, c)
        integer, intent(in) :: k, n
        type(cdual_t), intent(inout) :: c(0:)
        c(k) = cdual_const(cinterval(0.0_dp, 0.0_dp), n)
    end subroutine ts_zero_c

    pure subroutine ts_mul_c(k, a, b, c)
        integer, intent(in) :: k
        type(cdual_t), intent(in) :: a(0:), b(0:)
        type(cdual_t), intent(inout) :: c(0:)
        integer :: j
        c(k) = a(0)*b(k)
        do j = 1, k
            c(k) = c(k) + a(j)*b(k - j)
        end do
    end subroutine ts_mul_c

    pure subroutine ts_div_c(k, a, b, c)
        integer, intent(in) :: k
        type(cdual_t), intent(in) :: a(0:), b(0:)
        type(cdual_t), intent(inout) :: c(0:)
        type(cdual_t) :: acc
        integer :: j
        acc = a(k)
        do j = 1, k
            acc = acc - b(j)*c(k - j)
        end do
        c(k) = acc/b(0)
    end subroutine ts_div_c

    pure subroutine ts_sincos_c(k, a, s, c)
        integer, intent(in) :: k
        type(cdual_t), intent(in) :: a(0:)
        type(cdual_t), intent(inout) :: s(0:), c(0:)
        type(cdual_t) :: as_, ac
        integer :: j
        if (k == 0) then
            s(0) = sin(a(0))
            c(0) = cos(a(0))
            return
        end if
        as_ = real(1, dp)*(a(1)*c(k - 1))
        ac = real(1, dp)*(a(1)*s(k - 1))
        do j = 2, k
            as_ = as_ + real(j, dp)*(a(j)*c(k - j))
            ac = ac + real(j, dp)*(a(j)*s(k - j))
        end do
        s(k) = as_/real(k, dp)
        c(k) = -(ac/real(k, dp))
    end subroutine ts_sincos_c

    pure subroutine ts_exp_c(k, a, e)
        integer, intent(in) :: k
        type(cdual_t), intent(in) :: a(0:)
        type(cdual_t), intent(inout) :: e(0:)
        type(cdual_t) :: acc
        integer :: j
        if (k == 0) then
            e(0) = exp(a(0))
            return
        end if
        acc = real(1, dp)*(a(1)*e(k - 1))
        do j = 2, k
            acc = acc + real(j, dp)*(a(j)*e(k - j))
        end do
        e(k) = acc/real(k, dp)
    end subroutine ts_exp_c

    pure subroutine ts_sqrt_c(k, a, s)
        integer, intent(in) :: k
        type(cdual_t), intent(in) :: a(0:)
        type(cdual_t), intent(inout) :: s(0:)
        type(cdual_t) :: acc
        integer :: j
        if (k == 0) then
            s(0) = sqrt(a(0))
            return
        end if
        acc = a(k)
        do j = 1, k - 1
            acc = acc - s(j)*s(k - j)
        end do
        s(k) = acc/(real(2, dp)*s(0))
    end subroutine ts_sqrt_c

    pure subroutine ts_lin_c(k, alpha, a, beta, b, c)
        integer, intent(in) :: k
        real(dp), intent(in) :: alpha, beta
        type(cdual_t), intent(in) :: a(0:), b(0:)
        type(cdual_t), intent(inout) :: c(0:)
        c(k) = alpha*a(k) + beta*b(k)
    end subroutine ts_lin_c

end module fortnum_taylor_series
