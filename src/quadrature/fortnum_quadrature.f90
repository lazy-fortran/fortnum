module fortnum_quadrature
    ! Fixed-rule Gauss quadrature: node/weight generation and interval-mapped
    ! integration helpers.
    !
    ! DERIVATIVE POLICY (ad.md §1, §4): transparent w.r.t. integrand values.
    !   Fixed-rule quadrature evaluates I = sum_i w_i * f(x_i).  The map
    !   f -> I is linear with Jacobian J = w^T (a 1 x n row vector).
    !   Active: integrand values f(x_i) supplied by the caller.
    !   Inactive: n, a, b, nodes x_i, weights w_i (rule parameters).
    !
    !   Derivative products (ad.md §2):
    !     gauss_legendre_jvp  -- forward:  dI  = sum_i w_i v_i   (scalar result)
    !     gauss_legendre_vjp  -- reverse:  df_i = u * w_i         (u scalar costate)
    !     gauss_legendre_grad -- gradient: dI/df_i = w_i (u = 1 case of vjp)
    !
    !   HVP: deferred; I is linear in f so the Hessian is zero -- HVP would
    !   return zero and carries no information.  Not implemented.
    !
    ! Algorithms: Golub-Welsch (1969) via Newton iteration on Legendre recurrence.
    !   G. H. Golub, J. H. Welsch, Math. Comp. 23 (1969) 221-230.
    !   DLMF 18.9.8 (three-term recurrence), DLMF 18.9.17 (derivative),
    !   A&S 22.16.6 (asymptotic initial estimate), A&S 25.4.29 (weight formula).
    use fortnum_kinds, only: dp
    implicit none
    private

    public :: gauss_legendre
    public :: gauss_legendre_ab
    public :: gauss_legendre_jvp
    public :: gauss_legendre_vjp
    public :: gauss_legendre_grad

contains

    ! Gauss-Legendre nodes and weights on [-1, 1].
    ! x(n) receives nodes in ascending order, w(n) the corresponding weights.
    ! Nodes are roots of P_n converged by Newton iteration from the asymptotic
    ! estimate cos(pi*(i - 1/4)/(n + 1/2)) (A&S 22.16.6); weights from the
    ! closed form 2/((1 - x^2)*P_n'(x)^2) (A&S 25.4.29).  Symmetry is
    ! exploited: only ceil(n/2) Newton solves are performed.
    pure subroutine gauss_legendre(n, x, w)
        integer,  intent(in)  :: n
        real(dp), intent(out) :: x(n), w(n)

        real(dp), parameter :: pi = 3.14159265358979324_dp
        integer,  parameter :: max_iter = 16

        integer  :: m, i, it
        real(dp) :: xh((n + 1)/2), p((n + 1)/2), dpdx((n + 1)/2)
        real(dp) :: dx((n + 1)/2)

        if (n < 1) error stop "fortnum_quadrature: n must be >= 1"
        m = (n + 1)/2
        do i = 1, m
            xh(i) = cos(pi*(real(i, dp) - 0.25_dp)/(real(n, dp) + 0.5_dp))
        end do
        do it = 1, max_iter
            call legendre_eval(n, m, xh, p, dpdx)
            dx  = p/dpdx
            xh  = xh - dx
            if (maxval(abs(dx)) <= 4.0_dp*epsilon(1.0_dp)) exit
        end do
        call legendre_eval(n, m, xh, p, dpdx)
        do i = 1, m
            x(n + 1 - i) = xh(i)
            x(i)         = -xh(i)
            w(i)         = 2.0_dp/((1.0_dp - xh(i)*xh(i))*dpdx(i)*dpdx(i))
            w(n + 1 - i) = w(i)
        end do
    end subroutine gauss_legendre

    ! Three-term recurrence P_k (DLMF 18.9.8) and its derivative (DLMF 18.9.17),
    ! evaluated simultaneously for all m nodes to allow vectorization.
    pure subroutine legendre_eval(n, m, x, p, dpdx)
        integer,  intent(in)  :: n, m
        real(dp), intent(in)  :: x(m)
        real(dp), intent(out) :: p(m), dpdx(m)

        real(dp) :: p_prev(m), p_next(m), c_cur, c_prev
        integer  :: k

        p_prev = 1.0_dp
        p      = x
        do k = 1, n - 1
            c_cur  = real(2*k + 1, dp)/real(k + 1, dp)
            c_prev = real(k,       dp)/real(k + 1, dp)
            p_next = c_cur*x*p - c_prev*p_prev
            p_prev = p
            p      = p_next
        end do
        ! DLMF 18.9.17: P_n'(x) = n*(x*P_n(x) - P_{n-1}(x))/(x^2 - 1)
        dpdx = real(n, dp)*(x*p - p_prev)/(x*x - 1.0_dp)
    end subroutine legendre_eval

    ! Gauss-Legendre nodes and weights mapped to [a, b].
    ! x(n) and w(n) are the transformed rule; the quadrature sum
    !   sum_i w(i)*f(x(i))
    ! approximates integral_a^b f(t) dt.
    subroutine gauss_legendre_ab(n, a, b, x, w)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: a, b
        real(dp), intent(out) :: x(n), w(n)

        real(dp) :: half_length, midpoint

        call gauss_legendre(n, x, w)
        half_length = 0.5_dp*(b - a)
        midpoint    = 0.5_dp*(a + b)
        x = midpoint + half_length*x
        w = half_length*w
    end subroutine gauss_legendre_ab

    ! Forward product for the linear map f -> I = sum_i w_i f_i.
    ! Active inputs: f (sampled integrand values at the rule nodes), tangent v.
    ! Active output: jv = dI/df . v = sum_i w_i v_i  (scalar; returned in jv(1)).
    ! Inactive: w (the rule weights, produced by gauss_legendre or gauss_legendre_ab).
    ! w and v must have the same length n; jv must have size >= 1.
    pure subroutine gauss_legendre_jvp(w, v, jv)
        real(dp), intent(in)  :: w(:)   ! quadrature weights (inactive rule parameter)
        real(dp), intent(in)  :: v(:)   ! tangent for integrand values, size n
        real(dp), intent(out) :: jv(1)  ! forward product dI, scalar
        jv(1) = dot_product(w, v)
    end subroutine gauss_legendre_jvp

    ! Reverse product for the linear map f -> I = sum_i w_i f_i.
    ! Active inputs: f (sampled integrand values, inactive here -- only w needed),
    !   u (scalar output costate, size 1).
    ! Active output: jtu = (dI/df)^T u, i.e. jtu_i = u(1) * w_i.
    ! Inactive: w (rule weights).
    pure subroutine gauss_legendre_vjp(w, u, jtu)
        real(dp), intent(in)  :: w(:)   ! quadrature weights (inactive rule parameter)
        real(dp), intent(in)  :: u(1)   ! output costate (scalar), size 1
        real(dp), intent(out) :: jtu(:) ! reverse product, size n
        jtu = u(1) * w
    end subroutine gauss_legendre_vjp

    ! Gradient of I = sum_i w_i f_i w.r.t. the integrand sample values f_i.
    ! dI/df_i = w_i; this is the u = 1 case of gauss_legendre_vjp and exposes
    ! the weights directly as the gradient vector.
    ! Inactive: n (determined from size(w)); w (rule weights, already computed).
    ! Output: grad(n) receives the weight vector w.
    pure subroutine gauss_legendre_grad(w, grad)
        real(dp), intent(in)  :: w(:)    ! quadrature weights (inactive rule parameter)
        real(dp), intent(out) :: grad(:) ! dI/df_i = w_i, size n
        grad = w
    end subroutine gauss_legendre_grad

end module fortnum_quadrature
