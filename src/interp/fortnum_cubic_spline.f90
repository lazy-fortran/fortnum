module fortnum_cubic_spline
    ! Interpolating cubic splines with not-a-knot end conditions on strictly
    ! increasing nodes, for several value columns sharing one node set.
    !
    ! The spline is stored by its nodal slopes s_i (Hermite form). Interior
    ! rows enforce continuity of the second derivative,
    !   h_i s_{i-1} + 2 (h_{i-1} + h_i) s_i + h_{i-1} s_{i+1}
    !     = 3 (h_i d_{i-1} + h_{i-1} d_i),   h_i = x_{i+1} - x_i,
    !   d_i = (y_{i+1} - y_i)/h_i,
    ! and the first and last rows continuity of the third derivative at x_2
    ! and x_{n-1} (de Boor, A Practical Guide to Splines, Ch. IV). The system
    ! is tridiagonal and solved once per column, so a cubic polynomial is
    ! reproduced exactly. Needs n >= 4 nodes.
    use fortnum_kinds, only: dp
    implicit none
    private

    public :: cubic_spline_t
    public :: cubic_spline_init
    public :: cubic_spline_eval

    type :: cubic_spline_t
        real(dp), allocatable :: x(:) ! nodes (n)
        real(dp), allocatable :: y(:, :) ! values (n, columns)
        real(dp), allocatable :: slope(:, :) ! nodal first derivatives
    end type cubic_spline_t

contains

    ! status: 0 ok, 1 fewer than four nodes, 2 nodes not strictly increasing,
    ! 3 shape mismatch between x and y.
    subroutine cubic_spline_init(x, y, sp, status)
        real(dp), intent(in) :: x(:), y(:, :)
        type(cubic_spline_t), intent(out) :: sp
        integer, intent(out) :: status
        real(dp), allocatable :: h(:), lower(:), diag(:), upper(:), d(:, :), rhs(:, :)
        real(dp) :: w
        integer :: n, i

        n = size(x)
        status = 0
        if (n < 4) then
            status = 1
            return
        end if
        if (size(y, 1) /= n) then
            status = 3
            return
        end if
        if (any(x(2:) <= x(:n - 1))) then
            status = 2
            return
        end if
        h = x(2:) - x(:n - 1)
        allocate (d(n - 1, size(y, 2)))
        do i = 1, n - 1
            d(i, :) = (y(i + 1, :) - y(i, :))/h(i)
        end do
        allocate (lower(n), diag(n), upper(n), rhs(n, size(y, 2)))
        lower = 0; upper = 0
        diag(1) = h(2)
        upper(1) = h(1) + h(2)
        rhs(1, :) = ((h(1) + 2*(h(1) + h(2)))*h(2)*d(1, :) + h(1)**2*d(2, :)) &
            /(h(1) + h(2))
        do i = 2, n - 1
            lower(i) = h(i)
            diag(i) = 2*(h(i - 1) + h(i))
            upper(i) = h(i - 1)
            rhs(i, :) = 3*(h(i)*d(i - 1, :) + h(i - 1)*d(i, :))
        end do
        lower(n) = h(n - 1) + h(n - 2)
        diag(n) = h(n - 2)
        rhs(n, :) = (h(n - 1)**2*d(n - 2, :) &
            + (2*(h(n - 2) + h(n - 1)) + h(n - 1))*h(n - 2)*d(n - 1, :)) &
            /(h(n - 2) + h(n - 1))
        ! Thomas elimination; all rows are nonsingular for increasing nodes.
        do i = 2, n
            w = lower(i)/diag(i - 1)
            diag(i) = diag(i) - w*upper(i - 1)
            rhs(i, :) = rhs(i, :) - w*rhs(i - 1, :)
        end do
        rhs(n, :) = rhs(n, :)/diag(n)
        do i = n - 1, 1, -1
            rhs(i, :) = (rhs(i, :) - upper(i)*rhs(i + 1, :))/diag(i)
        end do
        sp%x = x
        sp%y = y
        sp%slope = rhs
    end subroutine cubic_spline_init

    ! Values and first derivatives of every column at xq. Outside the nodes
    ! the end cubic is extrapolated.
    pure subroutine cubic_spline_eval(sp, xq, value, deriv)
        type(cubic_spline_t), intent(in) :: sp
        real(dp), intent(in) :: xq
        real(dp), intent(out) :: value(:), deriv(:)
        real(dp) :: h, t, h00, h10, h01, h11, d00, d10, d01, d11
        integer :: n, lo, hi, mid

        n = size(sp%x)
        lo = 1; hi = n - 1
        do while (hi > lo)
            mid = (lo + hi + 1)/2
            if (xq >= sp%x(mid)) then
                lo = mid
            else
                hi = mid - 1
            end if
        end do
        h = sp%x(lo + 1) - sp%x(lo)
        t = (xq - sp%x(lo))/h
        h00 = (1 + 2*t)*(1 - t)**2
        h10 = t*(1 - t)**2
        h01 = t*t*(3 - 2*t)
        h11 = t*t*(t - 1)
        d00 = 6*t*(t - 1)/h
        d10 = (1 - t)*(1 - 3*t)
        d01 = -d00
        d11 = t*(3*t - 2)
        value = h00*sp%y(lo, :) + h*h10*sp%slope(lo, :) &
            + h01*sp%y(lo + 1, :) + h*h11*sp%slope(lo + 1, :)
        deriv = d00*sp%y(lo, :) + d10*sp%slope(lo, :) &
            + d01*sp%y(lo + 1, :) + d11*sp%slope(lo + 1, :)
    end subroutine cubic_spline_eval

end module fortnum_cubic_spline
