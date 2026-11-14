program test_fortnum_lapack_qr
    ! Householder QR against quad precision: Q orthogonal (unitary) and
    ! Q R = A, ORMQR in every side and transpose, column pivoting, and the
    ! least squares and minimum norm solutions of DGELS by their defining
    ! conditions.
    use fortnum_kinds, only: dp
    use fortnum_lapack_qr
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(qp), parameter :: eps = real(epsilon(1.0_dp), qp)

    call qr_real(90, 40)
    call qr_real(40, 40)
    call qr_real(7, 3)
    call qr_complex(70, 33)
    call apply_q(50, 20)
    call degenerate()
    call pivoted(60, 25)
    call least_squares(80, 30, 'N')
    call least_squares(30, 80, 'N')
    call least_squares(80, 30, 'T')
    call least_squares(30, 80, 'T')
    print '(a)', 'QR family matches quad precision'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    subroutine fill(x, seed)
        real(dp), intent(out) :: x(:, :)
        real(dp), intent(in) :: seed
        integer :: i, j
        do j = 1, size(x, 2)
            do i = 1, size(x, 1)
                x(i, j) = sin(seed*i + 1.91_dp*j)*exp(0.3_dp*cos(1.3_dp*i*j + seed))
            end do
        end do
    end subroutine fill

    function identity(n) result(e)
        integer, intent(in) :: n
        real(qp) :: e(n, n)
        integer :: i
        e = 0
        do i = 1, n
            e(i, i) = 1
        end do
    end function identity

    subroutine qr_real(m, n)
        integer, intent(in) :: m, n
        real(dp) :: a(m, n), f(m + 1, n), q(m, n), tau(min(m, n)), work(64*n)
        real(qp) :: r(n, n)
        integer :: info, k, i, j
        call fill(a, 0.37_dp)
        f = 5
        f(1:m, :) = a
        call fortnum_dgeqrf(m, n, f, m + 1, tau, work, size(work), info)
        call require(info == 0, 'dgeqrf info')
        k = min(m, n)
        r = 0
        do j = 1, n
            do i = 1, min(j, k)
                r(i, j) = f(i, j)
            end do
        end do
        q(:, 1:k) = f(1:m, 1:k)
        call fortnum_dorgqr(m, k, k, q, m, tau, work, size(work), info)
        call require(info == 0, 'dorgqr info')
        call require(maxval(abs(matmul(transpose(real(q(:, 1:k), qp)), real(q(:, 1:k), qp)) - identity(k))) <= &
            8*m*eps, 'Q orthogonal')
        call require(maxval(abs(matmul(real(q(:, 1:k), qp), r(1:k, :)) - real(a, qp))) <= &
            8*m*eps*maxval(abs(real(a, qp))), 'Q R = A')
    end subroutine qr_real

    subroutine qr_complex(m, n)
        integer, intent(in) :: m, n
        real(dp) :: re(m, n), im(m, n)
        complex(dp) :: a(m, n), f(m, n), q(m, n), tau(n), work(64*n), c(m, 2), c0(m, 2)
        complex(qp) :: r(n, n)
        integer :: info, i, j
        call fill(re, 0.3_dp)
        call fill(im, 0.8_dp)
        a = cmplx(re, im, dp)
        f = a
        call fortnum_zgeqrf(m, n, f, m, tau, work, size(work), info)
        r = 0
        do j = 1, n
            do i = 1, j
                r(i, j) = f(i, j)
            end do
            call require(aimag(f(j, j)) == 0, 'R has a real diagonal')
        end do
        q = f
        call fortnum_zungqr(m, n, n, q, m, tau, work, size(work), info)
        call require(info == 0 .and. maxval(abs(matmul(conjg(transpose(cmplx(q, kind=qp))), cmplx(q, kind=qp)) &
            - identity(n))) <= 16*m*eps, 'Q unitary')
        call require(maxval(abs(matmul(cmplx(q, kind=qp), r) - cmplx(a, kind=qp))) <= &
            16*m*eps*maxval(abs(cmplx(a, kind=qp))), 'complex Q R = A')
        c0(:, 1) = cmplx(re(:, 1), -im(:, 2), dp)
        c0(:, 2) = cmplx(im(:, 3), re(:, 4), dp)
        c = c0
        call fortnum_zunmqr('L', 'C', m, 2, n, f, m, tau, c, m, work, size(work), info)
        call fortnum_zunmqr('L', 'N', m, 2, n, f, m, tau, c, m, work, size(work), info)
        call require(maxval(abs(c - c0)) <= 64*m*epsilon(1.0_dp)*maxval(abs(c0)), 'Q Q^H c = c')
    end subroutine qr_complex

    ! every side and transpose of DORMQR against the explicit Q
    subroutine apply_q(m, k)
        integer, intent(in) :: m, k
        real(dp) :: a(m, k), f(m, k), tau(k), q(m, m), cl(m, 7), cr(7, m), work(64*m)
        real(qp) :: qq(m, m)
        integer :: info, j, s, t
        character :: side, tr
        call fill(a, 0.61_dp)
        f = a
        call fortnum_dgeqrf(m, k, f, m, tau, work, size(work), info)
        q = 0
        do j = 1, m
            q(j, j) = 1
        end do
        call fortnum_dormqr('L', 'N', m, m, k, f, m, tau, q, m, work, size(work), info)
        qq = real(q, qp)
        do s = 1, 2
            do t = 1, 2
                side = 'LR'(s:s)
                tr = 'NT'(t:t)
                if (side == 'L') then
                    call fill(cl, 0.2_dp)
                    block
                        real(dp) :: c0(m, 7)
                        c0 = cl
                        call fortnum_dormqr('L', tr, m, 7, k, f, m, tau, cl, m, work, size(work), info)
                        call require(maxval(abs(real(cl, qp) - matmul(merge(qq, transpose(qq), tr == 'N'), real(c0, qp)))) <= 64*m*eps, &
                            'dormqr L '//tr)
                    end block
                else
                    block
                        real(dp) :: c0(7, m)
                        call fill(cr, 0.9_dp)
                        c0 = cr
                        call fortnum_dormqr('R', tr, 7, m, k, f, m, tau, cr, 7, work, size(work), info)
                        call require(maxval(abs(real(cr, qp) - matmul(real(c0, qp), merge(qq, transpose(qq), tr == 'N')))) <= 64*m*eps, &
                            'dormqr R '//tr)
                    end block
                end if
            end do
        end do
    end subroutine apply_q

    subroutine degenerate()
        real(dp) :: alpha, x(3), tau
        alpha = -2
        x = 0
        call fortnum_dlarfg(4, alpha, x, 1, tau)
        call require(tau == 0 .and. alpha == -2, 'dlarfg leaves a vector with zero tail alone')
        alpha = 3
        x = [4.0_dp, 0.0_dp, 0.0_dp]
        call fortnum_dlarfg(4, alpha, x, 1, tau)
        call require(abs(alpha + 5) < 1e-14_dp .and. abs(tau - 1.6_dp) < 1e-14_dp, 'dlarfg beta and tau')
    end subroutine degenerate

    subroutine pivoted(m, n)
        integer, intent(in) :: m, n
        real(dp) :: a(m, n), f(m, n), q(m, n), tau(n), work(4*n + 1)
        real(qp) :: r(n, n), ap(m, n)
        integer :: jpvt(n), info, i, j
        call fill(a, 0.45_dp)
        a(:, 3) = 1e-3_dp*a(:, 3)
        f = a
        jpvt = 0
        jpvt(5) = 1
        call fortnum_dgeqp3(m, n, f, m, jpvt, tau, work, size(work), info)
        call require(info == 0 .and. jpvt(1) == 5, 'dgeqp3 keeps the leading column first')
        do j = 3, n
            call require(abs(f(j, j)) <= abs(f(j - 1, j - 1))*(1 + 1e-10_dp), 'pivoted |R(i,i)| decrease')
        end do
        r = 0
        do j = 1, n
            do i = 1, j
                r(i, j) = f(i, j)
            end do
            ap(:, j) = a(:, jpvt(j))
        end do
        q = f
        call fortnum_dorgqr(m, n, n, q, m, tau, work, size(work), info)
        call require(maxval(abs(matmul(real(q, qp), r) - ap)) <= 16*m*eps*maxval(abs(real(a, qp))), 'A P = Q R')
    end subroutine pivoted

    ! overdetermined: A^T (A x - b) = 0; underdetermined: A x = b and x in range(A^T)
    subroutine least_squares(m, n, tr)
        integer, intent(in) :: m, n
        character, intent(in) :: tr
        real(dp) :: a(m, n), f(m, n), b(max(m, n), 2), b0(max(m, n), 2), work(64*max(m, n))
        real(qp) :: op(merge(m, n, tr == 'N'), merge(n, m, tr == 'N')), x(merge(n, m, tr == 'N'), 2), r(merge(m, n, tr == 'N'), 2)
        integer :: info, rows, cols
        call fill(a, 0.53_dp)
        call fill(b0, 0.21_dp)
        if (tr == 'N') then
            op = real(a, qp)
        else
            op = transpose(real(a, qp))
        end if
        rows = size(op, 1)
        cols = size(op, 2)
        f = a
        b = b0
        call fortnum_dgels(tr, m, n, 2, f, m, b, max(m, n), work, size(work), info)
        call require(info == 0, 'dgels info')
        x = real(b(1:cols, :), qp)
        r = matmul(op, x) - real(b0(1:rows, :), qp)
        if (rows >= cols) then
            call require(maxval(abs(matmul(transpose(op), r))) <= &
                1e-11_qp*maxval(abs(op))**2*maxval(abs(x)), 'dgels normal equations '//tr)
        else
            call require(maxval(abs(r)) <= 64*rows*eps*maxval(abs(op))*maxval(abs(x)), 'dgels solves '//tr)
            ! minimum norm: x is in the row space, so x = op^T y for some y
            call require(in_row_space(op, x), 'dgels minimum norm '//tr)
        end if
    end subroutine least_squares

    logical function in_row_space(op, x)
        real(qp), intent(in) :: op(:, :), x(:, :)
        real(qp) :: g(size(op, 1), size(op, 1)), y(size(op, 1), size(x, 2)), resid
        ! y from the normal equations (op op^T) y = op x, then x - op^T y
        g = matmul(op, transpose(op))
        y = matmul(op, x)
        call gauss(g, y)
        resid = maxval(abs(x - matmul(transpose(op), y)))
        in_row_space = resid <= 1e-10_qp*maxval(abs(x))
    end function in_row_space

    subroutine gauss(g, y)
        real(qp), intent(inout) :: g(:, :), y(:, :)
        integer :: k, i, n
        n = size(g, 1)
        do k = 1, n
            do i = k + 1, n
                y(i, :) = y(i, :) - g(i, k)/g(k, k)*y(k, :)
                g(i, k:) = g(i, k:) - g(i, k)/g(k, k)*g(k, k:)
            end do
        end do
        do k = n, 1, -1
            y(k, :) = (y(k, :) - matmul(g(k, k + 1:n), y(k + 1:n, :)))/g(k, k)
        end do
    end subroutine gauss

end program test_fortnum_lapack_qr
