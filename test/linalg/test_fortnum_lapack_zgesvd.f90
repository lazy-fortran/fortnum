program test_fortnum_lapack_zgesvd
    ! Complex SVD by its defining properties: A = U S V^H with orthonormal
    ! factors and descending S for tall, wide and square shapes and every
    ! job, singular values known by construction from unitary factors, a
    ! rank-deficient matrix and argument errors.
    use fortnum_kinds, only: dp
    use fortnum_lapack_svd, only: fortnum_zgesvd
    implicit none
    real(dp), parameter :: eps = epsilon(1.0_dp)
    integer, parameter :: shapes(2, 7) = reshape([1, 1, 5, 3, 3, 5, 40, 25, 25, 40, 60, 60, 7, 70], [2, 7])
    integer :: k

    do k = 1, size(shapes, 2)
        call jobs(shapes(1, k), shapes(2, k))
    end do
    call known_values(30, 20)
    call known_values(20, 30)
    call rank_deficient(30, 20, 4)
    call argument_errors()
    print '(a)', 'complex SVD factors reproduce A'

contains

    subroutine jobs(m, n)
        integer, intent(in) :: m, n
        complex(dp) :: a(m, n), b(m, n), u(m, m), vt(n, n), us(m, min(m, n)), vts(min(m, n), n)
        complex(dp) :: d(1, 1), e(1, 1), q(1)
        complex(dp), allocatable :: work(:)
        real(dp) :: s(min(m, n)), s2(min(m, n)), rwork(5*min(m, n))
        integer :: info, mn
        mn = min(m, n)
        a = general(m, n, 0.37_dp + 0.01_dp*m)
        b = a
        call fortnum_zgesvd('A', 'A', m, n, b, m, s, u, m, vt, n, q, -1, rwork, info)
        call require(info == 0, 'query')
        allocate (work(int(real(q(1)))))
        call fortnum_zgesvd('A', 'A', m, n, b, m, s, u, m, vt, n, work, size(work), rwork, info)
        call require(info == 0, 'A A info')
        call require(unit_error(u) <= 20*max(m, n)*eps, 'A A: full U unitary')
        call require(unit_error(conjg(transpose(vt))) <= 20*max(m, n)*eps, 'A A: full V unitary')
        call check(a, s, u(:, 1:mn), vt(1:mn, :), 'A A')
        b = a
        call fortnum_zgesvd('S', 'S', m, n, b, m, s2, us, m, vts, mn, work, size(work), rwork, info)
        call require(info == 0, 'S S info')
        call check(a, s2, us, vts, 'S S')
        b = a
        call fortnum_zgesvd('O', 'S', m, n, b, m, s2, d, 1, vts, mn, work, size(work), rwork, info)
        call require(info == 0, 'O S info')
        call check(a, s2, b(:, 1:mn), vts, 'O S')
        b = a
        call fortnum_zgesvd('S', 'O', m, n, b, m, s2, us, m, d, 1, work, size(work), rwork, info)
        call require(info == 0, 'S O info')
        call check(a, s2, us, b(1:mn, :), 'S O')
        b = a
        call fortnum_zgesvd('N', 'N', m, n, b, m, s2, d, 1, e, 1, work, size(work), rwork, info)
        call require(info == 0, 'N N info')
        call require(maxval(abs(s2 - s)) <= 20*max(m, n)*eps*s(1), 'values alone')
    end subroutine jobs

    ! A = Q1 diag(sigma) Q2^H with unitary reflector products; sigma graded
    ! over eight orders of magnitude
    subroutine known_values(m, n)
        integer, intent(in) :: m, n
        complex(dp) :: q1(m, m), q2(n, n), x(m, n), a(m, n), u(1, 1), vt(1, 1), work(4*(m + n))
        real(dp) :: sig(min(m, n)), s(min(m, n)), rwork(5*min(m, n))
        integer :: i, info
        q1 = unitary(m, 0.3_dp)
        q2 = unitary(n, 0.7_dp)
        x = 0
        do i = 1, min(m, n)
            sig(i) = 10.0_dp**(-8*real(i - 1, dp)/(min(m, n) - 1))
            x(i, i) = sig(i)
        end do
        a = matmul(matmul(q1, x), conjg(transpose(q2)))
        call fortnum_zgesvd('N', 'N', m, n, a, m, s, u, 1, vt, 1, work, size(work), rwork, info)
        call require(info == 0, 'known values info')
        call require(maxval(abs(s - sig)) <= 50*max(m, n)*eps, 'singular values known by construction')
    end subroutine known_values

    subroutine rank_deficient(m, n, r)
        integer, intent(in) :: m, n, r
        complex(dp) :: a(m, n), b(m, n), u(m, m), vt(n, n), work(4*(m + n))
        real(dp) :: s(min(m, n)), rwork(5*min(m, n))
        integer :: info
        a = matmul(general(m, r, 0.21_dp), general(r, n, 0.53_dp))
        b = a
        call fortnum_zgesvd('A', 'A', m, n, b, m, s, u, m, vt, n, work, size(work), rwork, info)
        call require(info == 0, 'rank deficient info')
        call require(all(s(r + 1:) <= 50*max(m, n)*eps*s(1)), 'rank deficient: trailing values vanish')
        call require(unit_error(u) <= 20*max(m, n)*eps, 'rank deficient: U unitary')
        call check(a, s, u(:, 1:min(m, n)), vt(1:min(m, n), :), 'rank deficient')
    end subroutine rank_deficient

    subroutine argument_errors()
        complex(dp) :: a(3, 3), u(3, 3), vt(3, 3), work(20)
        real(dp) :: s(3), rwork(15)
        integer :: info
        a = 0
        call fortnum_zgesvd('X', 'N', 3, 3, a, 3, s, u, 3, vt, 3, work, 20, rwork, info)
        call require(info == -1, 'bad jobu')
        call fortnum_zgesvd('O', 'O', 3, 3, a, 3, s, u, 3, vt, 3, work, 20, rwork, info)
        call require(info == -2, 'jobu and jobvt both O')
        call fortnum_zgesvd('N', 'N', 3, 3, a, 2, s, u, 3, vt, 3, work, 20, rwork, info)
        call require(info == -6, 'bad lda')
        call fortnum_zgesvd('N', 'N', 3, 3, a, 3, s, u, 3, vt, 3, work, 2, rwork, info)
        call require(info == -13, 'short work')
    end subroutine argument_errors

    subroutine check(a, s, u, vt, what)
        complex(dp), intent(in) :: a(:, :), u(:, :), vt(:, :)
        real(dp), intent(in) :: s(:)
        character(*), intent(in) :: what
        complex(dp) :: us(size(u, 1), size(u, 2))
        integer :: j
        do j = 1, size(s)
            us(:, j) = u(:, j)*s(j)
        end do
        call require(all(s(:size(s) - 1) >= s(2:)), what//': descending')
        call require(maxval(abs(matmul(us, vt) - a)) <= 50*maxval(shape(a))*eps*max(s(1), tiny(1.0_dp)), &
            what//': U S V^H = A')
        call require(unit_error(conjg(transpose(vt))) <= 20*maxval(shape(a))*eps, what//': V columns orthonormal')
    end subroutine check

    real(dp) function unit_error(q)
        complex(dp), intent(in) :: q(:, :)
        complex(dp) :: g(size(q, 2), size(q, 2))
        integer :: i
        g = matmul(conjg(transpose(q)), q)
        do i = 1, size(g, 1)
            g(i, i) = g(i, i) - 1
        end do
        unit_error = maxval(abs(g))
    end function unit_error

    function general(m, n, seed) result(a)
        integer, intent(in) :: m, n
        real(dp), intent(in) :: seed
        complex(dp) :: a(m, n)
        integer :: i, j
        do j = 1, n
            do i = 1, m
                a(i, j) = cmplx(sin(seed*i + 1.3_dp*j + 0.07_dp*i*j), cos(0.9_dp*seed*j - 0.4_dp*i), dp)
            end do
        end do
    end function general

    ! product of three complex Householder reflectors I - 2 v v^H / v^H v
    function unitary(n, seed) result(q)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed
        complex(dp) :: q(n, n), v(n)
        integer :: i, k
        q = 0
        do i = 1, n
            q(i, i) = 1
        end do
        do k = 1, 3
            do i = 1, n
                v(i) = cmplx(cos(seed*k*i + 0.3_dp*k), sin(0.5_dp*seed*i - k), dp)
            end do
            q = q - 2*matmul(reshape(v, [n, 1]), matmul(reshape(conjg(v), [1, n]), q))/real(dot_product(v, v))
        end do
    end function unitary

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAIL: '//what
            error stop 1
        end if
    end subroutine require

end program test_fortnum_lapack_zgesvd
