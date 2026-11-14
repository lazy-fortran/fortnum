program test_fortnum_blas
    ! fortnum_dgemm and fortnum_zgemm against a direct sum in quad precision,
    ! for every transpose, sizes across the tile edges and the cache blocks,
    ! padded leading dimensions and the special values of alpha and beta.
    use fortnum_kinds, only: dp
    use fortnum_blas, only: fortnum_dgemm, fortnum_zgemm, fortnum_blas_tile
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    integer, parameter :: sizes(3, 9) = reshape([ &
        1, 1, 1, 7, 5, 3, 8, 6, 4, 25, 13, 9, 121, 9, 257, &
        33, 130, 300, 300, 7, 513, 0, 4, 3, 5, 0, 2], [3, 9])
    character, parameter :: real_ops(2) = ['N', 'T'], complex_ops(3) = ['N', 'T', 'C']
    integer :: s, ia, ib, mr, nr, isa, last
    character(len=8) :: mode

    call fortnum_blas_tile(mr, nr, isa)
    print '(a, i0, a, i0, a, i0)', 'tile ', mr, 'x', nr, ', instruction set ', isa
    ! "full" adds the products far beyond one cache block (fortnum_blas_slow)
    last = 5
    if (command_argument_count() > 0) then
        call get_command_argument(1, mode)
        if (mode == "full") last = size(sizes, 2)
    end if
    do s = 1, size(sizes, 2)
        if (s > last .and. s < 8) cycle
        do ia = 1, 2
            do ib = 1, 2
                call check_real(real_ops(ia), real_ops(ib), sizes(1, s), sizes(2, s), sizes(3, s), &
                    1.5_dp, -0.75_dp)
            end do
        end do
        do ia = 1, 3
            do ib = 1, 3
                call check_complex(complex_ops(ia), complex_ops(ib), sizes(1, s), sizes(2, s), &
                    sizes(3, s), (0.5_dp, -1.25_dp), (0.25_dp, 0.5_dp))
            end do
        end do
    end do
    call check_real('n', 't', 40, 30, 20, 2.0_dp, 0.0_dp)
    call check_real('N', 'N', 40, 30, 20, 0.0_dp, 3.0_dp)
    call check_real('N', 'N', 40, 30, 0, 2.0_dp, -1.0_dp)
    call check_complex('c', 'n', 20, 30, 40, (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp))
    call check_complex('N', 'N', 20, 30, 40, (0.0_dp, 0.0_dp), (0.0_dp, 2.0_dp))
    print '(a)', 'fortnum_dgemm and fortnum_zgemm match the quad-precision sums'

contains

    subroutine check_real(ta, tb, m, n, k, alpha, beta)
        character, intent(in) :: ta, tb
        integer, intent(in) :: m, n, k
        real(dp), intent(in) :: alpha, beta
        integer :: lda, ldb, ldc, ra, ca, rb, cb, i, j, p
        real(dp), allocatable :: a(:, :), b(:, :), c(:, :), c0(:, :)
        real(qp) :: exact, bound, x, y

        ra = merge(m, k, ta == 'N' .or. ta == 'n')
        ca = merge(k, m, ta == 'N' .or. ta == 'n')
        rb = merge(k, n, tb == 'N' .or. tb == 'n')
        cb = merge(n, k, tb == 'N' .or. tb == 'n')
        lda = max(1, ra) + 3
        ldb = max(1, rb) + 1
        ldc = max(1, m) + 2
        allocate (a(lda, max(1, ca)), b(ldb, max(1, cb)), c(ldc, max(1, n)))
        call fill(a)
        call fill(b)
        call fill(c)
        if (beta == 0) c(1:m, 1:n) = ieee_value(1.0_dp, ieee_quiet_nan)
        c0 = c
        call fortnum_dgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        do j = 1, n
            do i = 1, m
                exact = 0
                bound = 0
                do p = 1, k
                    x = real(ropa(ta, a, i, p), qp)
                    y = real(ropa(tb, b, p, j), qp)
                    exact = exact + x*y
                    bound = bound + abs(x*y)
                end do
                exact = alpha*exact
                bound = abs(alpha)*bound
                if (beta /= 0) then
                    exact = exact + beta*real(c0(i, j), qp)
                    bound = bound + abs(beta*real(c0(i, j), qp))
                end if
                if (ieee_is_nan(c(i, j)) .or. &
                    abs(real(c(i, j), qp) - exact) > 4*(k + 2)*real(epsilon(1.0_dp), qp)*bound) then
                    print '(a, 2a2, 3i5, 2i5)', 'dgemm mismatch ', ta, tb, m, n, k, i, j
                    error stop 1
                end if
            end do
        end do
        if (any(c(m + 1:, :) /= c0(m + 1:, :))) error stop 'dgemm wrote outside C'


    end subroutine check_real

    subroutine check_complex(ta, tb, m, n, k, alpha, beta)
        character, intent(in) :: ta, tb
        integer, intent(in) :: m, n, k
        complex(dp), intent(in) :: alpha, beta
        integer :: lda, ldb, ldc, ra, ca, rb, cb, i, j, p
        complex(dp), allocatable :: a(:, :), b(:, :), c(:, :), c0(:, :)
        real(dp), allocatable :: re(:, :), im(:, :)
        complex(qp) :: exact, x, y
        real(qp) :: bound

        ra = merge(m, k, ta == 'N' .or. ta == 'n')
        ca = merge(k, m, ta == 'N' .or. ta == 'n')
        rb = merge(k, n, tb == 'N' .or. tb == 'n')
        cb = merge(n, k, tb == 'N' .or. tb == 'n')
        lda = max(1, ra) + 1
        ldb = max(1, rb) + 2
        ldc = max(1, m) + 1
        allocate (a(lda, max(1, ca)), b(ldb, max(1, cb)), c(ldc, max(1, n)))
        allocate (re(lda, max(1, ca)), im(lda, max(1, ca)))
        call fill(re)
        call fill(im)
        a = cmplx(re, im, dp)
        deallocate (re, im)
        allocate (re(ldb, max(1, cb)), im(ldb, max(1, cb)))
        call fill(re)
        call fill(im)
        b = cmplx(re, -im, dp)
        deallocate (re, im)
        allocate (re(ldc, max(1, n)), im(ldc, max(1, n)))
        call fill(re)
        call fill(im)
        c = cmplx(im, re, dp)
        if (beta == (0.0_dp, 0.0_dp)) c(1:m, 1:n) = ieee_value(1.0_dp, ieee_quiet_nan)
        c0 = c
        call fortnum_zgemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        do j = 1, n
            do i = 1, m
                exact = 0
                bound = 0
                do p = 1, k
                    x = cmplx(zop(ta, a, i, p), kind=qp)
                    y = cmplx(zop(tb, b, p, j), kind=qp)
                    exact = exact + x*y
                    bound = bound + abs(x)*abs(y)
                end do
                exact = alpha*exact
                bound = abs(alpha)*bound
                if (beta /= (0.0_dp, 0.0_dp)) then
                    exact = exact + beta*cmplx(c0(i, j), kind=qp)
                    bound = bound + abs(beta)*abs(c0(i, j))
                end if
                if (ieee_is_nan(real(c(i, j))) .or. ieee_is_nan(aimag(c(i, j))) .or. &
                    abs(cmplx(c(i, j), kind=qp) - exact) > 8*(k + 2)*real(epsilon(1.0_dp), qp)*bound) then
                    print '(a, 2a2, 3i5, 2i5)', 'zgemm mismatch ', ta, tb, m, n, k, i, j
                    error stop 1
                end if
            end do
        end do
        if (any(c(m + 1:, :) /= c0(m + 1:, :))) error stop 'zgemm wrote outside C'


    end subroutine check_complex

    ! Element (i, j) of op(x) for a real operand.
    real(dp) function ropa(t, x, i, j)
        character, intent(in) :: t
        real(dp), intent(in) :: x(:, :)
        integer, intent(in) :: i, j
        if (t == 'N' .or. t == 'n') then
            ropa = x(i, j)
        else
            ropa = x(j, i)
        end if
    end function ropa

    ! Element (i, j) of op(x) for a complex operand.
    complex(dp) function zop(t, x, i, j)
        character, intent(in) :: t
        complex(dp), intent(in) :: x(:, :)
        integer, intent(in) :: i, j
        select case (t)
        case ('N', 'n')
            zop = x(i, j)
        case ('T', 't')
            zop = x(j, i)
        case default
            zop = conjg(x(j, i))
        end select
    end function zop

    ! Deterministic values of mixed sign and magnitude.
    subroutine fill(x)
        real(dp), intent(out) :: x(:, :)
        integer :: i, j
        do j = 1, size(x, 2)
            do i = 1, size(x, 1)
                x(i, j) = sin(0.37_dp*i + 1.91_dp*j)*exp(0.3_dp*cos(1.3_dp*i*j))
            end do
        end do
    end subroutine fill

end program test_fortnum_blas
