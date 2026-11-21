program test_fortnum_blas_workspace_slow
    ! GEMM and TRSM with caller-provided workspace: results against quad
    ! precision with exactly the queried workspace and with the workspace of
    ! one thread, workspace sizes nondecreasing in every dimension, and no
    ! call of the allocator where the C library's malloc can be counted.
    use, intrinsic :: iso_c_binding, only: c_int
    use fortnum_kinds, only: dp
    use fortnum_blas, only: fortnum_dgemm_ws, fortnum_zgemm_ws, fortnum_dgemm_lwork, &
        fortnum_zgemm_lwork, fortnum_dgemm
    use fortnum_blas3, only: fortnum_dtrsm_ws, fortnum_ztrsm_ws, fortnum_dtrsm_lwork, &
        fortnum_ztrsm_lwork
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    interface
        ! calls of malloc, calloc and realloc so far, or -1 where not counted
        integer(c_int) function allocations() bind(c, name='fortnum_test_allocations')
            import :: c_int
        end function allocations
    end interface
    character, parameter :: trans(2) = ['N', 'T'], ztrans(3) = ['N', 'T', 'C']
    character, parameter :: sides(2) = ['L', 'R'], uplos(2) = ['U', 'L']
    ! thin (n <= 96), few rows (m <= 96), blocked with two depth blocks (KC = 512), tiny
    integer, parameter :: shapes(3, 5) = reshape([150, 40, 50, 30, 160, 60, 110, 100, 520, &
        7, 5, 3, 1, 1, 1], [3, 5])
    real(dp), allocatable :: work(:)
    real(dp), allocatable :: warm(:, :)
    integer :: s, i, j, threads

    ! tile selection and the OpenMP thread pool are set up once, before counting
    allocate (warm(400, 400), source=1.0_dp)
    call fortnum_dgemm('N', 'N', 400, 400, 400, 1.0_dp, warm, 400, warm, 400, 0.0_dp, warm, 400)
    do threads = 1, 2
        do s = 1, size(shapes, 2)
            do i = 1, 2
                do j = 1, 2
                    call check_dgemm(trans(i), trans(j), shapes(1, s), shapes(2, s), shapes(3, s), threads)
                end do
            end do
            call check_zgemm('N', 'N', shapes(1, s), shapes(2, s), shapes(3, s), threads)
            do i = 1, 3
                call check_zgemm(ztrans(i), ztrans(1 + mod(i, 3)), shapes(1, s), shapes(2, s), shapes(3, s), threads)
            end do
            do i = 1, 2
                do j = 1, 2
                    call check_dtrsm(sides(i), uplos(j), trans(1 + mod(i + j, 2)), shapes(1, s), shapes(2, s), &
                        threads)
                    call check_ztrsm(sides(i), uplos(j), ztrans(1 + mod(i + j, 3)), shapes(1, s), shapes(2, s), &
                        threads)
                end do
            end do
        end do
    end do
    call check_monotone()
    if (allocations() < 0) then
        print '(a)', 'workspace GEMM and TRSM match quad precision (allocations not counted here)'
    else
        print '(a)', 'workspace GEMM and TRSM match quad precision and allocate nothing'
    end if

contains

    subroutine fill(x, seed)
        real(dp), intent(out) :: x(:, :)
        real(dp), intent(in) :: seed
        integer :: i, j
        do j = 1, size(x, 2)
            do i = 1, size(x, 1)
                x(i, j) = sin(seed + 0.37_dp*i + 1.91_dp*j)*exp(0.3_dp*cos(1.3_dp*i*j))
            end do
        end do
    end subroutine fill

    subroutine expect(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAIL: '//what
            error stop 1
        end if
    end subroutine expect

    integer function workspace(lwork) result(lw)
        integer, intent(in) :: lwork
        lw = lwork
        if (allocated(work)) deallocate (work)
        allocate (work(lw))
        work = huge(1.0_dp)
    end function workspace

    subroutine check_dgemm(ta, tb, m, n, k, threads)
        character, intent(in) :: ta, tb
        integer, intent(in) :: m, n, k, threads
        real(dp), allocatable :: a(:, :), b(:, :), c(:, :), c0(:, :)
        real(qp) :: exact, bound
        real(qp), allocatable :: ga(:, :), gb(:, :)
        integer :: i, j, l, lw, before
        real(dp), parameter :: alpha = 0.75_dp, beta = -1.25_dp
        allocate (a(merge(k, m, ta == 'T'), merge(m, k, ta == 'T')), &
            b(merge(n, k, tb == 'T'), merge(k, n, tb == 'T')), c(m, n), c0(m, n))
        call fill(a, 0.1_dp)
        call fill(b, 0.7_dp)
        call fill(c0, 1.3_dp)
        c = c0
        ! threads 1: the workspace of one thread; 2: the default query
        if (threads == 1) then
            lw = workspace(fortnum_dgemm_lwork(m, n, k, 1))
        else
            lw = workspace(fortnum_dgemm_lwork(m, n, k))
        end if
        before = allocations()
        call fortnum_dgemm_ws(ta, tb, m, n, k, alpha, a, size(a, 1), b, size(b, 1), beta, c, m, work, lw)
        call expect(allocations() == before, 'fortnum_dgemm_ws called the allocator')
        ga = real(a, qp)
        gb = real(b, qp)
        if (ta == 'T') ga = transpose(ga)
        if (tb == 'T') gb = transpose(gb)
        do j = 1, n
            do i = 1, m
                exact = real(beta, qp)*c0(i, j)
                bound = abs(real(beta, qp)*c0(i, j))
                do l = 1, k
                    exact = exact + real(alpha, qp)*ga(i, l)*gb(l, j)
                    bound = bound + abs(real(alpha, qp)*ga(i, l)*gb(l, j))
                end do
                call expect(abs(c(i, j) - exact) <= 4*(k + 2)*epsilon(1.0_dp)*bound + tiny(1.0_dp), &
                    'fortnum_dgemm_ws '//ta//tb)
            end do
        end do
    end subroutine check_dgemm

    subroutine check_zgemm(ta, tb, m, n, k, threads)
        character, intent(in) :: ta, tb
        integer, intent(in) :: m, n, k, threads
        complex(dp), allocatable :: a(:, :), b(:, :), c(:, :), c0(:, :)
        real(dp), allocatable :: re(:, :), im(:, :)
        complex(qp) :: exact
        complex(qp), allocatable :: ga(:, :), gb(:, :)
        real(qp), allocatable :: aa(:, :), ab(:, :)
        real(qp) :: bound
        integer :: i, j, l, lw, before
        complex(dp), parameter :: alpha = (0.75_dp, -0.5_dp), beta = (-1.25_dp, 0.25_dp)
        allocate (a(merge(k, m, ta /= 'N'), merge(m, k, ta /= 'N')), &
            b(merge(n, k, tb /= 'N'), merge(k, n, tb /= 'N')), c(m, n), c0(m, n))
        allocate (re(size(a, 1), size(a, 2)), im(size(a, 1), size(a, 2)))
        call fill(re, 0.2_dp)
        call fill(im, 0.9_dp)
        a = cmplx(re, im, dp)
        deallocate (re, im)
        allocate (re(size(b, 1), size(b, 2)), im(size(b, 1), size(b, 2)))
        call fill(re, 1.1_dp)
        call fill(im, 1.7_dp)
        b = cmplx(re, im, dp)
        deallocate (re, im)
        allocate (re(m, n), im(m, n))
        call fill(re, 2.3_dp)
        call fill(im, 2.9_dp)
        c0 = cmplx(re, im, dp)
        c = c0
        ! threads 1: the workspace of one thread; 2: the default query
        if (threads == 1) then
            lw = workspace(fortnum_zgemm_lwork(m, n, k, 1))
        else
            lw = workspace(fortnum_zgemm_lwork(m, n, k))
        end if
        before = allocations()
        call fortnum_zgemm_ws(ta, tb, m, n, k, alpha, a, size(a, 1), b, size(b, 1), beta, c, m, work, lw)
        call expect(allocations() == before, 'fortnum_zgemm_ws called the allocator')
        ga = cmplx(a, kind=qp)
        gb = cmplx(b, kind=qp)
        if (ta /= 'N') ga = transpose(ga)
        if (ta == 'C') ga = conjg(ga)
        if (tb /= 'N') gb = transpose(gb)
        if (tb == 'C') gb = conjg(gb)
        aa = abs(ga)
        ab = abs(gb)
        do j = 1, n
            do i = 1, m
                exact = cmplx(beta, kind=qp)*c0(i, j)
                bound = abs(exact)/abs(cmplx(alpha, kind=qp))
                do l = 1, k
                    exact = exact + cmplx(alpha, kind=qp)*ga(i, l)*gb(l, j)
                    bound = bound + aa(i, l)*ab(l, j)
                end do
                call expect(abs(c(i, j) - exact) <= 8*(k + 2)*epsilon(1.0_dp)*abs(cmplx(alpha, kind=qp))*bound &
                    + tiny(1.0_dp), 'fortnum_zgemm_ws '//ta//tb)
            end do
        end do
    end subroutine check_zgemm

    ! Well-conditioned triangle: diagonal about 2, off-diagonal small.
    subroutine triangle(a)
        real(dp), intent(out) :: a(:, :)
        integer :: i, j
        do j = 1, size(a, 2)
            do i = 1, size(a, 1)
                a(i, j) = 0.3_dp*sin(0.7_dp*i + 1.3_dp*j)/sqrt(real(size(a, 1), dp))
            end do
            a(j, j) = 2.0_dp + cos(real(j, dp))
        end do
    end subroutine triangle

    ! op(A) X = alpha B (L) or X op(A) = alpha B (R): backward error of X
    ! against quad precision, relative to |op(A)| |X| + |alpha B|.
    subroutine check_dtrsm(side, uplo, ta, m, n, threads)
        character, intent(in) :: side, uplo, ta
        integer, intent(in) :: m, n, threads
        real(dp), allocatable :: a(:, :), b(:, :), x(:, :)
        real(qp), allocatable :: t(:, :), r(:, :), s(:, :)
        integer :: na, lw, before, i, j
        real(dp), parameter :: alpha = 1.5_dp
        na = merge(m, n, side == 'L')
        allocate (a(na, na), b(m, n), t(na, na))
        call triangle(a)
        call fill(b, 0.4_dp)
        x = b
        ! threads 1: the workspace of one thread; 2: the default query
        if (threads == 1) then
            lw = workspace(fortnum_dtrsm_lwork(side, m, n, 1))
        else
            lw = workspace(fortnum_dtrsm_lwork(side, m, n))
        end if
        before = allocations()
        call fortnum_dtrsm_ws(side, uplo, ta, 'N', m, n, alpha, a, na, x, m, work, lw)
        call expect(allocations() == before, 'fortnum_dtrsm_ws called the allocator')
        t = 0
        do j = 1, na
            do i = 1, na
                if ((uplo == 'U' .and. i <= j) .or. (uplo == 'L' .and. i >= j)) t(i, j) = a(i, j)
            end do
        end do
        if (ta == 'T') t = transpose(t)
        if (side == 'L') then
            r = matmul(t, real(x, qp)) - alpha*real(b, qp)
            s = matmul(abs(t), abs(real(x, qp))) + abs(alpha*real(b, qp))
        else
            r = matmul(real(x, qp), t) - alpha*real(b, qp)
            s = matmul(abs(real(x, qp)), abs(t)) + abs(alpha*real(b, qp))
        end if
        call expect(all(abs(r) <= 4*(na + 2)*epsilon(1.0_dp)*s), 'fortnum_dtrsm_ws '//side//uplo//ta)
    end subroutine check_dtrsm

    subroutine check_ztrsm(side, uplo, ta, m, n, threads)
        character, intent(in) :: side, uplo, ta
        integer, intent(in) :: m, n, threads
        complex(dp), allocatable :: a(:, :), b(:, :), x(:, :)
        real(dp), allocatable :: re(:, :), im(:, :)
        complex(qp), allocatable :: t(:, :), r(:, :)
        real(qp), allocatable :: s(:, :)
        integer :: na, lw, before, i, j
        complex(dp), parameter :: alpha = (1.5_dp, -0.5_dp)
        na = merge(m, n, side == 'L')
        allocate (re(na, na), im(na, na))
        call triangle(re)
        call fill(im, 0.8_dp)
        a = cmplx(re, 0.1_dp*im, dp)
        deallocate (re, im)
        allocate (re(m, n), im(m, n))
        call fill(re, 0.4_dp)
        call fill(im, 1.4_dp)
        b = cmplx(re, im, dp)
        x = b
        allocate (t(na, na))
        ! threads 1: the workspace of one thread; 2: the default query
        if (threads == 1) then
            lw = workspace(fortnum_ztrsm_lwork(side, m, n, 1))
        else
            lw = workspace(fortnum_ztrsm_lwork(side, m, n))
        end if
        before = allocations()
        call fortnum_ztrsm_ws(side, uplo, ta, 'N', m, n, alpha, a, na, x, m, work, lw)
        call expect(allocations() == before, 'fortnum_ztrsm_ws called the allocator')
        t = 0
        do j = 1, na
            do i = 1, na
                if ((uplo == 'U' .and. i <= j) .or. (uplo == 'L' .and. i >= j)) t(i, j) = a(i, j)
            end do
        end do
        if (ta == 'T') t = transpose(t)
        if (ta == 'C') t = conjg(transpose(t))
        if (side == 'L') then
            r = matmul(t, cmplx(x, kind=qp)) - alpha*cmplx(b, kind=qp)
            s = matmul(abs(t), abs(cmplx(x, kind=qp))) + abs(alpha*cmplx(b, kind=qp))
        else
            r = matmul(cmplx(x, kind=qp), t) - alpha*cmplx(b, kind=qp)
            s = matmul(abs(cmplx(x, kind=qp)), abs(t)) + abs(alpha*cmplx(b, kind=qp))
        end if
        call expect(all(abs(r) <= 8*(na + 2)*epsilon(1.0_dp)*s), 'fortnum_ztrsm_ws '//side//uplo//ta)
    end subroutine check_ztrsm

    ! one query for the largest product covers every smaller one
    subroutine check_monotone()
        integer, parameter :: dims(9) = [0, 1, 7, 33, 96, 97, 500, 513, 5000]
        integer :: i, j, l, t
        do t = 1, 4
            do l = 1, size(dims)
                do j = 1, size(dims)
                    do i = 1, size(dims) - 1
                        call expect(fortnum_dgemm_lwork(dims(i + 1), dims(j), dims(l), t) >= &
                            fortnum_dgemm_lwork(dims(i), dims(j), dims(l), t) .and. &
                            fortnum_dgemm_lwork(dims(j), dims(i + 1), dims(l), t) >= &
                            fortnum_dgemm_lwork(dims(j), dims(i), dims(l), t) .and. &
                            fortnum_dgemm_lwork(dims(j), dims(l), dims(i + 1), t) >= &
                            fortnum_dgemm_lwork(dims(j), dims(l), dims(i), t) .and. &
                            fortnum_dgemm_lwork(dims(i), dims(j), dims(l), t + 1) >= &
                            fortnum_dgemm_lwork(dims(i), dims(j), dims(l), t), 'fortnum_dgemm_lwork monotone')
                        call expect(fortnum_zgemm_lwork(dims(i + 1), dims(j), dims(l), t) >= &
                            fortnum_zgemm_lwork(dims(i), dims(j), dims(l), t) .and. &
                            fortnum_zgemm_lwork(dims(j), dims(i + 1), dims(l), t) >= &
                            fortnum_zgemm_lwork(dims(j), dims(i), dims(l), t) .and. &
                            fortnum_zgemm_lwork(dims(j), dims(l), dims(i + 1), t) >= &
                            fortnum_zgemm_lwork(dims(j), dims(l), dims(i), t), 'fortnum_zgemm_lwork monotone')
                    end do
                    call expect(fortnum_dtrsm_lwork('L', dims(l) + 1, dims(j), t) >= &
                        fortnum_dtrsm_lwork('L', dims(l), dims(j), t) .and. &
                        fortnum_dtrsm_lwork('R', dims(j), dims(l) + 1, t) >= &
                        fortnum_dtrsm_lwork('R', dims(j), dims(l), t) .and. &
                        fortnum_ztrsm_lwork('L', dims(j), dims(l) + 1, t) >= &
                        fortnum_ztrsm_lwork('L', dims(j), dims(l), t), 'TRSM workspace monotone')
                end do
            end do
        end do
    end subroutine check_monotone

end program test_fortnum_blas_workspace_slow
