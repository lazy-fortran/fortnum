program test_fortnum_blas3_threaded
    ! TRSM split over the columns or rows of B and the triangle-masked SYRK,
    ! on sizes that take the threaded paths, by double-precision residuals:
    ! op(A) X = B (or X op(A) = B), SYRK against the full GEMM, and the other
    ! triangle of C untouched.
    use fortnum_kinds, only: dp
    use fortnum_blas, only: fortnum_dgemm
    use fortnum_blas3, only: fortnum_dtrsm, fortnum_dsyrk
    implicit none
    call trsm_case('L', 'L', 'T', 160, 256)
    call trsm_case('R', 'U', 'N', 256, 160)
    call syrk_case('L', 'N', 170, 150)
    call syrk_case('U', 'T', 170, 150)
    print '(a)', 'threaded TRSM and SYRK match their residuals'

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
                x(i, j) = sin(seed*i + 0.71_dp*j + 0.13_dp*i*j)
            end do
        end do
    end subroutine fill

    subroutine trsm_case(side, uplo, trans, m, n)
        character, intent(in) :: side, uplo, trans
        integer, intent(in) :: m, n
        real(dp), allocatable :: a(:, :), b(:, :), x(:, :), r(:, :)
        integer :: k, i
        k = merge(m, n, side == 'L')
        allocate (a(k, k), b(m, n), x(m, n), r(m, n))
        call fill(a, 0.31_dp)
        do i = 1, k
            a(i, i) = a(i, i) + 4
        end do
        call fill(b, 0.57_dp)
        x = b
        call fortnum_dtrsm(side, uplo, trans, 'N', m, n, 1.0_dp, a, k, x, m)
        ! with the zero triangle cleared, op(A) X or X op(A) by GEMM
        do i = 1, k
            if (uplo == 'L') a(1:i - 1, i) = 0
            if (uplo == 'U') a(i + 1:k, i) = 0
        end do
        if (side == 'L') then
            call fortnum_dgemm(trans, 'N', m, n, m, 1.0_dp, a, k, x, m, 0.0_dp, r, m)
        else
            call fortnum_dgemm('N', trans, m, n, n, 1.0_dp, x, m, a, k, 0.0_dp, r, m)
        end if
        call require(maxval(abs(r - b)) <= 1e-12_dp*maxval(abs(b)), 'trsm residual '//side//uplo//trans)
    end subroutine trsm_case

    subroutine syrk_case(uplo, trans, n, k)
        character, intent(in) :: uplo, trans
        integer, intent(in) :: n, k
        real(dp), allocatable :: a(:, :), c(:, :), c0(:, :), full(:, :)
        integer :: i, j
        if (trans == 'N') then
            allocate (a(n, k))
        else
            allocate (a(k, n))
        end if
        allocate (c(n, n), c0(n, n), full(n, n))
        call fill(a, 0.43_dp)
        call fill(c0, 0.29_dp)
        c = c0
        call fortnum_dsyrk(uplo, trans, n, k, -0.7_dp, a, size(a, 1), 1.3_dp, c, n)
        full = c0
        if (trans == 'N') then
            call fortnum_dgemm('N', 'T', n, n, k, -0.7_dp, a, n, a, n, 1.3_dp, full, n)
        else
            call fortnum_dgemm('T', 'N', n, n, k, -0.7_dp, a, k, a, k, 1.3_dp, full, n)
        end if
        do j = 1, n
            do i = 1, n
                if ((uplo == 'L' .and. i >= j) .or. (uplo == 'U' .and. i <= j)) then
                    call require(abs(c(i, j) - full(i, j)) <= 1e-12_dp*k, 'syrk triangle '//uplo//trans)
                else
                    call require(c(i, j) == c0(i, j), 'syrk wrote the other triangle '//uplo//trans)
                end if
            end do
        end do
    end subroutine syrk_case

end program test_fortnum_blas3_threaded
