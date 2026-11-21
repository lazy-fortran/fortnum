program test_fortnum_lapack_band_chol
    ! Band Cholesky against the dense Cholesky solve of the same matrix, both
    ! triangles, kd from 0 to beyond n, and the info codes.
    use fortnum_kinds, only: dp
    use fortnum_lapack_band_chol
    use fortnum_lapack_chol, only: fortnum_dposv
    implicit none
    integer :: kd

    do kd = 0, 9, 3
        call against_dense(40, kd, 'L')
        call against_dense(40, kd, 'U')
    end do
    call against_dense(5, 7, 'L')
    call against_dense(5, 7, 'U')
    call not_definite()
    call argument_errors()
    print '(a)', 'band Cholesky matches the dense solve'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    subroutine against_dense(n, kd, uplo)
        integer, intent(in) :: n, kd
        character, intent(in) :: uplo
        real(dp) :: a(n, n), ab(kd + 1, n), b(n, 2), x(n, 2)
        integer :: i, j, info
        a = 0
        do j = 1, n
            do i = max(1, j - kd), min(n, j + kd)
                if (i /= j) a(i, j) = sin(0.7_dp*(i + j)) + 0.1_dp*cos(1.3_dp*i*j)
            end do
        end do
        a = 0.5_dp*(a + transpose(a))
        do j = 1, n
            a(j, j) = 3*kd + 1 + 0.3_dp*j
        end do
        ab = -999
        do j = 1, n
            do i = max(1, j - kd), min(n, j + kd)
                if (uplo == 'L' .and. i >= j) ab(1 + i - j, j) = a(i, j)
                if (uplo == 'U' .and. i <= j) ab(kd + 1 + i - j, j) = a(i, j)
            end do
        end do
        do j = 1, 2
            do i = 1, n
                b(i, j) = cos(0.37_dp*i*j)
            end do
        end do
        x = b
        call fortnum_dpbsv(uplo, n, kd, 2, ab, kd + 1, x, n, info)
        call require(info == 0, 'dpbsv info')
        call fortnum_dposv(uplo, n, 2, a, n, b, n, info)
        call require(info == 0, 'dposv info')
        call require(maxval(abs(x - b)) <= 100*n*epsilon(1.0_dp)*maxval(abs(b)), &
            'band solve equals dense solve')
    end subroutine against_dense

    subroutine not_definite()
        real(dp) :: ab(2, 4), b(4, 1)
        integer :: info
        ab(1, :) = [4.0_dp, 4.0_dp, -1.0_dp, 4.0_dp]
        ab(2, :) = [1.0_dp, 1.0_dp, 1.0_dp, 0.0_dp]
        b = 1
        call fortnum_dpbsv('L', 4, 1, 1, ab, 2, b, 4, info)
        call require(info == 3, 'leading minor 3 not positive')
    end subroutine not_definite

    subroutine argument_errors()
        real(dp) :: ab(2, 3), b(3, 1)
        integer :: info
        ab = 1
        b = 1
        call fortnum_dpbtrf('X', 3, 1, ab, 2, info)
        call require(info == -1, 'uplo')
        call fortnum_dpbtrf('L', 3, 2, ab, 2, info)
        call require(info == -5, 'dpbtrf ldab')
        call fortnum_dpbsv('L', 3, 1, 1, ab, 1, b, 3, info)
        call require(info == -6, 'dpbsv ldab')
        call fortnum_dpbsv('L', 3, 1, 1, ab, 2, b, 2, info)
        call require(info == -8, 'dpbsv ldb')
    end subroutine argument_errors

end program test_fortnum_lapack_band_chol
