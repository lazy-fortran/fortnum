module fortnum_lapack_band_chol
    ! Cholesky factorization and solve of symmetric positive definite band
    ! matrices in LAPACK band storage: DPBTRF, DPBTRS, DPBSV with reference
    ! LAPACK arguments and info codes (prefix fortnum_). Column-oriented,
    ! O(n kd^2) work.
    use fortnum_kinds, only: dp
    implicit none
    private
    public :: fortnum_dpbtrf, fortnum_dpbtrs, fortnum_dpbsv

contains

    subroutine fortnum_dpbtrf(uplo, n, kd, ab, ldab, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, kd, ldab
        real(dp), intent(inout) :: ab(ldab, *)
        integer, intent(out) :: info
        real(dp) :: ajj
        integer :: j, kn, c, r
        logical :: upper
        call check(uplo, n, kd, 0, ldab, max(1, n), 5, upper, info)
        if (info /= 0 .or. n == 0) return
        do j = 1, n
            kn = min(kd, n - j)
            if (upper) then
                ! U(j, j + c) = ab(kd + 1 - c, j + c)
                ajj = ab(kd + 1, j)
                if (.not. ajj > 0) then
                    info = j
                    return
                end if
                ajj = sqrt(ajj)
                ab(kd + 1, j) = ajj
                do c = 1, kn
                    ab(kd + 1 - c, j + c) = ab(kd + 1 - c, j + c)/ajj
                end do
                do r = 1, kn
                    do c = 1, r
                        ab(kd + 1 + c - r, j + r) = ab(kd + 1 + c - r, j + r) - &
                            ab(kd + 1 - c, j + c)*ab(kd + 1 - r, j + r)
                    end do
                end do
            else
                ! L(j + r, j) = ab(1 + r, j)
                ajj = ab(1, j)
                if (.not. ajj > 0) then
                    info = j
                    return
                end if
                ajj = sqrt(ajj)
                ab(1, j) = ajj
                ab(2:kn + 1, j) = ab(2:kn + 1, j)/ajj
                do c = 1, kn
                    ab(1:kn - c + 1, j + c) = ab(1:kn - c + 1, j + c) - &
                        ab(1 + c:kn + 1, j)*ab(1 + c, j)
                end do
            end if
        end do
    end subroutine fortnum_dpbtrf

    subroutine fortnum_dpbtrs(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, kd, nrhs, ldab, ldb
        real(dp), intent(in) :: ab(ldab, *)
        real(dp), intent(inout) :: b(ldb, *)
        integer, intent(out) :: info
        integer :: j, k, kn, c
        logical :: upper
        call check(uplo, n, kd, nrhs, ldab, ldb, 6, upper, info)
        if (info /= 0 .or. n == 0) return
        do k = 1, nrhs
            if (upper) then
                ! U^T y = b, then U x = y
                do j = 1, n
                    kn = min(kd, j - 1)
                    do c = 1, kn
                        b(j, k) = b(j, k) - ab(kd + 1 - c, j)*b(j - c, k)
                    end do
                    b(j, k) = b(j, k)/ab(kd + 1, j)
                end do
                do j = n, 1, -1
                    b(j, k) = b(j, k)/ab(kd + 1, j)
                    kn = min(kd, j - 1)
                    do c = 1, kn
                        b(j - c, k) = b(j - c, k) - ab(kd + 1 - c, j)*b(j, k)
                    end do
                end do
            else
                ! L y = b, then L^T x = y
                do j = 1, n
                    kn = min(kd, n - j)
                    b(j, k) = b(j, k)/ab(1, j)
                    b(j + 1:j + kn, k) = b(j + 1:j + kn, k) - ab(2:kn + 1, j)*b(j, k)
                end do
                do j = n, 1, -1
                    kn = min(kd, n - j)
                    b(j, k) = (b(j, k) - dot_product(ab(2:kn + 1, j), b(j + 1:j + kn, k)))/ab(1, j)
                end do
            end if
        end do
    end subroutine fortnum_dpbtrs

    subroutine fortnum_dpbsv(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, kd, nrhs, ldab, ldb
        real(dp), intent(inout) :: ab(ldab, *), b(ldb, *)
        integer, intent(out) :: info
        logical :: upper
        call check(uplo, n, kd, nrhs, ldab, ldb, 6, upper, info)
        if (info /= 0) return
        call fortnum_dpbtrf(uplo, n, kd, ab, ldab, info)
        if (info == 0) call fortnum_dpbtrs(uplo, n, kd, nrhs, ab, ldab, b, ldb, info)
    end subroutine fortnum_dpbsv

    ! Reference argument codes; ldab is argument ldab_position, ldb follows it
    ! two places later (DPBTRS, DPBSV).
    subroutine check(uplo, n, kd, nrhs, ldab, ldb, ldab_position, upper, info)
        character, intent(in) :: uplo
        integer, intent(in) :: n, kd, nrhs, ldab, ldb, ldab_position
        logical, intent(out) :: upper
        integer, intent(out) :: info
        upper = uplo == 'U' .or. uplo == 'u'
        info = 0
        if (.not. (upper .or. uplo == 'L' .or. uplo == 'l')) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (kd < 0) then
            info = -3
        else if (nrhs < 0) then
            info = -4
        else if (ldab < kd + 1) then
            info = -ldab_position
        else if (ldb < max(1, n)) then
            info = -(ldab_position + 2)
        end if
    end subroutine check

end module fortnum_lapack_band_chol
