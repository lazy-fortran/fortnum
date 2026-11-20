module fortnum_lapack_gtsv
    ! Tridiagonal solves with reference LAPACK arguments and info codes
    ! (prefix fortnum_): DGTSV and ZGTSV, Gaussian elimination with partial
    ! pivoting by rows. On return d, du and dl hold the factor (dl the second
    ! superdiagonal of U for interchanged rows) and b the solution.
    use fortnum_kinds, only: dp
    implicit none
    private

    public :: fortnum_dgtsv, fortnum_zgtsv

contains

    subroutine fortnum_dgtsv(n, nrhs, dl, d, du, b, ldb, info)
        integer, intent(in) :: n, nrhs, ldb
        real(dp), intent(inout) :: dl(*), d(*), du(*), b(ldb, *)
        integer, intent(out) :: info
        real(dp) :: mult, temp
        integer :: k, j
        call check(n, nrhs, ldb, info)
        if (info /= 0 .or. n == 0) return
        do k = 1, n - 1
            if (dl(k) == 0) then
                if (d(k) == 0) then
                    info = k
                    return
                end if
            else if (abs(d(k)) >= abs(dl(k))) then
                mult = dl(k)/d(k)
                d(k + 1) = d(k + 1) - mult*du(k)
                b(k + 1, 1:nrhs) = b(k + 1, 1:nrhs) - mult*b(k, 1:nrhs)
                if (k < n - 1) dl(k) = 0
            else
                mult = d(k)/dl(k)
                d(k) = dl(k)
                temp = d(k + 1)
                d(k + 1) = du(k) - mult*temp
                if (k < n - 1) then
                    dl(k) = du(k + 1)
                    du(k + 1) = -mult*dl(k)
                end if
                du(k) = temp
                do j = 1, nrhs
                    temp = b(k, j)
                    b(k, j) = b(k + 1, j)
                    b(k + 1, j) = temp - mult*b(k + 1, j)
                end do
            end if
        end do
        if (d(n) == 0) then
            info = n
            return
        end if
        do j = 1, nrhs
            b(n, j) = b(n, j)/d(n)
            if (n > 1) b(n - 1, j) = (b(n - 1, j) - du(n - 1)*b(n, j))/d(n - 1)
            do k = n - 2, 1, -1
                b(k, j) = (b(k, j) - du(k)*b(k + 1, j) - dl(k)*b(k + 2, j))/d(k)
            end do
        end do
    end subroutine fortnum_dgtsv

    subroutine fortnum_zgtsv(n, nrhs, dl, d, du, b, ldb, info)
        integer, intent(in) :: n, nrhs, ldb
        complex(dp), intent(inout) :: dl(*), d(*), du(*), b(ldb, *)
        integer, intent(out) :: info
        complex(dp) :: mult, temp
        integer :: k, j
        call check(n, nrhs, ldb, info)
        if (info /= 0 .or. n == 0) return
        do k = 1, n - 1
            if (dl(k) == (0.0_dp, 0.0_dp)) then
                if (d(k) == (0.0_dp, 0.0_dp)) then
                    info = k
                    return
                end if
            else if (cabs1(d(k)) >= cabs1(dl(k))) then
                mult = dl(k)/d(k)
                d(k + 1) = d(k + 1) - mult*du(k)
                b(k + 1, 1:nrhs) = b(k + 1, 1:nrhs) - mult*b(k, 1:nrhs)
                if (k < n - 1) dl(k) = 0
            else
                mult = d(k)/dl(k)
                d(k) = dl(k)
                temp = d(k + 1)
                d(k + 1) = du(k) - mult*temp
                if (k < n - 1) then
                    dl(k) = du(k + 1)
                    du(k + 1) = -mult*dl(k)
                end if
                du(k) = temp
                do j = 1, nrhs
                    temp = b(k, j)
                    b(k, j) = b(k + 1, j)
                    b(k + 1, j) = temp - mult*b(k + 1, j)
                end do
            end if
        end do
        if (d(n) == (0.0_dp, 0.0_dp)) then
            info = n
            return
        end if
        do j = 1, nrhs
            b(n, j) = b(n, j)/d(n)
            if (n > 1) b(n - 1, j) = (b(n - 1, j) - du(n - 1)*b(n, j))/d(n - 1)
            do k = n - 2, 1, -1
                b(k, j) = (b(k, j) - du(k)*b(k + 1, j) - dl(k)*b(k + 2, j))/d(k)
            end do
        end do
    end subroutine fortnum_zgtsv

    pure subroutine check(n, nrhs, ldb, info)
        integer, intent(in) :: n, nrhs, ldb
        integer, intent(out) :: info
        info = 0
        if (n < 0) then
            info = -1
        else if (nrhs < 0) then
            info = -2
        else if (ldb < max(1, n)) then
            info = -7
        end if
    end subroutine check

    ! LAPACK's pivoting magnitude |Re z| + |Im z|
    pure real(dp) function cabs1(z)
        complex(dp), intent(in) :: z
        cabs1 = abs(real(z)) + abs(aimag(z))
    end function cabs1

end module fortnum_lapack_gtsv
