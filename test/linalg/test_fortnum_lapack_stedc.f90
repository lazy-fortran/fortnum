program test_fortnum_lapack_stedc
    ! Divide and conquer (DSTEDC, and DSYEV/DSYEVD, which use it) on
    ! tridiagonal matrices that exercise every deflation path: random, the
    ! Wilkinson matrix (pairs of nearly equal eigenvalues), the 1-2-1 matrix
    ! (known eigenvalues), glued blocks with tiny couplings, a constant
    ! diagonal with tiny off-diagonals and split blocks. Oracles: the
    ! residual and orthogonality in quad precision and the eigenvalues of
    ! the QL iteration (DSTEQR).
    use fortnum_kinds, only: dp
    use fortnum_lapack_syev, only: fortnum_dstedc, fortnum_dsteqr, fortnum_dsyev, fortnum_dsyevd
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(dp), parameter :: pi = acos(-1.0_dp)
    integer, parameter :: sizes(6) = [2, 25, 26, 51, 100, 257]
    integer :: kind, s, n

    do kind = 1, 6
        do s = 1, size(sizes)
            n = sizes(s)
            call check_tridiagonal(kind, n)
        end do
    end do
    call check_one_two_one(300)
    call check_compz_v(120)
    call check_dense(130)
    call check_query(100)
    print '(a)', 'divide and conquer matches quad-precision residuals and DSTEQR'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    subroutine fill(kind, n, d, e)
        integer, intent(in) :: kind, n
        real(dp), intent(out) :: d(n), e(n)
        integer :: i
        e = 0
        select case (kind)
        case (1)
            do i = 1, n
                d(i) = sin(1.3_dp*i) + 0.2_dp*cos(7.1_dp*i)
                e(i) = cos(0.7_dp*i**2)
            end do
        case (2)
            do i = 1, n
                d(i) = abs(real(n + 1, dp)/2 - i)
            end do
            e = 1
        case (3)
            ! Blocks of the 1-2-1 matrix glued by 1e-14.
            d = 2
            e = -1
            do i = 20, n - 1, 20
                e(i) = 1.0e-14_dp
            end do
        case (4)
            d = 1
            e = 1.0e-16_dp
        case (5)
            ! Exact splits: unreduced blocks of different sizes.
            do i = 1, n
                d(i) = cos(0.31_dp*i)
                e(i) = sin(0.53_dp*i)
            end do
            do i = 7, n - 1, 37
                e(i) = 0
            end do
        case default
            ! Graded entries.
            do i = 1, n
                d(i) = 10.0_dp**(-8*real(i, dp)/n)
                e(i) = 0.5_dp*10.0_dp**(-8*real(i, dp)/n)
            end do
        end select
    end subroutine fill

    subroutine check_tridiagonal(kind, n)
        integer, intent(in) :: kind, n
        real(dp), allocatable :: d(:), e(:), d0(:), e0(:), z(:, :), w(:), work(:), wq(:)
        integer, allocatable :: iwork(:)
        real(dp) :: unused(1, 1)
        integer :: info
        character(len=48) :: label
        allocate (d(n), e(n), z(n, n), work(1), iwork(1))
        call fill(kind, n, d, e)
        d0 = d
        e0 = e
        write (label, '(a,i0,a,i0)') 'kind ', kind, ' n ', n
        call fortnum_dstedc('I', n, d, e, z, n, work, -1, iwork, -1, info)
        call require(info == 0, 'workspace query '//trim(label))
        deallocate (work, iwork)
        allocate (work(int(work_size(n))), iwork(iwork_size(n)))
        call fortnum_dstedc('I', n, d, e, z, n, work, size(work), iwork, size(iwork), info)
        call require(info == 0, 'dstedc info '//trim(label))
        call require(all(d(2:n) >= d(1:n - 1)), 'ascending '//trim(label))
        call require(residual(d0, e0, d, z) <= 20*n*epsilon(1.0_dp)*max(maxval(abs(d0)), &
            maxval(abs(e0(1:n - 1)))), 'residual '//trim(label))
        call require(orthogonality(z) <= 20*n*epsilon(1.0_dp), 'orthogonality '//trim(label))
        ! Values against QL; eigenvalues only.
        wq = d0
        e = e0
        allocate (w(max(1, 2*n - 2)))
        call fortnum_dsteqr('N', n, wq, e, unused, 1, w, info)
        call require(maxval(abs(wq - d)) <= 20*n*epsilon(1.0_dp)*maxval(abs(wq)), &
            'eigenvalues against DSTEQR '//trim(label))
        w = d0
        e = e0
        call fortnum_dstedc('N', n, w, e, unused, 1, work, size(work), iwork, size(iwork), info)
        call require(info == 0 .and. maxval(abs(w(1:n) - d)) <= 20*n*epsilon(1.0_dp)* &
            maxval(abs(wq)), 'COMPZ N '//trim(label))
    end subroutine check_tridiagonal

    subroutine check_one_two_one(n)
        integer, intent(in) :: n
        real(dp) :: d(n), e(n), z(n, n), exact(n)
        real(dp), allocatable :: work(:)
        integer, allocatable :: iwork(:)
        integer :: info, k
        d = 2
        e = -1
        do k = 1, n
            exact(k) = 2 - 2*cos(k*pi/(n + 1))
        end do
        allocate (work(int(work_size(n))), iwork(iwork_size(n)))
        call fortnum_dstedc('I', n, d, e, z, n, work, size(work), iwork, size(iwork), info)
        call require(info == 0 .and. maxval(abs(d - exact)) <= 10*n*epsilon(1.0_dp), &
            'the 1-2-1 eigenvalues are 2 - 2 cos(k pi/(n + 1))')
    end subroutine check_one_two_one

    subroutine check_compz_v(n)
        ! COMPZ V returns Q times the tridiagonal vectors: Q a rotation of
        ! two coordinates.
        integer, intent(in) :: n
        real(dp) :: d(n), e(n), d0(n), e0(n), q(n, n), z(n, n), zi(n, n), di(n), ei(n), c, s
        real(dp), allocatable :: work(:)
        integer, allocatable :: iwork(:)
        integer :: info, j
        call fill(1, n, d, e)
        d0 = d
        e0 = e
        di = d
        ei = e
        c = cos(0.4_dp)
        s = sin(0.4_dp)
        q = 0
        do j = 1, n
            q(j, j) = 1
        end do
        q(1, 1) = c
        q(2, 2) = c
        q(1, 2) = -s
        q(2, 1) = s
        z = q
        allocate (work(1), iwork(1))
        call fortnum_dstedc('V', n, d, e, z, n, work, -1, iwork, -1, info)
        j = nint(work(1))
        info = iwork(1)
        deallocate (work, iwork)
        allocate (work(max(j, work_size(n))), iwork(max(info, iwork_size(n))))
        call fortnum_dstedc('V', n, d, e, z, n, work, size(work), iwork, size(iwork), info)
        call require(info == 0, 'COMPZ V info')
        call fortnum_dstedc('I', n, di, ei, zi, n, work, size(work), iwork, size(iwork), info)
        call require(maxval(abs(z - matmul(q, zi))) <= 100*n*epsilon(1.0_dp), 'COMPZ V is Q Z')
    end subroutine check_compz_v

    subroutine check_dense(n)
        ! DSYEV and DSYEVD on a dense matrix with clustered eigenvalues.
        integer, intent(in) :: n
        real(dp) :: a(n, n), v(n, n), w(n), w2(n)
        real(dp), allocatable :: work(:)
        integer, allocatable :: iwork(:)
        real(qp) :: r(n, n)
        integer :: i, j, info
        do j = 1, n
            do i = j, n
                a(i, j) = sin(0.37_dp*i*j) + merge(real(mod(i, 3), dp), 0.0_dp, i == j)
                a(j, i) = a(i, j)
            end do
        end do
        v = a
        allocate (work(3*n), iwork(1))
        call fortnum_dsyev('V', 'L', n, v, n, w, work, size(work), info)
        call require(info == 0, 'dsyev info')
        r = matmul(real(a, qp), real(v, qp)) - real(v, qp)*spread(real(w, qp), 1, n)
        call require(maxval(abs(r)) <= 20*n*epsilon(1.0_dp)*maxval(abs(w)), 'dsyev residual')
        call require(orthogonality(v) <= 20*n*epsilon(1.0_dp), 'dsyev orthogonality')
        deallocate (work, iwork)
        allocate (work(1 + 6*n + 2*n**2), iwork(3 + 5*n))
        v = a
        call fortnum_dsyevd('V', 'U', n, v, n, w2, work, size(work), iwork, size(iwork), info)
        call require(info == 0 .and. maxval(abs(w2 - w)) <= 20*n*epsilon(1.0_dp)*maxval(abs(w)), &
            'dsyevd eigenvalues')
        r = matmul(real(a, qp), real(v, qp)) - real(v, qp)*spread(real(w2, qp), 1, n)
        call require(maxval(abs(r)) <= 20*n*epsilon(1.0_dp)*maxval(abs(w)), 'dsyevd residual')
    end subroutine check_dense

    subroutine check_query(n)
        integer, intent(in) :: n
        real(dp) :: d(n), e(n), z(n, n), work(1)
        integer :: iwork(1), info
        call fortnum_dstedc('I', n, d, e, z, n, work, -1, iwork, -1, info)
        call require(info == 0 .and. nint(work(1)) == work_size(n) .and. iwork(1) == iwork_size(n), &
            'the workspace query returns LAPACK''s sizes')
        call fortnum_dstedc('I', n, d, e, z, n, work, 1, iwork, 1, info)
        call require(info == -8, 'a short workspace is refused')
        call fortnum_dstedc('X', n, d, e, z, n, work, 1, iwork, 1, info)
        call require(info == -1, 'an unknown COMPZ is refused')
    end subroutine check_query

    integer function work_size(n)
        integer, intent(in) :: n
        work_size = max(1, 1 + 4*n + n**2)
    end function work_size

    integer function iwork_size(n)
        integer, intent(in) :: n
        iwork_size = max(1, 3 + 5*n)
    end function iwork_size

    real(dp) function residual(d0, e0, w, z)
        real(dp), intent(in) :: d0(:), e0(:), w(:), z(:, :)
        real(qp) :: t
        integer :: n, i, j
        n = size(w)
        residual = 0
        do j = 1, n
            do i = 1, n
                t = real(d0(i), qp)*z(i, j) - real(w(j), qp)*z(i, j)
                if (i > 1) t = t + real(e0(i - 1), qp)*z(i - 1, j)
                if (i < n) t = t + real(e0(i), qp)*z(i + 1, j)
                residual = max(residual, real(abs(t), dp))
            end do
        end do
    end function residual

    real(dp) function orthogonality(z)
        real(dp), intent(in) :: z(:, :)
        real(qp), allocatable :: g(:, :)
        integer :: j
        g = matmul(transpose(real(z, qp)), real(z, qp))
        do j = 1, size(z, 2)
            g(j, j) = g(j, j) - 1
        end do
        orthogonality = real(maxval(abs(g)), dp)
    end function orthogonality
end program test_fortnum_lapack_stedc
