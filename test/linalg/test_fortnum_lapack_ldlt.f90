program test_fortnum_lapack_ldlt
    ! Bunch-Kaufman LDL^T: solves by their quad-precision backward error for
    ! both triangles, forced 2 x 2 pivots, the pivot encoding, a singular
    ! matrix, and DLASYF's contract (kb = nb - 1 or nb, the rest factored
    ! separately gives the same kind of factorization).
    use fortnum_kinds, only: dp
    use fortnum_lapack_ldlt
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(qp), parameter :: eps = real(epsilon(1.0_dp), qp)

    call solve_case(7, 'L', .false.)
    call solve_case(90, 'L', .false.)
    call solve_case(90, 'U', .false.)
    call solve_case(60, 'L', .true.)
    call solve_case(60, 'U', .true.)
    call panel_case(80, 16, 'L')
    call panel_case(80, 16, 'U')
    call panel_case(10, 32, 'L')
    call singular_case()
    print '(a)', 'LDL^T family matches quad precision'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    ! symmetric indefinite; zero_diagonal forces 2 x 2 pivots
    function sym(n, zero_diagonal) result(a)
        integer, intent(in) :: n
        logical, intent(in) :: zero_diagonal
        real(dp) :: a(n, n)
        integer :: i, j
        do j = 1, n
            do i = j, n
                a(i, j) = sin(0.37_dp*i + 1.91_dp*j)*exp(0.3_dp*cos(1.3_dp*i*j))
                a(j, i) = a(i, j)
            end do
            if (zero_diagonal) a(j, j) = 0
        end do
    end function sym

    subroutine check_ipiv(n, ipiv, upper)
        integer, intent(in) :: n, ipiv(:)
        logical, intent(in) :: upper
        integer :: k
        k = merge(n, 1, upper)
        do while (k >= 1 .and. k <= n)
            if (ipiv(k) > 0) then
                call require(ipiv(k) >= 1 .and. ipiv(k) <= n, 'ipiv in range')
                k = k + merge(-1, 1, upper)
            else
                if (upper) then
                    call require(k > 1 .and. ipiv(k - 1) == ipiv(k), 'upper 2 x 2 pair')
                    k = k - 2
                else
                    call require(k < n .and. ipiv(k + 1) == ipiv(k), 'lower 2 x 2 pair')
                    k = k + 2
                end if
            end if
        end do
    end subroutine check_ipiv

    real(qp) function backward(a, x, b)
        real(dp), intent(in) :: a(:, :), x(:, :), b(:, :)
        backward = maxval(abs(matmul(real(a, qp), real(x, qp)) - real(b, qp)))/ &
            (maxval(sum(abs(real(a, qp)), 2))*maxval(abs(real(x, qp))) + maxval(abs(real(b, qp))))
    end function backward

    subroutine solve_case(n, uplo, zero_diagonal)
        integer, intent(in) :: n
        character, intent(in) :: uplo
        logical, intent(in) :: zero_diagonal
        real(dp) :: a(n, n), f(n + 1, n), b(n, 3), x(n, 3), work(1)
        integer :: ipiv(n), info, i, j, two
        a = sym(n, zero_diagonal)
        do i = 1, 3
            do j = 1, n
                b(j, i) = cos(0.3_dp*i*j)
            end do
        end do
        f = 99
        f(1:n, :) = a
        call fortnum_dsytrf(uplo, n, f, n + 1, ipiv, work, -1, info)
        call require(info == 0 .and. work(1) >= n, 'dsytrf workspace query')
        call fortnum_dsytrf(uplo, n, f, n + 1, ipiv, work, 1, info)
        call require(info == 0, 'dsytrf info')
        call check_ipiv(n, ipiv, uplo == 'U')
        two = count(ipiv < 0)
        if (zero_diagonal) call require(two > 0, 'a zero diagonal forces 2 x 2 pivots')
        x = b
        call fortnum_dsytrs(uplo, n, 3, f, n + 1, ipiv, x, n, info)
        call require(info == 0 .and. backward(a, x, b) <= 16*n*eps, 'dsytrs backward error')
    end subroutine solve_case

    subroutine panel_case(n, nb, uplo)
        integer, intent(in) :: n, nb
        character, intent(in) :: uplo
        real(dp) :: a(n, n), f(n, n), b(n, 2), x(n, 2), w(n, nb)
        integer :: ipiv(n), info, kb, i
        a = sym(n, .true.)
        do i = 1, n
            a(i, i) = 0.05_dp*cos(real(i, dp))
        end do
        call random_like(b)
        f = a
        call fortnum_dlasyf(uplo, n, nb, kb, f, n, ipiv, w, n, info)
        call require(info == 0, 'dlasyf info')
        if (nb >= n) then
            call require(kb == n, 'dlasyf factors everything when n <= nb')
        else
            call require(kb == nb - 1 .or. kb == nb, 'dlasyf kb is nb - 1 or nb')
        end if
        if (kb < n) then
            if (uplo == 'L') then
                call fortnum_dsytf2('L', n - kb, f(kb + 1, kb + 1), n, ipiv(kb + 1), info)
                where (ipiv(kb + 1:) > 0) ipiv(kb + 1:) = ipiv(kb + 1:) + kb
                where (ipiv(kb + 1:) < 0) ipiv(kb + 1:) = ipiv(kb + 1:) - kb
            else
                call fortnum_dsytf2('U', n - kb, f, n, ipiv, info)
            end if
            call require(info == 0, 'dsytf2 on the rest')
        end if
        x = b
        call fortnum_dsytrs(uplo, n, 2, f, n, ipiv, x, n, info)
        call require(info == 0 .and. backward(a, x, b) <= 16*n*eps, 'panel plus rest solves A x = b')
    end subroutine panel_case

    subroutine random_like(x)
        real(dp), intent(out) :: x(:, :)
        integer :: i, j
        do j = 1, size(x, 2)
            do i = 1, size(x, 1)
                x(i, j) = sin(1.7_dp*i - 0.4_dp*j)
            end do
        end do
    end subroutine random_like

    subroutine singular_case()
        real(dp) :: a(3, 3), work(1)
        integer :: ipiv(3), info
        a = 0
        a(1, 1) = 1
        call fortnum_dsytrf('L', 3, a, 3, ipiv, work, 1, info)
        call require(info == 2, 'dsytrf reports the first zero pivot block')
    end subroutine singular_case

end program test_fortnum_lapack_ldlt
