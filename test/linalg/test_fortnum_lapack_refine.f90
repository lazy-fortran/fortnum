program test_fortnum_lapack_refine
    ! Norms, condition estimates, iterative refinement with error bounds and
    ! the expert driver, against exact values computed in quad precision.
    use fortnum_kinds, only: dp
    use fortnum_lapack_lu, only: fortnum_dgetrf, fortnum_dgetri, fortnum_zgetrf, fortnum_zgetri
    use fortnum_lapack_ldlt, only: fortnum_dsytrf
    use fortnum_lapack_refine
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)

    call norms()
    call condition(60)
    call refinement(50)
    call symmetric_refinement(40)
    call expert(30, .false.)
    call expert(30, .true.)
    call nearly_singular()
    print '(a)', 'condition estimates, refinement and DGESVX behave as LAPACK'

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

    subroutine norms()
        real(dp) :: a(3, 2), work(3)
        complex(dp) :: z(2, 2)
        a = reshape([1.0_dp, -2.0_dp, 3.0_dp, -4.0_dp, 5.0_dp, 0.5_dp], [3, 2])
        call require(fortnum_dlange('M', 3, 2, a, 3, work) == 5, 'max norm')
        call require(fortnum_dlange('1', 3, 2, a, 3, work) == 9.5_dp, '1-norm')
        call require(fortnum_dlange('I', 3, 2, a, 3, work) == 7, 'infinity norm')
        call require(abs(fortnum_dlange('F', 3, 2, a, 3, work) - sqrt(55.25_dp)) < 1e-14_dp, 'Frobenius norm')
        z = reshape([(3.0_dp, 4.0_dp), (0.0_dp, 1.0_dp), (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)], [2, 2])
        call require(abs(fortnum_zlange('1', 2, 2, z, 2, work) - 6) < 1e-14_dp, 'complex 1-norm')
    end subroutine norms

    subroutine condition(n)
        integer, intent(in) :: n
        real(dp) :: a(n, n), f(n, n), inv(n, n), work(4*n), rcond, anorm, re(n, n), im(n, n)
        complex(dp) :: z(n, n), zf(n, n), zinv(n, n), zwork(2*n)
        real(qp) :: exact
        integer :: ipiv(n), info, iwork(n), k
        character :: norm
        call fill(a, 0.41_dp)
        do k = 1, 2
            norm = '1I'(k:k)
            f = a
            call fortnum_dgetrf(n, n, f, n, ipiv, info)
            anorm = fortnum_dlange(norm, n, n, a, n, work)
            call fortnum_dgecon(norm, n, f, n, anorm, rcond, work, iwork, info)
            inv = f
            call fortnum_dgetri(n, inv, n, ipiv, work, 4*n, info)
            if (norm == '1') then
                exact = 1/(real(anorm, qp)*maxval(sum(abs(real(inv, qp)), 1)))
            else
                exact = 1/(real(anorm, qp)*maxval(sum(abs(real(inv, qp)), 2)))
            end if
            call require(info == 0 .and. rcond >= exact*(1 - 1e-8_qp) .and. rcond <= 10*exact, 'dgecon '//norm)
        end do
        call fill(re, 0.2_dp)
        call fill(im, 0.9_dp)
        z = cmplx(re, im, dp)
        zf = z
        call fortnum_zgetrf(n, n, zf, n, ipiv, info)
        anorm = fortnum_zlange('1', n, n, z, n, work)
        call fortnum_zgecon('1', n, zf, n, anorm, rcond, zwork, work, info)
        zinv = zf
        call fortnum_zgetri(n, zinv, n, ipiv, zwork, 2*n, info)
        exact = 1/(real(anorm, qp)*maxval(sum(abs(cmplx(zinv, kind=qp)), 1)))
        call require(info == 0 .and. rcond >= exact*(1 - 1e-8_qp) .and. rcond <= 10*exact, 'zgecon')
    end subroutine condition

    ! from a perturbed solution, refinement must reach a backward error of a
    ! few ulp and bound the true forward error
    subroutine refinement(n)
        integer, intent(in) :: n
        real(dp) :: a(n, n), f(n, n), b(n, 2), x(n, 2), ferr(2), berr(2), work(3*n)
        real(qp) :: exact(n, 2), aq(n, n), t
        integer :: ipiv(n), info, iwork(n), j, k
        character :: tr
        call fill(a, 0.33_dp)
        call fill(b, 0.7_dp)
        aq = real(a, qp)
        do k = 1, 2
            tr = 'NT'(k:k)
            if (tr == 'T') aq = transpose(real(a, qp))
            exact = solve_q(aq, real(b, qp))
            f = a
            call fortnum_dgetrf(n, n, f, n, ipiv, info)
            x = real(exact, dp)*(1 + 1e-6_dp)
            call fortnum_dgerfs(tr, n, 2, a, n, f, n, ipiv, b, n, x, n, ferr, berr, work, iwork, info)
            call require(info == 0, 'dgerfs info')
            do j = 1, 2
                call require(berr(j) <= 4*epsilon(1.0_dp), 'dgerfs reaches a backward error of a few ulp '//tr)
                t = maxval(abs(real(x(:, j), qp) - exact(:, j)))/maxval(abs(exact(:, j)))
                call require(t <= ferr(j) .and. ferr(j) <= 1e-8_dp, 'dgerfs forward bound holds '//tr)
            end do
        end do
    end subroutine refinement

    subroutine symmetric_refinement(n)
        integer, intent(in) :: n
        real(dp) :: a(n, n), f(n, n), b(n, 1), x(n, 1), ferr(1), berr(1), work(64*n)
        real(qp) :: exact(n, 1), t
        integer :: ipiv(n), info, iwork(n), i, j
        call fill(a, 0.47_dp)
        do j = 1, n
            do i = j + 1, n
                a(j, i) = a(i, j)
            end do
        end do
        call fill(b, 0.8_dp)
        exact = solve_q(real(a, qp), real(b, qp))
        f = a
        call fortnum_dsytrf('L', n, f, n, ipiv, work, 64*n, info)
        x = real(exact, dp)*(1 - 1e-7_dp)
        call fortnum_dsyrfs('L', n, 1, a, n, f, n, ipiv, b, n, x, n, ferr, berr, work, iwork, info)
        t = maxval(abs(real(x(:, 1), qp) - exact(:, 1)))/maxval(abs(exact(:, 1)))
        call require(info == 0 .and. berr(1) <= 4*epsilon(1.0_dp) .and. t <= ferr(1), 'dsyrfs')
    end subroutine symmetric_refinement

    ! badly scaled rows: fact E equilibrates; fact F reuses the factors
    subroutine expert(n, scaled)
        integer, intent(in) :: n
        logical, intent(in) :: scaled
        real(dp) :: a(n, n), a0(n, n), af(n, n), b(n, 1), b0(n, 1), x(n, 1), r(n), c(n), ferr(1), berr(1)
        real(dp) :: work(4*n), rcond
        real(qp) :: exact(n, 1), t
        integer :: ipiv(n), info, iwork(n), i
        character :: equed
        call fill(a0, 0.29_dp)
        if (scaled) then
            do i = 1, n
                a0(i, :) = a0(i, :)*10.0_dp**(mod(i, 7)*3)
            end do
        end if
        call fill(b0, 0.6_dp)
        if (scaled) b0(:, 1) = b0(:, 1)*10.0_dp**(mod([(i, i = 1, n)], 7)*3)
        exact = solve_q(real(a0, qp), real(b0, qp))
        a = a0
        b = b0
        call fortnum_dgesvx(merge('E', 'N', scaled), 'N', n, 1, a, n, af, n, ipiv, equed, r, c, b, n, x, n, &
            rcond, ferr, berr, work, iwork, info)
        call require(info == 0, 'dgesvx info')
        if (scaled) call require(equed == 'R' .or. equed == 'B', 'dgesvx equilibrates the rows')
        if (.not. scaled) call require(equed == 'N', 'dgesvx leaves a balanced matrix alone')
        t = maxval(abs(real(x(:, 1), qp) - exact(:, 1)))/maxval(abs(exact(:, 1)))
        call require(t <= ferr(1) .and. ferr(1) < 1e-6_dp, 'dgesvx forward bound')
        call require(work(1) > 0 .and. work(1) <= 1 + 1e-12_dp, 'dgesvx pivot growth')
        ! reuse: fact F with the same equilibration and factors
        b = b0
        call fortnum_dgesvx('F', 'N', n, 1, a, n, af, n, ipiv, equed, r, c, b, n, x, n, &
            rcond, ferr, berr, work, iwork, info)
        t = maxval(abs(real(x(:, 1), qp) - exact(:, 1)))/maxval(abs(exact(:, 1)))
        call require(info == 0 .and. t <= ferr(1), 'dgesvx with given factors')
    end subroutine expert

    subroutine nearly_singular()
        real(dp) :: a(3, 3), af(3, 3), b(3, 1), x(3, 1), r(3), c(3), ferr(1), berr(1), work(12), rcond
        integer :: ipiv(3), iwork(3), info
        character :: equed
        a = reshape([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp + 1e-15_dp], [3, 3])
        b = 1
        call fortnum_dgesvx('N', 'N', 3, 1, a, 3, af, 3, ipiv, equed, r, c, b, 3, x, 3, &
            rcond, ferr, berr, work, iwork, info)
        call require(info == 4 .or. info == 3, 'dgesvx flags a matrix singular to working precision')
    end subroutine nearly_singular

    ! Gaussian elimination with partial pivoting in quad precision
    function solve_q(a, b) result(x)
        real(qp), intent(in) :: a(:, :), b(:, :)
        real(qp) :: x(size(b, 1), size(b, 2)), m(size(a, 1), size(a, 1) + size(b, 2)), row(size(a, 1) + size(b, 2))
        integer :: n, k, p, i
        n = size(a, 1)
        m(:, 1:n) = a
        m(:, n + 1:) = b
        do k = 1, n
            p = k - 1 + maxloc(abs(m(k:n, k)), 1)
            row = m(k, :)
            m(k, :) = m(p, :)
            m(p, :) = row
            do i = k + 1, n
                m(i, k:) = m(i, k:) - m(i, k)/m(k, k)*m(k, k:)
            end do
        end do
        do k = n, 1, -1
            m(k, n + 1:) = (m(k, n + 1:) - matmul(m(k, k + 1:n), m(k + 1:n, n + 1:)))/m(k, k)
        end do
        x = m(:, n + 1:)
    end function solve_q

end program test_fortnum_lapack_refine
