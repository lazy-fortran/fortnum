program test_fortnum_lapack_geev
    ! Nonsymmetric eigensolvers by their defining equations in quad
    ! precision: A v = lambda v and u^H A = lambda u^H with LAPACK's
    ! normalization, Schur forms A Z = Z T with orthogonal or unitary Z,
    ! spectra known by construction, balancing of badly scaled and
    ! reducible matrices, extreme magnitudes and argument errors.
    use fortnum_kinds, only: dp
    use fortnum_lapack_geev
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(dp), parameter :: eps = epsilon(1.0_dp)

    call real_general(1, 0.3_dp, 1.0_dp)
    call real_general(2, 0.7_dp, 1.0_dp)
    call real_general(7, 1.1_dp, 1.0_dp)
    call real_general(60, 0.37_dp, 1.0_dp)
    call real_general(120, 0.91_dp, 1.0_dp)
    call real_general(30, 0.53_dp, 1.0e-300_dp)
    call real_general(30, 0.53_dp, 1.0e300_dp)
    call real_known_spectrum(24)
    call real_badly_scaled(40)
    call real_reducible()
    call real_special()
    call real_schur(50)
    call complex_general(1, 0.3_dp)
    call complex_general(9, 1.7_dp)
    call complex_general(70, 0.29_dp)
    call complex_known_spectrum(30)
    call complex_badly_scaled(30)
    call complex_schur(40)
    call argument_errors()
    print '(a)', 'nonsymmetric eigensolvers match quad precision'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    function general(n, seed) result(a)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed
        real(dp) :: a(n, n)
        integer :: i, j
        do j = 1, n
            do i = 1, n
                a(i, j) = sin(seed*i + 1.91_dp*j + 0.37_dp*i*j)*exp(0.3_dp*cos(1.3_dp*i - 0.7_dp*j))
            end do
        end do
    end function general

    function cgeneral(n, seed) result(a)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed
        complex(dp) :: a(n, n)
        a = cmplx(general(n, seed), general(n, seed + 0.61_dp), dp)
    end function cgeneral

    ! eigenvalues of the real driver as complex numbers
    function lambdas(n, wr, wi) result(w)
        integer, intent(in) :: n
        real(dp), intent(in) :: wr(n), wi(n)
        complex(qp) :: w(n)
        w = cmplx(wr, wi, qp)
    end function lambdas

    ! right and left eigenvectors of the real driver as complex vectors
    function real_vectors(n, wi, v) result(x)
        integer, intent(in) :: n
        real(dp), intent(in) :: wi(n), v(n, n)
        complex(qp) :: x(n, n)
        integer :: j
        j = 1
        do while (j <= n)
            if (wi(j) == 0) then
                x(:, j) = cmplx(v(:, j), 0, qp)
                j = j + 1
            else
                x(:, j) = cmplx(v(:, j), v(:, j + 1), qp)
                x(:, j + 1) = conjg(x(:, j))
                j = j + 2
            end if
        end do
    end function real_vectors

    ! max_j |A x_j - lambda_j x_j| / (|A|_F |x_j|) and the same for the left
    ! vectors u_j^H A = lambda_j u_j^H
    subroutine check_pairs(n, a, w, xr, xl, tol, what)
        integer, intent(in) :: n
        complex(qp), intent(in) :: a(n, n), w(n), xr(n, n), xl(n, n)
        real(dp), intent(in) :: tol
        character(*), intent(in) :: what
        real(qp) :: anorm, res
        integer :: j
        anorm = sqrt(sum(abs(a)**2))
        if (anorm == 0) anorm = 1
        do j = 1, n
            res = sqrt(sum(abs(matmul(a, xr(:, j)) - w(j)*xr(:, j))**2))
            call require(res <= tol*anorm*sqrt(sum(abs(xr(:, j))**2)), what//': right residual')
            res = sqrt(sum(abs(matmul(conjg(transpose(a)), xl(:, j)) - conjg(w(j))*xl(:, j))**2))
            call require(res <= tol*anorm*sqrt(sum(abs(xl(:, j))**2)), what//': left residual')
        end do
    end subroutine check_pairs

    ! LAPACK's normalization: norm 1, largest component real
    subroutine check_normalized(n, x, what)
        integer, intent(in) :: n
        complex(qp), intent(in) :: x(n, n)
        character(*), intent(in) :: what
        integer :: j, imax
        do j = 1, n
            call require(abs(sqrt(sum(abs(x(:, j))**2)) - 1) <= 10*eps, what//': norm 1')
            imax = maxloc(abs(x(:, j)), 1)
            call require(abs(aimag(x(imax, j))) <= 10*eps*abs(x(imax, j)), &
                what//': largest component real')
        end do
    end subroutine check_normalized

    ! every lambda of w1 within tol of some lambda of w2, and vice versa
    subroutine check_same_spectrum(n, w1, w2, tol, what)
        integer, intent(in) :: n
        complex(qp), intent(in) :: w1(n), w2(n)
        real(qp), intent(in) :: tol
        character(*), intent(in) :: what
        integer :: j
        do j = 1, n
            call require(minval(abs(w2 - w1(j))) <= tol, what//': spectrum')
            call require(minval(abs(w1 - w2(j))) <= tol, what//': spectrum')
        end do
    end subroutine check_same_spectrum

    subroutine real_driver(n, a, wr, wi, vl, vr, what)
        integer, intent(in) :: n
        real(dp), intent(in) :: a(n, n)
        real(dp), intent(out) :: wr(n), wi(n), vl(n, n), vr(n, n)
        character(*), intent(in) :: what
        real(dp), allocatable :: work(:)
        real(dp) :: b(n, n), q(1)
        integer :: info, j
        b = a
        call fortnum_dgeev('V', 'V', n, b, n, wr, wi, vl, n, vr, n, q, -1, info)
        call require(info == 0 .and. q(1) >= 4*n, what//': workspace query')
        allocate (work(int(q(1))))
        call fortnum_dgeev('V', 'V', n, b, n, wr, wi, vl, n, vr, n, work, size(work), info)
        call require(info == 0, what//': info')
        j = 1
        do while (j <= n)
            if (wi(j) /= 0) then
                call require(j < n, what//': pair complete')
                call require(wi(j) > 0 .and. wi(j + 1) == -wi(j) .and. wr(j + 1) == wr(j), &
                    what//': conjugate pair, positive part first')
                j = j + 2
            else
                j = j + 1
            end if
        end do
    end subroutine real_driver

    subroutine real_general(n, seed, magnitude)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed, magnitude
        real(dp) :: a(n, n), wr(n), wi(n), vl(n, n), vr(n, n), wr2(n), wi2(n), d(1, 1), e(1, 1)
        real(dp) :: work(3*n + 1)
        complex(qp) :: w(n)
        integer :: info
        a = magnitude*general(n, seed)
        call real_driver(n, a, wr, wi, vl, vr, 'real general')
        w = lambdas(n, wr, wi)
        call check_pairs(n, cmplx(a, 0, qp), w, real_vectors(n, wi, vr), real_vectors(n, wi, vl), &
            20*n*eps, 'real general')
        call check_normalized(n, real_vectors(n, wi, vr), 'real right vectors')
        call check_normalized(n, real_vectors(n, wi, vl), 'real left vectors')
        ! eigenvalues alone
        vl = a
        call fortnum_dgeev('N', 'N', n, vl, n, wr2, wi2, d, 1, e, 1, work, size(work), info)
        call require(info == 0, 'real eigenvalues only: info')
        call check_same_spectrum(n, w, lambdas(n, wr2, wi2), &
            real(1.0e-9_dp*magnitude*n, qp), 'real eigenvalues only')
    end subroutine real_general

    ! Q T Q^T with an orthogonal Q (a product of reflectors) and T upper
    ! triangular with 2 x 2 rotation blocks: eigenvalues a(k) +- i b(k)
    ! and the real diagonal entries
    subroutine real_known_spectrum(n)
        integer, intent(in) :: n
        real(qp) :: t(n, n), q(n, n), v(n)
        real(dp) :: a(n, n), wr(n), wi(n), vl(n, n), vr(n, n)
        complex(qp) :: expected(n)
        integer :: i, j, k
        t = 0
        do j = 1, n
            do i = 1, j - 1
                t(i, j) = 0.3_qp*sin(real(i + 2*j, qp))
            end do
        end do
        k = 0
        j = 1
        do while (j <= n)
            k = k + 1
            if (mod(k, 3) == 0 .and. j < n) then
                t(j, j) = -1 + 0.25_qp*k
                t(j + 1, j + 1) = t(j, j)
                t(j, j + 1) = 0.5_qp + 0.1_qp*k
                t(j + 1, j) = -(0.5_qp + 0.1_qp*k)
                expected(j) = cmplx(t(j, j), t(j, j + 1), qp)
                expected(j + 1) = conjg(expected(j))
                j = j + 2
            else
                t(j, j) = 0.75_qp*k - 4
                expected(j) = t(j, j)
                j = j + 1
            end if
        end do
        q = 0
        do i = 1, n
            q(i, i) = 1
        end do
        do k = 1, 3
            do i = 1, n
                v(i) = cos(real(k*i, qp) + 0.3_qp)
            end do
            q = q - 2*spread(v, 2, n)*spread(matmul(v, q), 1, n)/dot_product(v, v)
        end do
        a = real(matmul(matmul(q, t), transpose(q)), dp)
        call real_driver(n, a, wr, wi, vl, vr, 'real known spectrum')
        call check_same_spectrum(n, lambdas(n, wr, wi), expected, 1.0e-12_qp, 'real known spectrum')
    end subroutine real_known_spectrum

    ! D A D^-1 with D spanning 2^-40 .. 2^40: same spectrum, which only
    ! balancing recovers to full accuracy
    subroutine real_badly_scaled(n)
        integer, intent(in) :: n
        real(dp) :: a(n, n), b(n, n), wr(n), wi(n), vl(n, n), vr(n, n), d(n)
        complex(qp) :: w0(n)
        integer :: i, j
        a = general(n, 0.83_dp)
        call real_driver(n, a, wr, wi, vl, vr, 'real unscaled')
        w0 = lambdas(n, wr, wi)
        do i = 1, n
            d(i) = 2.0_dp**(nint(80.0_dp*(i - 1)/(n - 1)) - 40)
        end do
        do j = 1, n
            do i = 1, n
                b(i, j) = d(i)*a(i, j)/d(j)
            end do
        end do
        call real_driver(n, b, wr, wi, vl, vr, 'real badly scaled')
        call check_same_spectrum(n, w0, lambdas(n, wr, wi), 1.0e-11_qp*n, 'real badly scaled')
        call check_pairs(n, cmplx(b, 0, qp), lambdas(n, wr, wi), real_vectors(n, wi, vr), &
            real_vectors(n, wi, vl), 1.0e3_dp*n*eps, 'real badly scaled')
    end subroutine real_badly_scaled

    ! rows and columns that isolate eigenvalues, so balancing permutes
    subroutine real_reducible()
        integer, parameter :: n = 9
        real(dp) :: a(n, n), wr(n), wi(n), vl(n, n), vr(n, n), scl(n)
        integer :: ilo, ihi, info
        a = general(n, 0.44_dp)
        a(4, [1, 2, 3, 5, 6, 7, 8, 9]) = 0
        a([1, 2, 3, 4, 5, 6, 8, 9], 7) = 0
        a(2, [1, 3, 4, 5, 6, 7, 8, 9]) = 0
        vl = a
        call fortnum_dgebal('B', n, vl, n, ilo, ihi, scl, info)
        call require(info == 0 .and. ilo > 1 .and. ihi < n, 'balancing isolates eigenvalues')
        call real_driver(n, a, wr, wi, vl, vr, 'real reducible')
        call check_pairs(n, cmplx(a, 0, qp), lambdas(n, wr, wi), real_vectors(n, wi, vr), &
            real_vectors(n, wi, vl), 20*n*eps, 'real reducible')
        call require(any(abs(wr - a(4, 4)) == 0 .and. wi == 0), 'isolated eigenvalue exact')
    end subroutine real_reducible

    subroutine real_special()
        real(dp) :: a2(2, 2), wr(2), wi(2), vl(2, 2), vr(2, 2)
        real(dp) :: a3(3, 3), wr3(3), wi3(3), vl3(3, 3), vr3(3, 3), d(1, 1), e(1, 1), work(12)
        integer :: info
        ! rotation: +-i
        a2 = reshape([0.0_dp, -1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
        call real_driver(2, a2, wr, wi, vl, vr, 'rotation')
        call require(abs(wr(1)) <= eps .and. abs(wi(1) - 1) <= 2*eps, 'rotation eigenvalues +-i')
        ! defective Jordan block
        a2 = reshape([2.0_dp, 0.0_dp, 1.0_dp, 2.0_dp], [2, 2])
        call real_driver(2, a2, wr, wi, vl, vr, 'Jordan block')
        call require(all(wr == 2) .and. all(wi == 0), 'Jordan block eigenvalues')
        call check_pairs(2, cmplx(a2, 0, qp), lambdas(2, wr, wi), real_vectors(2, wi, vr), &
            real_vectors(2, wi, vl), 10*eps, 'Jordan block')
        ! zero and identity
        a3 = 0
        call real_driver(3, a3, wr3, wi3, vl3, vr3, 'zero matrix')
        call require(all(wr3 == 0) .and. all(wi3 == 0), 'zero matrix eigenvalues')
        a3 = reshape([1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp], [3, 3])
        call real_driver(3, a3, wr3, wi3, vl3, vr3, 'identity')
        call require(all(wr3 == 1) .and. all(wi3 == 0), 'identity eigenvalues')
        ! n = 0
        call fortnum_dgeev('N', 'N', 0, a3, 1, wr3, wi3, d, 1, e, 1, work, 1, info)
        call require(info == 0, 'n = 0')
    end subroutine real_special

    ! DGEHRD, DORGHR and DHSEQR: A Z = Z T, Z orthogonal, T quasi-triangular
    ! with standardized 2 x 2 blocks
    subroutine real_schur(n)
        integer, intent(in) :: n
        real(dp) :: a(n, n), h(n, n), z(n, n), tau(n), wr(n), wi(n), work(n)
        real(qp) :: e(n, n)
        integer :: info, i, j
        a = general(n, 0.71_dp)
        h = a
        call fortnum_dgehrd(n, 1, n, h, n, tau, work, n, info)
        call require(info == 0, 'dgehrd info')
        z = h
        call fortnum_dorghr(n, 1, n, z, n, tau, work, n, info)
        call require(info == 0, 'dorghr info')
        do j = 1, n - 2
            h(j + 2:, j) = 0
        end do
        e = matmul(real(a, qp), real(z, qp)) - matmul(real(z, qp), real(h, qp))
        call require(maxval(abs(e)) <= 50*n*eps*maxval(abs(a)), 'Hessenberg similarity')
        call fortnum_dhseqr('S', 'V', n, 1, n, h, n, wr, wi, z, n, work, n, info)
        call require(info == 0, 'dhseqr info')
        e = matmul(real(a, qp), real(z, qp)) - matmul(real(z, qp), real(h, qp))
        call require(maxval(abs(e)) <= 100*n*eps*maxval(abs(a)), 'real Schur form')
        e = matmul(transpose(real(z, qp)), real(z, qp))
        do i = 1, n
            e(i, i) = e(i, i) - 1
        end do
        call require(maxval(abs(e)) <= 20*n*eps, 'Schur vectors orthonormal')
        do j = 1, n - 1
            call require(all(h(j + 2:, j) == 0), 'quasi-triangular')
            if (h(j + 1, j) /= 0) then
                call require(h(j, j) == h(j + 1, j + 1) .and. h(j + 1, j)*h(j, j + 1) < 0, &
                    'standardized 2 x 2 block')
                call require(abs(wi(j)) > 0, 'block holds a complex pair')
                if (j < n - 1) call require(h(j + 2, j + 1) == 0, 'blocks do not overlap')
            end if
        end do
    end subroutine real_schur

    subroutine complex_driver(n, a, w, vl, vr, what)
        integer, intent(in) :: n
        complex(dp), intent(in) :: a(n, n)
        complex(dp), intent(out) :: w(n), vl(n, n), vr(n, n)
        character(*), intent(in) :: what
        complex(dp), allocatable :: work(:)
        complex(dp) :: b(n, n), q(1)
        real(dp) :: rwork(2*n)
        integer :: info
        b = a
        call fortnum_zgeev('V', 'V', n, b, n, w, vl, n, vr, n, q, -1, rwork, info)
        call require(info == 0 .and. real(q(1)) >= 2*n, what//': workspace query')
        allocate (work(int(real(q(1)))))
        call fortnum_zgeev('V', 'V', n, b, n, w, vl, n, vr, n, work, size(work), rwork, info)
        call require(info == 0, what//': info')
    end subroutine complex_driver

    subroutine complex_general(n, seed)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed
        complex(dp) :: a(n, n), w(n), w2(n), vl(n, n), vr(n, n), d(1, 1), e(1, 1), work(2*n)
        real(dp) :: rwork(2*n)
        integer :: info
        a = cgeneral(n, seed)
        call complex_driver(n, a, w, vl, vr, 'complex general')
        call check_pairs(n, cmplx(a, kind=qp), cmplx(w, kind=qp), cmplx(vr, kind=qp), cmplx(vl, kind=qp), &
            20*n*eps, 'complex general')
        call check_normalized(n, cmplx(vr, kind=qp), 'complex right vectors')
        call check_normalized(n, cmplx(vl, kind=qp), 'complex left vectors')
        vl = a
        call fortnum_zgeev('N', 'N', n, vl, n, w2, d, 1, e, 1, work, size(work), rwork, info)
        call require(info == 0, 'complex eigenvalues only: info')
        call check_same_spectrum(n, cmplx(w, kind=qp), cmplx(w2, kind=qp), real(1.0e-9_dp*n, qp), &
            'complex eigenvalues only')
    end subroutine complex_general

    ! U T U^H with a unitary U and triangular T: eigenvalues diag(T)
    subroutine complex_known_spectrum(n)
        integer, intent(in) :: n
        complex(qp) :: t(n, n), u(n, n), v(n)
        complex(dp) :: a(n, n), w(n), vl(n, n), vr(n, n)
        integer :: i, j, k
        t = 0
        do j = 1, n
            t(j, j) = cmplx(cos(real(j, qp)), 0.5_qp*j - 3, qp)
            do i = 1, j - 1
                t(i, j) = cmplx(0.2_qp*sin(real(i + j, qp)), 0.1_qp*cos(real(i*j, qp)), qp)
            end do
        end do
        u = 0
        do i = 1, n
            u(i, i) = 1
        end do
        do k = 1, 3
            do i = 1, n
                v(i) = cmplx(cos(real(k*i, qp)), sin(real(k + i, qp)), qp)
            end do
            u = u - 2*spread(v, 2, n)*spread(matmul(conjg(v), u), 1, n)/real(dot_product(v, v), qp)
        end do
        a = cmplx(matmul(matmul(u, t), conjg(transpose(u))), kind=dp)
        call complex_driver(n, a, w, vl, vr, 'complex known spectrum')
        call check_same_spectrum(n, cmplx(w, kind=qp), [(t(j, j), j=1, n)], 1.0e-12_qp, &
            'complex known spectrum')
    end subroutine complex_known_spectrum

    subroutine complex_badly_scaled(n)
        integer, intent(in) :: n
        complex(dp) :: a(n, n), b(n, n), w(n), w0(n), vl(n, n), vr(n, n)
        real(dp) :: d(n)
        integer :: i, j
        a = cgeneral(n, 0.47_dp)
        call complex_driver(n, a, w0, vl, vr, 'complex unscaled')
        do i = 1, n
            d(i) = 2.0_dp**(nint(80.0_dp*(i - 1)/(n - 1)) - 40)
        end do
        do j = 1, n
            do i = 1, n
                b(i, j) = d(i)*a(i, j)/d(j)
            end do
        end do
        call complex_driver(n, b, w, vl, vr, 'complex badly scaled')
        call check_same_spectrum(n, cmplx(w0, kind=qp), cmplx(w, kind=qp), 1.0e-11_qp*n, &
            'complex badly scaled')
    end subroutine complex_badly_scaled

    subroutine complex_schur(n)
        integer, intent(in) :: n
        complex(dp) :: a(n, n), h(n, n), z(n, n), tau(n), w(n), work(n)
        complex(qp) :: e(n, n)
        integer :: info, i, j
        a = cgeneral(n, 0.19_dp)
        h = a
        call fortnum_zgehrd(n, 1, n, h, n, tau, work, n, info)
        z = h
        call fortnum_zunghr(n, 1, n, z, n, tau, work, n, info)
        do j = 1, n - 2
            h(j + 2:, j) = 0
        end do
        call fortnum_zhseqr('S', 'V', n, 1, n, h, n, w, z, n, work, n, info)
        call require(info == 0, 'zhseqr info')
        e = matmul(cmplx(a, kind=qp), cmplx(z, kind=qp)) - matmul(cmplx(z, kind=qp), cmplx(h, kind=qp))
        call require(maxval(abs(e)) <= 100*n*eps*maxval(abs(a)), 'complex Schur form')
        e = matmul(conjg(transpose(cmplx(z, kind=qp))), cmplx(z, kind=qp))
        do i = 1, n
            e(i, i) = e(i, i) - 1
        end do
        call require(maxval(abs(e)) <= 20*n*eps, 'Schur vectors unitary')
        do j = 1, n - 1
            call require(all(h(j + 1:, j) == 0), 'triangular')
            call require(w(j) == h(j, j), 'eigenvalues on the diagonal')
        end do
    end subroutine complex_schur

    subroutine argument_errors()
        real(dp) :: a(3, 3), wr(3), wi(3), v(3, 3), u(3, 3), work(12)
        complex(dp) :: c(3, 3), w(3), cv(3, 3), cu(3, 3), cwork(6)
        real(dp) :: rwork(6)
        integer :: info
        a = 1
        call fortnum_dgeev('X', 'N', 3, a, 3, wr, wi, v, 3, u, 3, work, 12, info)
        call require(info == -1, 'dgeev jobvl')
        call fortnum_dgeev('N', 'N', 3, a, 2, wr, wi, v, 3, u, 3, work, 12, info)
        call require(info == -5, 'dgeev lda')
        call fortnum_dgeev('V', 'N', 3, a, 3, wr, wi, v, 2, u, 3, work, 12, info)
        call require(info == -9, 'dgeev ldvl')
        call fortnum_dgeev('N', 'V', 3, a, 3, wr, wi, v, 3, u, 3, work, 11, info)
        call require(info == -13, 'dgeev lwork')
        c = 1
        call fortnum_zgeev('N', 'V', 3, c, 3, w, cv, 3, cu, 2, cwork, 6, rwork, info)
        call require(info == -10, 'zgeev ldvr')
        call fortnum_zgeev('N', 'N', 3, c, 3, w, cv, 3, cu, 3, cwork, 5, rwork, info)
        call require(info == -12, 'zgeev lwork')
    end subroutine argument_errors

end program test_fortnum_lapack_geev
