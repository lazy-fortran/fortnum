module test_fortnum_lapack_gees_select
    ! eigenvalue selections for ZGEES
    use fortnum_kinds, only: dp
    implicit none
contains
    logical function left_half(w)
        complex(dp), intent(in) :: w
        left_half = real(w, dp) < 0
    end function left_half

    logical function none(w)
        complex(dp), intent(in) :: w
        none = .false.
    end function none
end module test_fortnum_lapack_gees_select

program test_fortnum_lapack_gees
    ! Complex Schur factorization by its defining equations in quad
    ! precision: A Z = Z T with Z unitary and T upper triangular, spectra
    ! known by construction, eigenvalue ordering, swaps of diagonal entries
    ! as similarities, eigenvectors of triangular matrices and of A by their
    ! residuals, plane rotations by their definition, and argument errors.
    use fortnum_kinds, only: dp
    use fortnum_lapack_gees
    use test_fortnum_lapack_gees_select, only: left_half, none
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    real(dp), parameter :: eps = epsilon(1.0_dp)

    call rotations()
    call schur_general(1, 0.3_dp, 1.0_dp)
    call schur_general(2, 0.8_dp, 1.0_dp)
    call schur_general(9, 1.7_dp, 1.0_dp)
    call schur_general(70, 0.29_dp, 1.0_dp)
    call schur_general(30, 0.53_dp, 1.0e-300_dp)
    call schur_general(30, 0.53_dp, 1.0e300_dp)
    call schur_known_spectrum(25)
    call schur_reducible()
    call schur_sorted(40)
    call swaps(12)
    call triangular_vectors(15)
    call triangular_overflow(8)
    call schur_vectors(20)
    call argument_errors()
    print '(a)', 'complex Schur factorization matches quad precision'

contains

    subroutine require(ok, what)
        logical, intent(in) :: ok
        character(*), intent(in) :: what
        if (.not. ok) then
            print '(a)', 'FAILED: '//what
            error stop 1
        end if
    end subroutine require

    function cgeneral(n, seed) result(a)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed
        complex(dp) :: a(n, n)
        integer :: i, j
        do j = 1, n
            do i = 1, n
                a(i, j) = cmplx(sin(seed*i + 1.91_dp*j + 0.37_dp*i*j), &
                    cos(0.7_dp*seed*i - 0.53_dp*j + 0.11_dp*i*j), dp)
            end do
        end do
    end function cgeneral

    function fro(a) result(s)
        complex(qp), intent(in) :: a(:, :)
        real(qp) :: s
        s = sqrt(sum(abs(a)**2))
    end function fro

    ! A Z = Z T, Z^H Z = I, T upper triangular with w on its diagonal
    subroutine check_schur(n, a, t, z, w, tol, what)
        integer, intent(in) :: n
        complex(dp), intent(in) :: a(n, n), t(n, n), z(n, n), w(n)
        real(dp), intent(in) :: tol
        character(*), intent(in) :: what
        complex(qp) :: aq(n, n), tq(n, n), zq(n, n), id(n, n)
        integer :: i, j
        aq = a
        tq = t
        zq = z
        id = 0
        do i = 1, n
            id(i, i) = 1
        end do
        call require(fro(matmul(aq, zq) - matmul(zq, tq)) <= tol*max(fro(aq), tiny(1.0_qp)), &
            what//': A Z = Z T')
        call require(fro(matmul(conjg(transpose(zq)), zq) - id) <= tol, what//': Z unitary')
        do j = 1, n
            do i = j + 1, n
                call require(t(i, j) == 0, what//': T upper triangular')
            end do
            call require(w(j) == t(j, j), what//': eigenvalues on the diagonal')
        end do
    end subroutine check_schur

    subroutine schur_general(n, seed, scale)
        integer, intent(in) :: n
        real(dp), intent(in) :: seed, scale
        complex(dp) :: a(n, n), t(n, n), z(n, n), w(n), w2(n), t2(n, n), zdum(1, 1), work(2*n)
        real(dp) :: rwork(n)
        logical :: bwork(1)
        integer :: sdim, info
        a = cgeneral(n, seed)*scale
        t = a
        call fortnum_zgees('V', 'N', none, n, t, n, sdim, w, z, n, work, size(work), rwork, &
            bwork, info)
        call require(info == 0 .and. sdim == 0, 'zgees general: info')
        call check_schur(n, a, t, z, w, 30*n*eps, 'zgees general')
        ! without Schur vectors: the same triangular factor
        t2 = a
        call fortnum_zgees('N', 'N', none, n, t2, n, sdim, w2, zdum, 1, work, size(work), &
            rwork, bwork, info)
        call require(info == 0, 'zgees without vectors: info')
        call require(all(abs(t2 - t) <= 100*n*eps*maxval(abs(t))), 'zgees without vectors: T')
    end subroutine schur_general

    ! A = U D U^H + strictly upper part, U unitary: eigenvalues D known
    subroutine schur_known_spectrum(n)
        integer, intent(in) :: n
        complex(dp) :: a(n, n), t(n, n), z(n, n), w(n), d(n), u(n, n), work(2*n)
        complex(qp) :: r(n, n)
        real(dp) :: rwork(n)
        logical :: bwork(1)
        integer :: sdim, info, i, j
        do i = 1, n
            d(i) = cmplx(0.5_dp*i - 4, 0.3_dp*(-1)**i*i, dp)
        end do
        r = 0
        do j = 1, n
            r(j, j) = d(j)
            do i = 1, j - 1
                r(i, j) = cmplx(0.2_dp*sin(1.0_dp*i*j), 0.1_dp*cos(1.0_dp*(i + j)), qp)
            end do
        end do
        u = householder(n)
        a = cmplx(matmul(matmul(cmplx(u, kind=qp), r), conjg(transpose(cmplx(u, kind=qp)))), kind=dp)
        t = a
        call fortnum_zgees('V', 'N', none, n, t, n, sdim, w, z, n, work, size(work), rwork, &
            bwork, info)
        call require(info == 0, 'zgees known spectrum: info')
        call check_schur(n, a, t, z, w, 30*n*eps, 'zgees known spectrum')
        do i = 1, n
            call require(minval(abs(w - d(i))) <= 1.0e-10_dp, 'zgees known spectrum: eigenvalue')
        end do
    end subroutine schur_known_spectrum

    ! unitary I - 2 v v^H / v^H v
    function householder(n) result(u)
        integer, intent(in) :: n
        complex(dp) :: u(n, n)
        complex(qp) :: v(n)
        integer :: i, j
        do i = 1, n
            v(i) = cmplx(cos(0.7_qp*i), sin(1.3_qp*i), qp)
        end do
        do j = 1, n
            do i = 1, n
                u(i, j) = cmplx(-2*v(i)*conjg(v(j))/sum(abs(v)**2), kind=dp)
            end do
            u(j, j) = u(j, j) + 1
        end do
    end function householder

    ! upper triangular up to a permutation: balancing isolates every eigenvalue
    subroutine schur_reducible()
        integer, parameter :: n = 6
        integer, parameter :: p(n) = [4, 1, 6, 2, 5, 3]
        complex(dp) :: r(n, n), a(n, n), t(n, n), z(n, n), w(n), work(2*n)
        real(dp) :: rwork(n)
        logical :: bwork(1)
        integer :: sdim, info, i, j
        r = 0
        do j = 1, n
            do i = 1, j
                r(i, j) = cmplx(i + 0.5_dp*j, i - j, dp)
            end do
        end do
        do j = 1, n
            do i = 1, n
                a(p(i), p(j)) = r(i, j)
            end do
        end do
        t = a
        call fortnum_zgees('V', 'N', none, n, t, n, sdim, w, z, n, work, size(work), rwork, &
            bwork, info)
        call require(info == 0, 'zgees reducible: info')
        call check_schur(n, a, t, z, w, 30*n*eps, 'zgees reducible')
        do i = 1, n
            call require(minval(abs(w - r(i, i))) <= 10*eps*abs(r(i, i)), &
                'zgees reducible: eigenvalue')
        end do
    end subroutine schur_reducible

    ! sort 'S': the eigenvalues in the left half plane lead the diagonal
    subroutine schur_sorted(n)
        integer, intent(in) :: n
        complex(dp) :: a(n, n), t(n, n), z(n, n), w(n), w0(n), t0(n, n), z0(n, n), work(2*n)
        real(dp) :: rwork(n)
        logical :: bwork(n)
        integer :: sdim, info, i
        a = cgeneral(n, 0.77_dp) - 0.2_dp
        t0 = a
        call fortnum_zgees('V', 'N', none, n, t0, n, sdim, w0, z0, n, work, size(work), rwork, &
            bwork, info)
        t = a
        call fortnum_zgees('V', 'S', left_half, n, t, n, sdim, w, z, n, work, size(work), &
            rwork, bwork, info)
        call require(info == 0, 'zgees sorted: info')
        call check_schur(n, a, t, z, w, 60*n*eps, 'zgees sorted')
        call require(sdim == count(real(w0, dp) < 0) .and. sdim > 0 .and. sdim < n, &
            'zgees sorted: sdim')
        call require(all(real(w(1:sdim), dp) < 0) .and. all(real(w(sdim + 1:n), dp) >= 0), &
            'zgees sorted: order')
        call require(all(bwork(1:n) .eqv. real(w0, dp) < 0), 'zgees sorted: bwork')
        do i = 1, n
            call require(minval(abs(w - w0(i))) <= 1.0e-12_dp*maxval(abs(w0)), &
                'zgees sorted: spectrum kept')
        end do
    end subroutine schur_sorted

    ! ztrexc: T1 = Q^H T0 Q with Q unitary, T1 upper triangular and its
    ! diagonal the diagonal of T0 with entry ifst moved to ilst
    subroutine swaps(n)
        integer, intent(in) :: n
        integer, parameter :: moves(2, 4) = reshape([1, 12, 12, 1, 3, 9, 10, 2], [2, 4])
        complex(dp) :: t0(n, n), t(n, n), q(n, n), q2(n, n), t2(n, n), d(n)
        complex(qp) :: id(n, n)
        integer :: i, j, k, info, ifst, ilst
        t0 = 0
        do j = 1, n
            do i = 1, j
                t0(i, j) = cmplx(sin(1.3_dp*i + 0.4_dp*j), cos(0.9_dp*i*j), dp)
            end do
        end do
        id = 0
        do i = 1, n
            id(i, i) = 1
        end do
        do k = 1, size(moves, 2)
            ifst = moves(1, k)
            ilst = moves(2, k)
            t = t0
            q = cmplx(id, kind=dp)
            call fortnum_ztrexc('V', n, t, n, q, n, ifst, ilst, info)
            call require(info == 0, 'ztrexc: info')
            call require(fro(matmul(matmul(conjg(transpose(cmplx(q, kind=qp))), &
                cmplx(t0, kind=qp)), cmplx(q, kind=qp)) - t) <= 30*n*eps*fro(cmplx(t0, kind=qp)), &
                'ztrexc: similarity')
            call require(fro(matmul(conjg(transpose(cmplx(q, kind=qp))), cmplx(q, kind=qp)) &
                - id) <= 30*n*eps, 'ztrexc: Q unitary')
            do j = 1, n
                call require(all(t(j + 1:n, j) == 0), 'ztrexc: triangular')
            end do
            d = [(t0(i, i), i = 1, n)]
            if (ifst < ilst) then
                d(ifst:ilst) = [d(ifst + 1:ilst), d(ifst)]
            else
                d(ilst:ifst) = [d(ifst), d(ilst:ifst - 1)]
            end if
            call require(all(abs([(t(i, i), i = 1, n)] - d) <= 30*n*eps), 'ztrexc: diagonal')
            ! compq 'N' gives the same T
            t2 = t0
            call fortnum_ztrexc('N', n, t2, n, q2, 1, ifst, ilst, info)
            call require(info == 0 .and. all(t2 == t), 'ztrexc without Q')
        end do
    end subroutine swaps

    ! ztrevc on an upper triangular T: T x = lambda x and y^H T = lambda y^H,
    ! for all ('A') and selected ('S') vectors, largest |re| + |im| equal 1
    subroutine triangular_vectors(n)
        integer, intent(in) :: n
        integer, parameter :: picked(3) = [2, 7, 11]
        complex(dp) :: t(n, n), vr(n, n), vl(n, n), work(2*n)
        complex(qp) :: tq(n, n)
        real(dp) :: rwork(n)
        logical :: sel(n)
        integer :: i, j, m, info, k
        t = 0
        do j = 1, n
            do i = 1, j
                t(i, j) = cmplx(cos(0.8_dp*i + 1.1_dp*j), sin(0.3_dp*i*j + 0.2_dp), dp)
            end do
            t(j, j) = cmplx(j, 0.5_dp*j, dp)
        end do
        ! a repeated eigenvalue: the divisor is perturbed to ulp |t(4,4)|
        t(4, 4) = t(3, 3)
        tq = t
        call fortnum_ztrevc('B', 'A', sel, n, t, n, vl, n, vr, n, n, m, work, rwork, info)
        call require(info == 0 .and. m == n, 'ztrevc all: info')
        do j = 1, n
            call require(all(vr(j + 1:n, j) == 0) .and. all(vl(1:j - 1, j) == 0), &
                'ztrevc all: triangular vectors')
            call check_vector(tq, cmplx(vr(:, j), kind=qp), tq(j, j), .false., 'ztrevc all right')
            call check_vector(tq, cmplx(vl(:, j), kind=qp), tq(j, j), .true., 'ztrevc all left')
        end do
        sel = .false.
        sel(picked) = .true.
        vr = 0
        vl = 0
        call fortnum_ztrevc('R', 'S', sel, n, t, n, vl, 1, vr, n, 3, m, work, rwork, info)
        call require(info == 0 .and. m == 3, 'ztrevc selected: info')
        do k = 1, 3
            j = picked(k)
            call check_vector(tq, cmplx(vr(:, k), kind=qp), tq(j, j), .false., 'ztrevc selected')
        end do
        call require(all(vr(:, 4:) == 0), 'ztrevc selected: m columns')
    end subroutine triangular_vectors

    ! ztrevc with Q = I and howmny 'B' on a triangular matrix whose
    ! eigenvectors overflow without rescaling: off-diagonal 1e150, gaps 1
    subroutine triangular_overflow(n)
        integer, intent(in) :: n
        complex(dp) :: t(n, n), vr(n, n), vl(n, n), work(2*n)
        complex(qp) :: tq(n, n)
        real(dp) :: rwork(n)
        logical :: sel(1)
        integer :: i, j, m, info
        t = 0
        do j = 1, n
            do i = 1, j - 1
                t(i, j) = cmplx(1.0e150_dp, 0.5e150_dp*sin(1.0_dp*i*j), dp)
            end do
            t(j, j) = j
        end do
        tq = t
        vr = 0
        vl = 0
        do i = 1, n
            vr(i, i) = 1
            vl(i, i) = 1
        end do
        call fortnum_ztrevc('B', 'B', sel, n, t, n, vl, n, vr, n, n, m, work, rwork, info)
        call require(info == 0 .and. m == n, 'ztrevc overflow: info')
        do j = 1, n
            call check_vector(tq, cmplx(vr(:, j), kind=qp), tq(j, j), .false., 'ztrevc overflow right')
            call check_vector(tq, cmplx(vl(:, j), kind=qp), tq(j, j), .true., 'ztrevc overflow left')
        end do
    end subroutine triangular_overflow

    ! |M x - lambda x| small and x normalized (left: M^H y = conj(lambda) y)
    subroutine check_vector(mq, x, lambda, left, what)
        complex(qp), intent(in) :: mq(:, :), x(:), lambda
        logical, intent(in) :: left
        character(*), intent(in) :: what
        real(qp) :: res
        if (left) then
            res = sqrt(sum(abs(matmul(conjg(transpose(mq)), x) - conjg(lambda)*x)**2))
        else
            res = sqrt(sum(abs(matmul(mq, x) - lambda*x)**2))
        end if
        call require(res <= 100*size(x)*eps*fro(mq)*sqrt(sum(abs(x)**2)), what//': residual')
        call require(abs(maxval(abs(real(x, qp)) + abs(aimag(x))) - 1) <= 4*eps, &
            what//': normalized')
    end subroutine check_vector

    ! ztrevc 'B' after zgees: eigenvectors of A itself
    subroutine schur_vectors(n)
        integer, intent(in) :: n
        complex(dp) :: a(n, n), t(n, n), z(n, n), w(n), work(2*n), vl(n, n), vr(n, n)
        real(dp) :: rwork(n)
        logical :: bwork(1), sel(1)
        integer :: sdim, info, m, j
        a = cgeneral(n, 1.23_dp)
        t = a
        call fortnum_zgees('V', 'N', none, n, t, n, sdim, w, z, n, work, size(work), rwork, &
            bwork, info)
        call require(info == 0, 'zgees for vectors: info')
        vr = z
        vl = z
        call fortnum_ztrevc('B', 'B', sel, n, t, n, vl, n, vr, n, n, m, work, rwork, info)
        call require(info == 0 .and. m == n, 'ztrevc back-transformed: info')
        do j = 1, n
            call check_vector(cmplx(a, kind=qp), cmplx(vr(:, j), kind=qp), cmplx(w(j), kind=qp), &
                .false., 'ztrevc back-transformed right')
            call check_vector(cmplx(a, kind=qp), cmplx(vl(:, j), kind=qp), cmplx(w(j), kind=qp), &
                .true., 'ztrevc back-transformed left')
        end do
    end subroutine schur_vectors

    ! [c s; -conj(s) c] [f; g] = [r; 0], c real, c^2 + |s|^2 = 1
    subroutine rotations()
        complex(dp), parameter :: fs(5) = [(3.0_dp, -4.0_dp), (0.0_dp, 0.0_dp), (1.0_dp, 2.0_dp), &
            (1.0e300_dp, 1.0e300_dp), (1.0e-300_dp, 0.0_dp)]
        complex(dp), parameter :: gs(5) = [(1.0_dp, 1.0_dp), (0.0_dp, -2.0_dp), (0.0_dp, 0.0_dp), &
            (-1.0e300_dp, 2.0e300_dp), (0.0_dp, 3.0e-300_dp)]
        complex(dp) :: s, r
        complex(qp) :: f, g
        real(dp) :: c
        real(qp) :: scale
        integer :: k
        do k = 1, size(fs)
            call fortnum_zlartg(fs(k), gs(k), c, s, r)
            f = fs(k)
            g = gs(k)
            scale = max(abs(f), abs(g))
            call require(c >= 0 .and. abs(c**2 + abs(s)**2 - 1) <= 4*eps, 'zlartg: unitary')
            call require(abs(c*f + s*g - r) <= 4*eps*scale, 'zlartg: first row')
            call require(abs(-conjg(cmplx(s, kind=qp))*f + c*g) <= 4*eps*scale, 'zlartg: zero')
        end do
    end subroutine rotations

    subroutine argument_errors()
        complex(dp) :: a(2, 2), z(2, 2), w(2), work(4)
        real(dp) :: rwork(2)
        logical :: bwork(2)
        integer :: sdim, info, m
        a = 0
        call fortnum_zgees('X', 'N', none, 2, a, 2, sdim, w, z, 2, work, 4, rwork, bwork, info)
        call require(info == -1, 'zgees: jobvs')
        call fortnum_zgees('V', 'X', none, 2, a, 2, sdim, w, z, 2, work, 4, rwork, bwork, info)
        call require(info == -2, 'zgees: sort')
        call fortnum_zgees('V', 'N', none, -1, a, 2, sdim, w, z, 2, work, 4, rwork, bwork, info)
        call require(info == -4, 'zgees: n')
        call fortnum_zgees('V', 'N', none, 2, a, 1, sdim, w, z, 2, work, 4, rwork, bwork, info)
        call require(info == -6, 'zgees: lda')
        call fortnum_zgees('V', 'N', none, 2, a, 2, sdim, w, z, 1, work, 4, rwork, bwork, info)
        call require(info == -10, 'zgees: ldvs')
        call fortnum_zgees('V', 'N', none, 2, a, 2, sdim, w, z, 2, work, 3, rwork, bwork, info)
        call require(info == -12, 'zgees: lwork')
        call fortnum_zgees('V', 'N', none, 2, a, 2, sdim, w, z, 2, work, -1, rwork, bwork, info)
        call require(info == 0 .and. int(real(work(1), dp)) >= 4, 'zgees: workspace query')
        call fortnum_ztrexc('X', 2, a, 2, z, 2, 1, 2, info)
        call require(info == -1, 'ztrexc: compq')
        call fortnum_ztrexc('V', 2, a, 2, z, 1, 1, 2, info)
        call require(info == -6, 'ztrexc: ldq')
        call fortnum_ztrexc('V', 2, a, 2, z, 2, 3, 2, info)
        call require(info == -7, 'ztrexc: ifst')
        call fortnum_ztrexc('V', 2, a, 2, z, 2, 1, 0, info)
        call require(info == -8, 'ztrexc: ilst')
        call fortnum_ztrevc('X', 'A', bwork, 2, a, 2, z, 2, z, 2, 2, m, work, rwork, info)
        call require(info == -1, 'ztrevc: side')
        call fortnum_ztrevc('R', 'X', bwork, 2, a, 2, z, 2, z, 2, 2, m, work, rwork, info)
        call require(info == -2, 'ztrevc: howmny')
        call fortnum_ztrevc('L', 'A', bwork, 2, a, 2, z, 1, z, 2, 2, m, work, rwork, info)
        call require(info == -8, 'ztrevc: ldvl')
        call fortnum_ztrevc('R', 'A', bwork, 2, a, 2, z, 2, z, 1, 2, m, work, rwork, info)
        call require(info == -10, 'ztrevc: ldvr')
        call fortnum_ztrevc('R', 'A', bwork, 2, a, 2, z, 2, z, 2, 1, m, work, rwork, info)
        call require(info == -11, 'ztrevc: mm')
    end subroutine argument_errors

end program test_fortnum_lapack_gees
