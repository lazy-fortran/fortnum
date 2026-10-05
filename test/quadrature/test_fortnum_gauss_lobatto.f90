program test_fortnum_gauss_lobatto
    ! Behavioral tests for gauss_lobatto_legendre and the barycentric
    ! Lagrange differentiation matrix.
    ! (1) Closed-form rules for n = 2..6 (A&S Table 25.6).
    ! (2) Exact rational integrals of x^k, k <= 2n - 3; x^(2n-2) is not exact.
    ! (3) Nodes and weights against an independent quad-precision Newton
    !     refinement of the zeros of P_N' (n up to 201).
    ! (4) Exact symmetry and endpoint values.
    ! (5) Differentiation matrix against a quad-precision product formula with
    !     the direct diagonal sum_k 1/(x_i - x_k) on the same nodes, against the
    !     GLL closed form D_ij = P_N(x_i)/(P_N(x_j)(x_i - x_j)),
    !     D_11 = -N(N+1)/4 (within the O(N^2 eps) node conditioning), and
    !     exact differentiation of polynomials on nonsymmetric nodes.
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_quadrature, only: gauss_lobatto_legendre
    use fortnum_polynomial, only: barycentric_weights, &
        lagrange_differentiation_matrix
    implicit none

    real(dp), parameter :: eps = epsilon(1.0_dp)
    integer, parameter :: big_cases(6) = [10, 33, 64, 100, 200, 201]
    logical :: ok
    integer :: j

    ok = .true.
    call check_closed_forms()
    do j = 2, 40
        call check_exactness(j)
    end do
    call check_exactness(120)
    call check_exactness(200)
    do j = 1, size(big_cases)
        call check_quad_reference(big_cases(j))
    end do
    call check_derivative_matrix(5)
    call check_derivative_matrix(17)
    call check_derivative_matrix(64)
    call check_derivative_matrix(200)
    call check_general_nodes()

    if (.not. ok) then
        write (error_unit, "(a)") "test_fortnum_gauss_lobatto FAILED"
        error stop 1
    end if
    write (*, "(a)") "test_fortnum_gauss_lobatto PASSED"

contains

    subroutine expect(cond, label, err)
        logical, intent(in) :: cond
        character(*), intent(in) :: label
        real(dp), intent(in) :: err
        if (.not. cond) then
            write (error_unit, "(a,a,es12.4)") "FAIL ", label, err
            ok = .false.
        end if
    end subroutine expect

    subroutine check_closed_forms()
        real(dp) :: x2(2), w2(2), x3(3), w3(3), x4(4), w4(4), x5(5), w5(5)
        real(dp) :: x6(6), w6(6), r(6), rw(6), t
        call gauss_lobatto_legendre(2, x2, w2)
        call expect(all(x2 == [-1.0_dp, 1.0_dp]) .and. all(w2 == 1.0_dp), &
                    "n=2", 0.0_dp)
        call gauss_lobatto_legendre(3, x3, w3)
        t = maxval(abs(w3 - [1.0_dp, 4.0_dp, 1.0_dp]/3.0_dp))
        call expect(x3(2) == 0.0_dp .and. t <= 2*eps, "n=3", t)
        call gauss_lobatto_legendre(4, x4, w4)
        t = max(abs(x4(3) - sqrt(0.2_dp)), &
                maxval(abs(w4 - [1.0_dp, 5.0_dp, 5.0_dp, 1.0_dp]/6.0_dp)))
        call expect(t <= 2*eps, "n=4", t)
        call gauss_lobatto_legendre(5, x5, w5)
        t = max(abs(x5(4) - sqrt(3.0_dp/7.0_dp)), maxval(abs(w5 - &
                [0.1_dp, 49.0_dp/90.0_dp, 32.0_dp/45.0_dp, 49.0_dp/90.0_dp, 0.1_dp])))
        call expect(x5(3) == 0.0_dp .and. t <= 2*eps, "n=5", t)
        call gauss_lobatto_legendre(6, x6, w6)
        r = [-1.0_dp, -sqrt(1.0_dp/3 + 2*sqrt(7.0_dp)/21), &
             -sqrt(1.0_dp/3 - 2*sqrt(7.0_dp)/21), &
             sqrt(1.0_dp/3 - 2*sqrt(7.0_dp)/21), &
             sqrt(1.0_dp/3 + 2*sqrt(7.0_dp)/21), 1.0_dp]
        rw = [1.0_dp/15, (14 - sqrt(7.0_dp))/30, (14 + sqrt(7.0_dp))/30, &
              (14 + sqrt(7.0_dp))/30, (14 - sqrt(7.0_dp))/30, 1.0_dp/15]
        t = max(maxval(abs(x6 - r)), maxval(abs(w6 - rw)))
        call expect(t <= 4*eps, "n=6", t)
    end subroutine check_closed_forms

    subroutine check_exactness(n)
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), q, exact, err, pn, pm
        integer :: k
        call gauss_lobatto_legendre(n, x, w)
        do k = 0, 2*n - 3
            q = sum(w*x**k)
            exact = 0.0_dp
            if (mod(k, 2) == 0) exact = 2.0_dp/real(k + 1, dp)
            err = abs(q - exact)
            call expect(err <= 8*eps*real(n, dp)*max(1.0_dp, abs(exact)), &
                        "exactness", err)
        end do
        ! Degree 2n - 2 is not exact: the rule gives sum w P_N^2 = 2/N while
        ! int P_N^2 = 2/(2N + 1) (Canuto et al., Spectral Methods, eq. 2.3.13).
        q = 0
        do k = 1, n
            call legendre_dp(n - 1, x(k), pn, pm)
            q = q + w(k)*pn*pn
        end do
        err = abs(q - 2.0_dp/real(n - 1, dp))
        call expect(err <= 8*eps*real(n, dp), "P_N^2 discrete norm", err)
        call expect(x(1) == -1.0_dp .and. x(n) == 1.0_dp, "endpoints", 0.0_dp)
        call expect(all(x(1:n) == -x(n:1:-1)) .and. all(w(1:n) == w(n:1:-1)), &
                    "symmetry", 0.0_dp)
        call expect(all(x(2:n) > x(1:n - 1)) .and. all(w > 0), "order", 0.0_dp)
    end subroutine check_exactness

    ! Independent oracle: quad-precision Newton on P_N' from each node.
    subroutine check_quad_reference(n)
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), ex, ew
        real(qp) :: xq, pn, pnm1, d1, d2, wq
        integer :: i, it, deg
        deg = n - 1
        call gauss_lobatto_legendre(n, x, w)
        ex = 0; ew = 0
        do i = 2, n - 1
            xq = real(x(i), qp)
            do it = 1, 6
                call legendre_qp(deg, xq, pn, pnm1)
                d1 = deg*(xq*pn - pnm1)/(xq*xq - 1)
                d2 = (2*xq*d1 - deg*(deg + 1)*pn)/(1 - xq*xq)
                xq = xq - d1/d2
            end do
            call legendre_qp(deg, xq, pn, pnm1)
            wq = 2/(real(deg, qp)*(deg + 1)*pn*pn)
            ex = max(ex, real(abs(real(x(i), qp) - xq), dp))
            ew = max(ew, real(abs(real(w(i), qp) - wq)/wq, dp))
        end do
        call expect(ex <= 4*eps, "node vs quad", ex)
        ! P_N by recurrence carries O(N eps) relative error; weights ~ P_N^-2.
        call expect(ew <= 8*eps*real(n, dp), "weight vs quad", ew)
        ew = abs(w(1) - 2.0_dp/(real(deg, dp)*(deg + 1)))/w(1)
        call expect(ew <= 2*eps, "endpoint weight", ew)
    end subroutine check_quad_reference

    pure subroutine legendre_qp(n, x, p, pm)
        integer, intent(in) :: n
        real(qp), intent(in) :: x
        real(qp), intent(out) :: p, pm
        real(qp) :: pn
        integer :: k
        pm = 1; p = x
        do k = 1, n - 1
            pn = ((2*k + 1)*x*p - k*pm)/(k + 1)
            pm = p; p = pn
        end do
    end subroutine legendre_qp

    subroutine check_derivative_matrix(n)
        integer, intent(in) :: n
        real(dp) :: x(n), w(n), d(n, n), pn(n), pm(n), ref, err, errq, c
        real(qp) :: lam(n), refq
        integer :: i, j, k, deg
        deg = n - 1
        call gauss_lobatto_legendre(n, x, w)
        call lagrange_differentiation_matrix(n, x, d)
        do i = 1, n
            call legendre_dp(deg, x(i), pn(i), pm(i))
            lam(i) = 1
            do k = 1, n
                if (k /= i) lam(i) = lam(i)*(real(x(i), qp) - real(x(k), qp))
            end do
        end do
        c = real(deg, dp)*real(deg + 1, dp)/4
        err = 0; errq = 0
        do j = 1, n
            do i = 1, n
                if (i /= j) then
                    ref = pn(i)/(pn(j)*(x(i) - x(j)))
                    refq = (lam(i)/lam(j))/(real(x(i), qp) - real(x(j), qp))
                else
                    ref = 0
                    if (i == 1) ref = -c
                    if (i == n) ref = c
                    refq = 0
                    do k = 1, n
                        if (k /= i) refq = refq + 1/(real(x(i), qp) - real(x(k), qp))
                    end do
                end if
                err = max(err, abs(d(i, j) - ref))
                errq = max(errq, real(abs(real(d(i, j), qp) - refq), dp))
            end do
        end do
        call expect(errq <= eps*c*real(n, dp), "D vs quad product", errq)
        call expect(err <= eps*c*real(n, dp)**2/4, "GLL D closed form", err)
    end subroutine check_derivative_matrix

    pure subroutine legendre_dp(n, x, p, pm)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp), intent(out) :: p, pm
        real(dp) :: pn
        integer :: k
        pm = 1; p = x
        do k = 1, n - 1
            pn = ((2*k + 1)*x*p - k*pm)/(k + 1)
            pm = p; p = pn
        end do
    end subroutine legendre_dp

    ! Nonsymmetric mapped nodes on [0.3, 2.1]: D x^k = k x^(k-1), k < n, and
    ! barycentric interpolation reproduces a polynomial off the nodes.
    subroutine check_general_nodes()
        integer, parameter :: n = 9
        real(dp) :: xi(n), w(n), x(n), d(n, n), bw(n), err, xe, num, den, f
        integer :: k
        call gauss_lobatto_legendre(n, xi, w)
        x = 0.3_dp + 0.9_dp*(xi + 1)**1.3_dp/2**0.3_dp
        call lagrange_differentiation_matrix(n, x, d)
        err = 0
        do k = 0, n - 1
            err = max(err, maxval(abs(matmul(d, x**k) - k*x**max(k - 1, 0))) &
                      /max(1.0_dp, k*2.1_dp**(k - 1)))
        end do
        call expect(err <= 1.0e-12_dp, "D on mapped nodes", err)
        call barycentric_weights(n, x, bw)
        xe = 1.234_dp
        num = sum(bw*(x**8 - 3*x)/(xe - x))
        den = sum(bw/(xe - x))
        f = xe**8 - 3*xe
        err = abs(num/den - f)/abs(f)
        call expect(err <= 1.0e-13_dp, "barycentric interpolation", err)
    end subroutine check_general_nodes

end program test_fortnum_gauss_lobatto
