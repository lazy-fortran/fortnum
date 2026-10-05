program test_fortnum_special_legendre_normalized
    ! Behavioral tests for legendre_p_normalized_table.
    ! (1) Closed forms Pbar_0^0, Pbar_1^0, Pbar_1^1, Pbar_2^m, Pbar_3^3 and the
    !     Condon-Shortley sign flip.
    ! (2) Agreement with the independently normalized Ferrers legendre_p.
    ! (3) Orthonormality by exact Gauss-Legendre quadrature (lmax = 100).
    ! (4) Unsold sum rule Pbar_l0^2 + 2 sum_{m>0} Pbar_lm^2 = (2l+1)/2, i.e.
    !     sum_m |Y_lm|^2 = (2l+1)/(4 pi), for all l <= 2000 at x where the
    !     unscaled seed would underflow.
    ! (5) Finite, bounded tables to lmax = 2000 near the poles; NaN outside.
    ! (6) Cross-check of Y_lm against spherical_harmonic.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan
    use fortnum_special, only: legendre_p, legendre_p_normalized_table, &
        spherical_harmonic
    use fortnum_quadrature, only: gauss_legendre
    implicit none

    real(dp), parameter :: eps = epsilon(1.0_dp)
    real(dp), parameter :: pi = acos(-1.0_dp)
    logical :: ok

    ok = .true.
    call check_closed_forms()
    call check_ferrers()
    call check_orthonormality(100)
    call check_sum_rule(2000, [0.0_dp, 0.3_dp, cos(asin(exp(-1.0_dp))), &
                               0.93_dp, 0.999999_dp, 1.0_dp - 1.0e-12_dp, &
                               -0.7_dp, 1.0_dp])
    call check_range()
    call check_harmonic()

    if (.not. ok) then
        write (error_unit, "(a)") "test_fortnum_special_legendre_normalized FAILED"
        error stop 1
    end if
    write (*, "(a)") "test_fortnum_special_legendre_normalized PASSED"

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
        real(dp) :: p(0:3, 0:3), q(0:3, 0:3), x, s, ref(0:3, 0:3), err
        integer :: k, m
        do k = 0, 8
            x = -1.0_dp + 0.25_dp*k
            s = sqrt(1 - x*x)
            ref = 0
            ref(0, 0) = sqrt(0.5_dp)
            ref(1, 0) = sqrt(1.5_dp)*x
            ref(1, 1) = sqrt(0.75_dp)*s
            ref(2, 0) = sqrt(2.5_dp)*(3*x*x - 1)/2
            ref(2, 1) = sqrt(15.0_dp/4)*x*s
            ref(2, 2) = sqrt(15.0_dp/16)*s*s
            ref(3, 0) = sqrt(3.5_dp)*(5*x**3 - 3*x)/2
            ref(3, 3) = sqrt(35.0_dp/32)*s**3
            call legendre_p_normalized_table(3, 3, x, p)
            call legendre_p_normalized_table(3, 3, x, q, condon_shortley=.true.)
            err = max(maxval(abs(p(0:2, 0:2) - ref(0:2, 0:2))), &
                      abs(p(3, 0) - ref(3, 0)), abs(p(3, 3) - ref(3, 3)))
            call expect(err <= 8*eps, "closed form", err)
            do m = 0, 3
                err = maxval(abs(q(:, m) - (-1)**m*p(:, m)))
                call expect(err == 0, "Condon-Shortley sign", err)
            end do
            call expect(p(0, 1) == 0 .and. p(1, 2) == 0 .and. p(2, 3) == 0, &
                        "l < m zero", 0.0_dp)
        end do
    end subroutine check_closed_forms

    ! Pbar_l^m = sqrt((2l+1)/2 (l-m)!/(l+m)!) P_l^m; legendre_p includes the
    ! Condon-Shortley phase and is evaluated through a separate recurrence.
    subroutine check_ferrers()
        integer, parameter :: lmax = 20
        real(dp) :: p(0:lmax, 0:lmax), x, f, ref, err
        integer :: k, l, m, i
        err = 0
        do k = 1, 9
            x = -0.95_dp + 0.21_dp*k
            call legendre_p_normalized_table(lmax, lmax, x, p, &
                                             condon_shortley=.true.)
            do l = 0, lmax
                do m = 0, l
                    f = real(2*l + 1, dp)/2
                    do i = l - m + 1, l + m
                        f = f/real(i, dp)
                    end do
                    ref = sqrt(f)*legendre_p(l, m, x)
                    err = max(err, abs(p(l, m) - ref)/max(1.0_dp, abs(ref)))
                end do
            end do
        end do
        call expect(err <= 1.0e-12_dp, "vs legendre_p", err)
    end subroutine check_ferrers

    subroutine check_orthonormality(lmax)
        integer, intent(in) :: lmax
        integer :: nq, q, m, l, k
        real(dp), allocatable :: x(:), w(:), p(:, :, :)
        real(dp) :: g, err
        nq = lmax + 1
        allocate (x(nq), w(nq), p(0:lmax, 0:lmax, nq))
        call gauss_legendre(nq, x, w)
        do q = 1, nq
            call legendre_p_normalized_table(lmax, lmax, x(q), p(:, :, q))
        end do
        err = 0
        do m = 0, lmax
            do l = m, lmax
                do k = m, l
                    g = sum(w*p(l, m, :)*p(k, m, :))
                    if (k == l) g = g - 1
                    err = max(err, abs(g))
                end do
            end do
        end do
        call expect(err <= 1.0e-12_dp, "orthonormality", err)
    end subroutine check_orthonormality

    subroutine check_sum_rule(lmax, xs)
        integer, intent(in) :: lmax
        real(dp), intent(in) :: xs(:)
        real(dp), allocatable :: p(:, :)
        real(dp) :: t, err
        integer :: i, l
        allocate (p(0:lmax, 0:lmax))
        do i = 1, size(xs)
            call legendre_p_normalized_table(lmax, lmax, xs(i), p)
            err = 0
            do l = 0, lmax
                t = p(l, 0)**2 + 2*sum(p(l, 1:l)**2)
                err = max(err, abs(t/(real(2*l + 1, dp)/2) - 1))
            end do
            ! Observed O(l eps) growth of the degree recurrence near the poles.
            call expect(err <= 25*eps*real(lmax, dp), "sum rule", err)
        end do
    end subroutine check_sum_rule

    subroutine check_range()
        integer, parameter :: lmax = 2000
        real(dp), allocatable :: p(:, :)
        real(dp) :: xs(5)
        integer :: i
        allocate (p(0:lmax, 0:lmax))
        xs = [1.0_dp - 1.0e-15_dp, cos(1.0e-9_dp), 0.5_dp, -1.0_dp, &
              -1.0_dp + 1.0e-300_dp]
        do i = 1, size(xs)
            call legendre_p_normalized_table(lmax, lmax, xs(i), p)
            ! |Pbar_l^m| <= sqrt((2l+1)/2) follows from the sum rule.
            call expect(all(ieee_is_finite(p)) .and. &
                        maxval(abs(p)) <= sqrt(real(2*lmax + 1, dp)/2), &
                        "finite bounded", maxval(abs(p)))
        end do
        call legendre_p_normalized_table(lmax, 0, -1.0_dp, p(:, 0:0))
        call expect(abs(p(lmax, 0) - sqrt(real(2*lmax + 1, dp)/2)) <= &
                    1.0e-12_dp*sqrt(real(lmax, dp)), "pole value", p(lmax, 0))
        call legendre_p_normalized_table(4, 2, 1.0_dp + 1.0e-12_dp, p(0:4, 0:2))
        call expect(all(ieee_is_nan(p(0:4, 0:2))), "NaN outside", 0.0_dp)
    end subroutine check_range

    subroutine check_harmonic()
        integer, parameter :: lmax = 12
        real(dp) :: p(0:lmax, 0:lmax), theta, phi, err
        complex(dp) :: y
        integer :: l, m
        theta = 1.1_dp; phi = 0.7_dp
        call legendre_p_normalized_table(lmax, lmax, cos(theta), p, .true.)
        err = 0
        do l = 0, lmax
            do m = 0, l
                y = p(l, m)*exp(cmplx(0.0_dp, m*phi, dp))/sqrt(2*pi)
                err = max(err, abs(y - spherical_harmonic(l, m, theta, phi)))
            end do
        end do
        call expect(err <= 1.0e-13_dp, "vs spherical_harmonic", err)
    end subroutine check_harmonic

end program test_fortnum_special_legendre_normalized
