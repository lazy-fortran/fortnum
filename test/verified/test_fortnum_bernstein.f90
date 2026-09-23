!> Bernstein evaluation on real and complex boxes, checked against the exact
!> rational Bernstein sum evaluated in real128 (real and complex128
!> arguments) at many sampled points, which is a different, direct
!> implementation of the same mathematical object than the forward-
!> difference/local-Taylor algorithm under test. `bernstein_hull` is checked
!> the same way, plus its derivative bound against the exact derivative
!> formula. `bernstein_lsq_fit` is checked by recovering the coefficients of
!> a known polynomial from noise-free samples.
program test_fortnum_bernstein
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_interval, only: interval_t, cinterval_t, interval, cinterval, &
        real_part, imag_part
    use fortnum_bernstein, only: bernstein_eval_real_box, &
        bernstein_eval_complex_box, bernstein_hull, bernstein_grid_t, &
        bernstein_grid_build, bernstein_grid_eval_real, bernstein_lsq_fit
    implicit none

    integer, parameter :: n = 4
    real(dp) :: coef(0:n)
    integer :: i, nfail, nbad, seed_size
    integer, allocatable :: seed(:)

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 71
    call random_seed(put=seed)
    call random_number(coef)
    coef = 2.0_dp*coef - 1.0_dp

    nbad = 0
    do i = 1, 2000
        call check_real_box(coef, nbad)
    end do
    call require(nbad == 0, "real-box evaluation and hull bound contain "// &
        "sampled exact rational values and derivatives", nfail)

    nbad = 0
    do i = 1, 2000
        call check_complex_box(coef, nbad)
    end do
    call require(nbad == 0, "complex-box evaluation contains sampled "// &
        "exact rational values", nfail)

    nbad = 0
    do i = 1, 300
        call check_grid(coef, nbad)
    end do
    call require(nbad == 0, "grid-cached evaluation agrees with the "// &
        "uncached local re-expansion", nfail)

    call check_lsq_fit(nfail)

    deallocate (seed)
    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " bernstein test(s) failed"
        error stop 1
    end if
    print '(a)', "fortnum_bernstein: all tests passed"

contains

    !> Exact Bernstein sum p(s) = sum_j coef(j) C(n,j) s^j (1-s)^(n-j) at a
    !> real128 point, and its exact derivative, both by direct evaluation of
    !> the defining formula (no forward differences, no basis recursion).
    pure subroutine exact_real(coef, sq, p, dp_)
        real(dp), intent(in) :: coef(0:n)
        real(qp), intent(in) :: sq
        real(qp), intent(out) :: p, dp_
        real(qp) :: bino, dbino
        integer :: j

        p = 0.0_qp
        dp_ = 0.0_qp
        bino = 1.0_qp
        do j = 0, n
            if (j > 0) bino = bino*real(n - j + 1, qp)/real(j, qp)
            p = p + real(coef(j), qp)*bino*sq**j*(1.0_qp - sq)**(n - j)
        end do
        dbino = 1.0_qp
        do j = 0, n - 1
            if (j > 0) dbino = dbino*real(n - 1 - j + 1, qp)/real(j, qp)
            dp_ = dp_ + real(n, qp)*(real(coef(j + 1), qp) - real(coef(j), qp)) &
                *dbino*sq**j*(1.0_qp - sq)**(n - 1 - j)
        end do
    end subroutine exact_real

    pure function exact_complex(coef, zq) result(p)
        real(dp), intent(in) :: coef(0:n)
        complex(qp), intent(in) :: zq
        complex(qp) :: p
        real(qp) :: bino
        integer :: j

        p = (0.0_qp, 0.0_qp)
        bino = 1.0_qp
        do j = 0, n
            if (j > 0) bino = bino*real(n - j + 1, qp)/real(j, qp)
            p = p + real(coef(j), qp)*bino*zq**j*(1.0_qp - zq)**(n - j)
        end do
    end function exact_complex

    subroutine check_real_box(coef, nbad)
        real(dp), intent(in) :: coef(0:n)
        integer, intent(inout) :: nbad
        real(dp) :: r(3), a, b, c
        type(interval_t) :: box, value, vhull, dhull
        real(qp) :: sq, p, dp_
        integer :: j

        call random_number(r)
        a = 0.02_dp + 0.9_dp*r(1)
        b = a + 0.001_dp + 0.05_dp*r(2)
        b = min(b, 1.0_dp)
        c = a + r(3)*(b - a)
        box = interval(a, b)

        call bernstein_eval_real_box(n, coef, c, box, value)
        call bernstein_hull(n, coef, a, b, vhull, dhull)
        do j = 0, 4
            sq = real(a, qp) + real(b - a, qp)*real(j, qp)/4.0_qp
            call exact_real(coef, sq, p, dp_)
            if (.not. encl(value, p)) nbad = nbad + 1
            if (.not. encl(vhull, p)) nbad = nbad + 1
            if (.not. encl(dhull, dp_)) nbad = nbad + 1
        end do
    end subroutine check_real_box

    subroutine check_complex_box(coef, nbad)
        real(dp), intent(in) :: coef(0:n)
        integer, intent(inout) :: nbad
        real(dp) :: r(4)
        real(dp) :: xa, xb, ya, yb, c
        type(cinterval_t) :: box, value
        complex(qp) :: zq, p
        integer :: j

        call random_number(r)
        xa = 0.05_dp + 0.8_dp*r(1)
        xb = xa + 0.001_dp + 0.03_dp*r(2)
        ya = -0.05_dp - 0.05_dp*r(3)
        yb = 0.05_dp + 0.05_dp*r(4)
        c = 0.5_dp*(xa + xb)
        box = cinterval(interval(xa, xb), interval(ya, yb))

        call bernstein_eval_complex_box(n, coef, c, box, value)
        do j = 0, 4
            zq = cmplx(real(xa, qp) + real(xb - xa, qp)*real(j, qp)/4.0_qp, &
                real(ya, qp) + real(yb - ya, qp)*real(4 - j, qp)/4.0_qp, qp)
            p = exact_complex(coef, zq)
            if (.not. encl_c(value, p)) nbad = nbad + 1
        end do
    end subroutine check_complex_box

    subroutine check_grid(coef, nbad)
        real(dp), intent(in) :: coef(0:n)
        integer, intent(inout) :: nbad
        real(dp) :: r(2), a, b
        type(interval_t) :: box, value_grid
        type(bernstein_grid_t) :: grid
        real(qp) :: sq, p, dp_
        integer :: j

        call bernstein_grid_build(n, coef, 64, grid)
        call random_number(r)
        a = 0.05_dp + 0.8_dp*r(1)
        b = a + 0.001_dp + 0.02_dp*r(2)
        box = interval(a, b)
        call bernstein_grid_eval_real(grid, box, value_grid)
        do j = 0, 4
            sq = real(a, qp) + real(b - a, qp)*real(j, qp)/4.0_qp
            call exact_real(coef, sq, p, dp_)
            if (.not. encl(value_grid, p)) nbad = nbad + 1
        end do
    end subroutine check_grid

    subroutine check_lsq_fit(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: true_coef(0:n), s(50), y(50), fit(0:n), rms, t(0:n), dt(0:n)
        integer :: i, info, j
        logical :: ok

        true_coef = [0.3_dp, -0.7_dp, 1.1_dp, 0.2_dp, -0.4_dp]
        do i = 1, 50
            s(i) = real(i - 1, dp)/49.0_dp
            call basis_at(n, s(i), t)
            y(i) = dot_product(t, true_coef)
        end do
        call bernstein_lsq_fit(n, s, y, fit, rms, info)
        ok = info == 0 .and. rms < 1.0e-8_dp
        do j = 0, n
            ok = ok .and. abs(fit(j) - true_coef(j)) < 1.0e-6_dp
        end do
        call require(ok, "least-squares Bernstein fit recovers a known "// &
            "polynomial's coefficients from noise-free samples", nfail)
    end subroutine check_lsq_fit

    pure subroutine basis_at(deg, s, t)
        integer, intent(in) :: deg
        real(dp), intent(in) :: s
        real(dp), intent(out) :: t(0:deg)
        real(dp) :: bino
        integer :: j
        bino = 1.0_dp
        do j = 0, deg
            if (j > 0) bino = bino*real(deg - j + 1, dp)/real(j, dp)
            t(j) = bino*s**j*(1.0_dp - s)**(deg - j)
        end do
    end subroutine basis_at

    logical function encl(v, e)
        type(interval_t), intent(in) :: v
        real(qp), intent(in) :: e
        encl = real(v%lo, qp) <= e .and. e <= real(v%hi, qp)
    end function encl

    logical function encl_c(v, e)
        type(cinterval_t), intent(in) :: v
        complex(qp), intent(in) :: e
        type(interval_t) :: re, im
        re = real_part(v)
        im = imag_part(v)
        encl_c = real(re%lo, qp) <= real(e, qp) .and. &
            real(e, qp) <= real(re%hi, qp) .and. &
            real(im%lo, qp) <= aimag(e) .and. aimag(e) <= real(im%hi, qp)
    end function encl_c

    subroutine require(cond, msg, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(inout) :: nfail
        if (.not. cond) then
            write (error_unit, '(a,a)') "FAIL: ", msg
            nfail = nfail + 1
        end if
    end subroutine require
end program test_fortnum_bernstein
