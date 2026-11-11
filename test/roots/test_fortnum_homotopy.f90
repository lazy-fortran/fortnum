program test_fortnum_homotopy
    ! Behavioral tests for fortnum_homotopy.
    !
    ! Oracles: hand-written values and Jacobians of the test polynomials;
    ! closed-form roots of x^k = a; the four solutions of x^2 + y^2 = 5,
    ! x y = 2; the permutations of prescribed roots for the Vieta system
    ! e_k(x) = e_k(r) with r spanning twelve orders of magnitude, whose
    ! coefficients are exact sums of positive terms; the Bezout number of a
    ! dense system with an all-pairs distinctness check and residuals
    ! evaluated independently of the module.
    !$ use omp_lib, only: omp_set_num_threads, omp_get_max_threads
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_homotopy, only: poly_system_t, poly_system_init, &
        poly_system_add_term, poly_system_finalize, poly_system_eval, &
        homotopy_options_t, homotopy_stats_t, homotopy_track, homotopy_plant, &
        homotopy_monodromy, homotopy_move
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    integer, parameter :: NO_VARS(0) = [integer ::]
    integer :: nfail

    nfail = 0
    call test_eval(nfail)
    call test_plant(nfail)
    call test_univariate(nfail)
    call test_two_by_two(nfail)
    call test_vieta_ill_scaled(nfail)
    call test_dense_bezout(nfail)
    call test_thread_determinism(nfail)
    call test_endpoint_resolution(nfail)
    call test_ill_conditioned_point(nfail)
    call test_quad_fallback(nfail)
    call test_coinciding_endpoints(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine fail(label, nfail)
        character(*), intent(in) :: label
        integer, intent(inout) :: nfail

        nfail = nfail + 1
        write (error_unit, "(a,a,a)") "FAIL [", label, "]"
    end subroutine fail

    ! F1 = c1 x1^3 x2^2 x3 + c2 x2 + c3, F2 = c4 x1 x1 + c5 x3, F3 = c6 x2 x3 + c7,
    ! with the second term of F2 given as repeated indices.
    subroutine test_eval(nfail)
        integer, intent(inout) :: nfail

        type(poly_system_t) :: sys
        complex(dp) :: x(3), c(7), dc(7), f(3), jac(3, 3), ft(3)
        complex(dp) :: fe(3), je(3, 3), fte(3), t(7)
        real(dp) :: rs(3), rse(3)
        integer :: k

        call poly_system_init(sys, 3)
        call poly_system_add_term(sys, 1, [1, 2, 1, 3, 2, 1], 1)
        call poly_system_add_term(sys, 1, [2], 2)
        call poly_system_add_term(sys, 1, NO_VARS, 3)
        call poly_system_add_term(sys, 2, [1, 1], 4)
        call poly_system_add_term(sys, 2, [3], 5)
        call poly_system_add_term(sys, 3, [2, 3], 6)
        call poly_system_add_term(sys, 3, NO_VARS, 7)
        call poly_system_finalize(sys)
        if (sys%ncoef /= 7) call fail("eval: ncoef", nfail)

        x = [cmplx(0.7_dp, -0.4_dp, dp), cmplx(-1.3_dp, 0.2_dp, dp), &
            cmplx(0.5_dp, 0.9_dp, dp)]
        do k = 1, 7
            c(k) = cmplx(0.3_dp * k - 1.0_dp, 0.1_dp * k * k, dp)
            dc(k) = cmplx(1.0_dp / k, -0.5_dp, dp)
        end do
        call poly_system_eval(sys, x, c, dc, f, jac, ft, rs)

        t(1) = x(1)**3 * x(2)**2 * x(3)
        t(2) = x(2)
        t(3) = (1.0_dp, 0.0_dp)
        t(4) = x(1)**2
        t(5) = x(3)
        t(6) = x(2) * x(3)
        t(7) = (1.0_dp, 0.0_dp)
        fe(1) = c(1) * t(1) + c(2) * t(2) + c(3)
        fe(2) = c(4) * t(4) + c(5) * t(5)
        fe(3) = c(6) * t(6) + c(7)
        fte(1) = dc(1) * t(1) + dc(2) * t(2) + dc(3)
        fte(2) = dc(4) * t(4) + dc(5) * t(5)
        fte(3) = dc(6) * t(6) + dc(7)
        je = (0.0_dp, 0.0_dp)
        je(1, 1) = c(1) * 3.0_dp * x(1)**2 * x(2)**2 * x(3)
        je(1, 2) = c(1) * 2.0_dp * x(1)**3 * x(2) * x(3) + c(2)
        je(1, 3) = c(1) * x(1)**3 * x(2)**2
        je(2, 1) = c(4) * 2.0_dp * x(1)
        je(2, 3) = c(5)
        je(3, 2) = c(6) * x(3)
        je(3, 3) = c(6) * x(2)
        rse(1) = max(abs(c(1) * t(1)), abs(c(2) * t(2)), abs(c(3)))
        rse(2) = max(abs(c(4) * t(4)), abs(c(5) * t(5)))
        rse(3) = max(abs(c(6) * t(6)), abs(c(7)))
        if (maxval(abs(f - fe)) > 1.0e-14_dp) call fail("eval: f", nfail)
        if (maxval(abs(ft - fte)) > 1.0e-14_dp) call fail("eval: ft", nfail)
        if (maxval(abs(jac - je)) > 1.0e-14_dp) call fail("eval: jac", nfail)
        if (maxval(abs(rs - rse)) > 1.0e-14_dp) call fail("eval: rowscale", nfail)
    end subroutine test_eval

    subroutine test_plant(nfail)
        integer, intent(inout) :: nfail

        type(poly_system_t) :: sys
        complex(dp) :: x0(2), c0(5), f(2), jac(2, 2), ft(2)
        real(dp) :: rs(2)

        call two_by_two_system(sys)
        call homotopy_plant(sys, x0, c0, 7)
        call poly_system_eval(sys, x0, c0, c0, f, jac, ft, rs)
        if (maxval(abs(f) / rs) > 1.0e-14_dp) call fail("plant: residual", nfail)
        if (maxval(abs(abs(x0) - 1.0_dp)) > 1.0e-14_dp) &
            call fail("plant: x0 on the unit torus", nfail)
    end subroutine test_plant

    ! x^k = a: every computed root matches exactly one closed-form root.
    ! x = c and y = 1 + c with c moved from 1 to 1e-30: the endpoint x = 1e-30
    ! lies far below the resolution of t near 1, so it is reached only when
    ! the tracker resolves the remaining distance to the target.
    subroutine test_endpoint_resolution(nfail)
        integer, intent(inout) :: nfail

        type(poly_system_t) :: sys
        type(homotopy_options_t) :: opts
        complex(dp) :: x(2), ca(5), cb(5)
        real(dp) :: acc
        logical :: ok
        integer :: nsteps

        call poly_system_init(sys, 2)
        call poly_system_add_term(sys, 1, [1], 1)
        call poly_system_add_term(sys, 1, NO_VARS, 2)
        call poly_system_add_term(sys, 2, [2], 3)
        call poly_system_add_term(sys, 2, [1], 4)
        call poly_system_add_term(sys, 2, NO_VARS, 5)
        call poly_system_finalize(sys)
        ! x - c = 0 and y - x - 1 = 0
        ca = [(1.0_dp, 0.0_dp), (-1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), &
            (-1.0_dp, 0.0_dp), (-1.0_dp, 0.0_dp)]
        cb = ca
        cb(2) = cmplx(-1.0e-30_dp, 0.0_dp, dp)
        x = [(1.0_dp, 0.0_dp), (2.0_dp, 0.0_dp)]
        call homotopy_track(sys, x, ca, cb, opts, ok, nsteps, acc)
        if (.not. ok) then
            call fail("endpoint resolution: path failed", nfail)
        else if (abs(x(1) - 1.0e-30_dp) > 1.0e-12_dp * 1.0e-30_dp .or. &
                 abs(x(2) - 1.0_dp) > 1.0e-12_dp) then
            call fail("endpoint resolution: wrong endpoint", nfail)
        end if
    end subroutine test_endpoint_resolution

    ! x + y = 2, x + (1 + delta) y = 2 + delta with delta = 2**-30: exact
    ! coefficients, the unique solution (1, 1), and a condition number near
    ! 4e9. With the residual in double, Newton stalls near 5e-7 relative; the
    ! corrector must still reach (1, 1) to double accuracy.
    subroutine test_ill_conditioned_point(nfail)
        integer, intent(inout) :: nfail

        type(poly_system_t) :: sys
        type(homotopy_options_t) :: opts
        complex(dp) :: x(2), c(6)
        real(dp) :: acc, delta
        logical :: ok
        integer :: nsteps

        call poly_system_init(sys, 2)
        call poly_system_add_term(sys, 1, [1], 1)
        call poly_system_add_term(sys, 1, [2], 2)
        call poly_system_add_term(sys, 1, NO_VARS, 3)
        call poly_system_add_term(sys, 2, [1], 4)
        call poly_system_add_term(sys, 2, [2], 5)
        call poly_system_add_term(sys, 2, NO_VARS, 6)
        call poly_system_finalize(sys)
        delta = 2.0_dp**(-30)
        c = [(1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), (-2.0_dp, 0.0_dp), &
            (1.0_dp, 0.0_dp), cmplx(1.0_dp + delta, 0.0_dp, dp), &
            cmplx(-(2.0_dp + delta), 0.0_dp, dp)]
        x = [cmplx(1.0_dp + 1.0e-6_dp, 0.0_dp, dp), cmplx(1.0_dp - 1.0e-6_dp, 0.0_dp, dp)]
        call homotopy_track(sys, x, c, c, opts, ok, nsteps, acc)
        if (.not. ok) then
            call fail("ill-conditioned point: not converged", nfail)
        else if (maxval(abs(x - (1.0_dp, 0.0_dp))) > 1.0e-12_dp) then
            call fail("ill-conditioned point: inaccurate", nfail)
        end if
    end subroutine test_ill_conditioned_point

    ! From x**2 + y = 6, x**2 + 2 y = 8 to x**2 + y = 6,
    ! x**2 + y + delta y = 6 + delta with delta = 2**-60 in separate terms: the
    ! solution moves from (2, 2) to (sqrt 5, 1). In double the Jacobian row
    ! 2 x, 1 + delta rounds to the first row, while in quad precision the
    ! condition number is near 1e18, so quad Newton reaches the irrational
    ! endpoint only to about 1e-16. The path must reach (sqrt 5, 1) with the
    ! quad-precision fallback and fail without it.
    subroutine test_quad_fallback(nfail)
        integer, intent(inout) :: nfail

        type(poly_system_t) :: sys
        type(homotopy_options_t) :: opts
        complex(dp) :: sols(2, 1), out(2, 1), ca(8), cb(8), want(2)
        real(dp) :: acc(1), delta
        logical :: ok(1)

        call poly_system_init(sys, 2)
        call poly_system_add_term(sys, 1, [1, 1], 1)
        call poly_system_add_term(sys, 1, [2], 2)
        call poly_system_add_term(sys, 1, NO_VARS, 3)
        call poly_system_add_term(sys, 2, [1, 1], 4)
        call poly_system_add_term(sys, 2, [2], 5)
        call poly_system_add_term(sys, 2, [2], 6)
        call poly_system_add_term(sys, 2, NO_VARS, 7)
        call poly_system_add_term(sys, 2, NO_VARS, 8)
        call poly_system_finalize(sys)
        delta = 2.0_dp**(-60)
        ca = [(1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), (-6.0_dp, 0.0_dp), &
            (1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), &
            (-8.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
        cb = [(1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), (-6.0_dp, 0.0_dp), &
            (1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), cmplx(delta, 0.0_dp, dp), &
            (-6.0_dp, 0.0_dp), cmplx(-delta, 0.0_dp, dp)]
        sols(:, 1) = [(2.0_dp, 0.0_dp), (2.0_dp, 0.0_dp)]
        want = [cmplx(sqrt(5.0_dp), 0.0_dp, dp), (1.0_dp, 0.0_dp)]
        call homotopy_move(sys, sols, 1, ca, cb, opts, out, ok, acc)
        if (.not. ok(1)) then
            call fail("quad fallback: path failed", nfail)
        else if (maxval(abs(out(:, 1) - want)) > 1.0e-12_dp * 3.0_dp) then
            call fail("quad fallback: inaccurate endpoint", nfail)
        end if
        opts%quad_fallback = .false.
        call homotopy_move(sys, sols, 1, ca, cb, opts, out, ok, acc)
        if (ok(1)) then
            if (maxval(abs(out(:, 1) - want)) <= 1.0e-12_dp * 3.0_dp) &
                call fail("quad fallback: double alone already succeeds", nfail)
        end if
    end subroutine test_quad_fallback

    ! The same start point twice: both paths end at one root of x**2 = 2, so
    ! homotopy_move must report both endpoints as coinciding, and the other
    ! root once when it is started separately.
    subroutine test_coinciding_endpoints(nfail)
        integer, intent(inout) :: nfail

        type(poly_system_t) :: sys
        type(homotopy_options_t) :: opts
        type(homotopy_stats_t) :: stats
        complex(dp) :: sols(1, 3), out(1, 3), ca(2), cb(2)
        real(dp) :: acc(3)
        logical :: ok(3)

        call poly_system_init(sys, 1)
        call poly_system_add_term(sys, 1, [1, 1], 1)
        call poly_system_add_term(sys, 1, NO_VARS, 2)
        call poly_system_finalize(sys)
        ca = [(1.0_dp, 0.0_dp), (-1.0_dp, 0.0_dp)]
        cb = [(1.0_dp, 0.0_dp), (-2.0_dp, 0.0_dp)]
        sols(1, :) = [(1.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), (-1.0_dp, 0.0_dp)]
        call homotopy_move(sys, sols, 3, ca, cb, opts, out, ok, acc, stats)
        if (.not. all(ok)) call fail("coinciding endpoints: path failed", nfail)
        if (stats%jumps /= 2) call fail("coinciding endpoints: not reported", nfail)
        if (abs(out(1, 3) + sqrt(2.0_dp)) > 1.0e-12_dp) &
            call fail("coinciding endpoints: wrong third root", nfail)
    end subroutine test_coinciding_endpoints

    subroutine test_univariate(nfail)
        integer, intent(inout) :: nfail

        type(poly_system_t) :: sys
        type(homotopy_options_t) :: opts
        type(fortnum_status_t) :: status
        complex(dp), allocatable :: sols(:, :), out(:, :), roots(:)
        logical, allocatable :: ok(:), used(:)
        complex(dp) :: x0(1), c0(2), ct(2), a
        real(dp), parameter :: pi = 4.0_dp * atan(1.0_dp)
        character(len=32) :: label
        integer :: k, j, m, nsol, best
        integer :: ones(8)
        real(dp) :: err, berr

        ones = 1
        a = cmplx(2.0_dp, -1.5_dp, dp)
        do k = 1, 8
            write (label, "(a,i0)") "univariate k=", k
            call poly_system_init(sys, 1)
            call poly_system_add_term(sys, 1, ones(1:k), 1)
            call poly_system_add_term(sys, 1, NO_VARS, 2)
            call poly_system_finalize(sys)
            call homotopy_plant(sys, x0, c0, 100 + k)
            ! Known count: do not stop on loops that find nothing new.
            opts = homotopy_options_t(target=k, seed=k, stagnant_loops=50)
            call homotopy_monodromy(sys, x0, c0, opts, sols, nsol, status)
            if (status%code /= FORTNUM_OK .or. nsol /= k) then
                call fail(trim(label)//": count", nfail)
                cycle
            end if
            ct(1) = (1.0_dp, 0.0_dp)
            ct(2) = -a
            allocate (out(1, nsol), ok(nsol), roots(k), used(k))
            call homotopy_move(sys, sols, nsol, c0, ct, opts, out, ok)
            do m = 1, k
                roots(m) = abs(a)**(1.0_dp / k) &
                    * exp(cmplx(0.0_dp, (atan2(aimag(a), real(a)) &
                    + 2.0_dp * pi * (m - 1)) / k, dp))
            end do
            used = .false.
            do j = 1, nsol
                berr = huge(1.0_dp)
                best = 0
                do m = 1, k
                    err = abs(out(1, j) - roots(m)) / abs(roots(m))
                    if (err < berr) then
                        berr = err
                        best = m
                    end if
                end do
                if (.not. ok(j) .or. berr > 1.0e-13_dp) then
                    call fail(trim(label)//": root accuracy", nfail)
                else if (used(best)) then
                    call fail(trim(label)//": root found twice", nfail)
                else
                    used(best) = .true.
                end if
            end do
            deallocate (out, ok, roots, used)
        end do
    end subroutine test_univariate

    ! c1 x^2 + c2 y^2 + c3 = 0, c4 x y + c5 = 0
    subroutine two_by_two_system(sys)
        type(poly_system_t), intent(out) :: sys

        call poly_system_init(sys, 2)
        call poly_system_add_term(sys, 1, [1, 1], 1)
        call poly_system_add_term(sys, 1, [2, 2], 2)
        call poly_system_add_term(sys, 1, NO_VARS, 3)
        call poly_system_add_term(sys, 2, [1, 2], 4)
        call poly_system_add_term(sys, 2, NO_VARS, 5)
        call poly_system_finalize(sys)
    end subroutine two_by_two_system

    ! The target count is unknown here, so monodromy stops on stagnation.
    subroutine test_two_by_two(nfail)
        integer, intent(inout) :: nfail

        type(poly_system_t) :: sys
        type(homotopy_options_t) :: opts
        type(homotopy_stats_t) :: stats
        type(fortnum_status_t) :: status
        complex(dp), allocatable :: sols(:, :), out(:, :)
        logical, allocatable :: ok(:)
        real(dp), allocatable :: acc(:)
        complex(dp) :: x0(2), c0(5), ct(5), xt(2)
        real(dp) :: ref(2, 4)
        logical :: hit(4), track_ok
        integer :: nsol, j, m, nsteps

        call two_by_two_system(sys)
        call homotopy_plant(sys, x0, c0, 3)
        opts = homotopy_options_t(target=0, stagnant_loops=5, seed=11)
        call homotopy_monodromy(sys, x0, c0, opts, sols, nsol, status, stats)
        if (status%code /= FORTNUM_OK .or. nsol /= 4) then
            call fail("2x2: count", nfail)
            return
        end if
        if (stats%jumps /= 0 .or. stats%failed /= 0) call fail("2x2: paths", nfail)
        ct = [1.0_dp, 1.0_dp, -5.0_dp, 1.0_dp, -2.0_dp]
        allocate (out(2, nsol), ok(nsol), acc(nsol))
        call homotopy_move(sys, sols, nsol, c0, ct, opts, out, ok, acc)
        ref = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp, -1.0_dp, -2.0_dp, &
            -2.0_dp, -1.0_dp], [2, 4])
        hit = .false.
        do j = 1, nsol
            if (.not. ok(j) .or. acc(j) > 1.0e-12_dp) cycle
            do m = 1, 4
                if (maxval(abs(out(:, j) - ref(:, m))) < 1.0e-13_dp) hit(m) = .true.
            end do
        end do
        if (.not. all(hit)) call fail("2x2: solutions", nfail)

        ! The single-path tracker follows the first solution the same way.
        xt = sols(:, 1)
        call homotopy_track(sys, xt, c0, ct, opts, track_ok, nsteps)
        if (.not. track_ok .or. nsteps < 1) call fail("2x2: track status", nfail)
        if (maxval(abs(xt - out(:, 1))) > 1.0e-13_dp) call fail("2x2: track", nfail)
    end subroutine test_two_by_two

    ! Vieta system e_k(x) = e_k(r), k = 1..4, with r = 1, 1e4, 1e8, 1e12. Its
    ! 24 solutions are the permutations of r.
    subroutine test_vieta_ill_scaled(nfail)
        integer, intent(inout) :: nfail

        integer, parameter :: n = 4, nperm = 24
        type(poly_system_t) :: sys
        type(homotopy_options_t) :: opts
        type(homotopy_stats_t) :: stats
        type(fortnum_status_t) :: status
        complex(dp), allocatable :: sols(:, :), out(:, :), x0(:), c0(:), ct(:)
        logical, allocatable :: ok(:), used(:)
        real(dp) :: r(n), e(0:n), err, berr
        integer :: perms(n, nperm), k, j, m, nsol, ncoef, best, buf(n)

        r = [1.0_dp, 1.0e4_dp, 1.0e8_dp, 1.0e12_dp]
        ! e_k(r) by the product recurrence: only sums of positive terms.
        e = 0.0_dp
        e(0) = 1.0_dp
        do j = 1, n
            do k = j, 1, -1
                e(k) = e(k) + r(j) * e(k - 1)
            end do
        end do
        call poly_system_init(sys, n)
        ncoef = 0
        allocate (ct(64))
        do k = 1, n
            call add_subsets(sys, n, k, 1, buf, 0, ncoef, ct)
            ncoef = ncoef + 1
            call poly_system_add_term(sys, k, NO_VARS, ncoef)
            ct(ncoef) = -e(k)
        end do
        call poly_system_finalize(sys)
        allocate (x0(n), c0(ncoef))
        call homotopy_plant(sys, x0, c0, 5)
        opts = homotopy_options_t(target=nperm, seed=2, stagnant_loops=50)
        call homotopy_monodromy(sys, x0, c0, opts, sols, nsol, status, stats)
        if (status%code /= FORTNUM_OK .or. nsol /= nperm) then
            call fail("vieta: count", nfail)
            return
        end if
        allocate (out(n, nsol), ok(nsol), used(nperm))
        call homotopy_move(sys, sols, nsol, c0, ct(1:ncoef), opts, out, ok)
        call permutations(n, perms)
        used = .false.
        do j = 1, nsol
            berr = huge(1.0_dp)
            best = 0
            do m = 1, nperm
                err = maxval(abs(out(:, j) - r(perms(:, m))) / r(perms(:, m)))
                if (err < berr) then
                    berr = err
                    best = m
                end if
            end do
            if (.not. ok(j) .or. berr > 1.0e-12_dp) then
                write (error_unit, "(a,es10.2)") "vieta relative error ", berr
                call fail("vieta: accuracy", nfail)
            else if (used(best)) then
                call fail("vieta: permutation found twice", nfail)
            else
                used(best) = .true.
            end if
        end do

    end subroutine test_vieta_ill_scaled

    ! Add the product of every kleft-subset of first..n, appended to
    ! buf(1:depth), to equation depth + kleft with coefficient 1 in ct.
    recursive subroutine add_subsets(sys, n, kleft, first, buf, depth, ncoef, ct)
        type(poly_system_t), intent(inout) :: sys
        integer, intent(in) :: n, kleft, first, depth
        integer, intent(inout) :: buf(:)
        integer, intent(inout) :: ncoef
        complex(dp), intent(inout) :: ct(:)

        integer :: i

        if (kleft == 0) then
            ncoef = ncoef + 1
            call poly_system_add_term(sys, depth, buf(1:depth), ncoef)
            ct(ncoef) = (1.0_dp, 0.0_dp)
            return
        end if
        do i = first, n - kleft + 1
            buf(depth + 1) = i
            call add_subsets(sys, n, kleft - 1, i + 1, buf, depth + 1, ncoef, ct)
        end do
    end subroutine add_subsets

    subroutine permutations(n, p)
        integer, intent(in) :: n
        integer, intent(out) :: p(:, :)

        integer :: a(n), i, j, cnt

        cnt = 0
        do i = 0, n**n - 1
            do j = 1, n
                a(j) = mod(i / n**(j - 1), n) + 1
            end do
            if (distinct(a)) then
                cnt = cnt + 1
                p(:, cnt) = a
            end if
        end do
    end subroutine permutations

    logical function distinct(a)
        integer, intent(in) :: a(:)

        integer :: i

        distinct = .true.
        do i = 2, size(a)
            if (any(a(1:i - 1) == a(i))) distinct = .false.
        end do
    end function distinct

    ! Dense system in three unknowns with degrees 2, 3, 3: every monomial up to
    ! the degree has its own coefficient, so the generic count is the Bezout
    ! number 18.
    subroutine dense_system(sys)
        type(poly_system_t), intent(out) :: sys

        integer :: eq, deg(3), a, b, cexp, ncoef, i, v(3)

        deg = [2, 3, 3]
        call poly_system_init(sys, 3)
        ncoef = 0
        do eq = 1, 3
            do a = 0, deg(eq)
                do b = 0, deg(eq) - a
                    do cexp = 0, deg(eq) - a - b
                        ncoef = ncoef + 1
                        ! Monomial x1^a x2^b x3^cexp as indices with multiplicity.
                        do i = 1, a + b + cexp
                            v(i) = merge(1, merge(2, 3, i <= a + b), i <= a)
                        end do
                        call poly_system_add_term(sys, eq, v(1:a + b + cexp), ncoef)
                    end do
                end do
            end do
        end do
        call poly_system_finalize(sys)
    end subroutine dense_system

    subroutine test_dense_bezout(nfail)
        integer, intent(inout) :: nfail

        type(poly_system_t) :: sys
        type(homotopy_options_t) :: opts
        type(homotopy_stats_t) :: stats
        type(fortnum_status_t) :: status
        complex(dp), allocatable :: sols(:, :), x0(:), c0(:)
        integer :: nsol, i, j

        call dense_system(sys)
        allocate (x0(3), c0(sys%ncoef))
        call homotopy_plant(sys, x0, c0, 9)
        opts = homotopy_options_t(target=0, stagnant_loops=4, seed=4)
        call homotopy_monodromy(sys, x0, c0, opts, sols, nsol, status, stats)
        if (status%code /= FORTNUM_OK .or. nsol /= 18) then
            write (error_unit, "(a,i0)") "dense: found ", nsol
            call fail("dense: Bezout count", nfail)
        end if
        if (stats%jumps /= 0) call fail("dense: path jumps", nfail)
        do j = 1, nsol
            if (dense_residual(sols(:, j), c0) > 1.0e-12_dp) &
                call fail("dense: residual", nfail)
            do i = 1, j - 1
                if (maxval(abs(sols(:, i) - sols(:, j)) / abs(sols(:, j))) &
                    < 1.0e-4_dp) call fail("dense: duplicate solution", nfail)
            end do
        end do
    end subroutine test_dense_bezout

    ! max_i |F_i| / max_k |c_k m_k| from the same monomial ordering as
    ! dense_system, evaluated here without the module.
    real(dp) function dense_residual(x, c) result(res)
        complex(dp), intent(in) :: x(3), c(:)

        integer :: eq, deg(3), a, b, cexp, k
        complex(dp) :: f, term
        real(dp) :: scale

        deg = [2, 3, 3]
        k = 0
        res = 0.0_dp
        do eq = 1, 3
            f = (0.0_dp, 0.0_dp)
            scale = 0.0_dp
            do a = 0, deg(eq)
                do b = 0, deg(eq) - a
                    do cexp = 0, deg(eq) - a - b
                        k = k + 1
                        term = c(k) * x(1)**a * x(2)**b * x(3)**cexp
                        f = f + term
                        scale = max(scale, abs(term))
                    end do
                end do
            end do
            res = max(res, abs(f) / scale)
        end do
    end function dense_residual

    ! Path results do not depend on the thread count, and the deduplication
    ! is serial, so 1 and 4 threads give identical solution sets.
    subroutine test_thread_determinism(nfail)
        integer, intent(inout) :: nfail

        type(poly_system_t) :: sys
        type(homotopy_options_t) :: opts
        type(homotopy_stats_t) :: s1, s4
        type(fortnum_status_t) :: status
        complex(dp), allocatable :: sols1(:, :), sols4(:, :), x0(:), c0(:)
        integer :: n1, n4
        !$ integer :: nthreads

        call dense_system(sys)
        allocate (x0(3), c0(sys%ncoef))
        call homotopy_plant(sys, x0, c0, 21)
        opts = homotopy_options_t(target=0, stagnant_loops=3, seed=8)
        !$ nthreads = omp_get_max_threads()
        !$ call omp_set_num_threads(1)
        call homotopy_monodromy(sys, x0, c0, opts, sols1, n1, status, s1)
        !$ call omp_set_num_threads(4)
        call homotopy_monodromy(sys, x0, c0, opts, sols4, n4, status, s4)
        !$ call omp_set_num_threads(nthreads)
        if (n1 /= n4 .or. s1%tracks /= s4%tracks .or. s1%steps /= s4%steps) then
            call fail("threads: counts", nfail)
            return
        end if
        if (maxval(abs(sols1 - sols4)) > 0.0_dp) call fail("threads: solutions", nfail)
    end subroutine test_thread_determinism

end program test_fortnum_homotopy
