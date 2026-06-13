program test_concurrency_roots_interp
    ! M5.1 concurrency: roots (bisect/newton/brent) and interp/polynomial.
    ! These routines hold no caller-owned mutable state across calls and write
    ! only caller-owned scalars/arrays; thread safety rests on each thread
    ! owning its output. This test solves a family of bracketed problems and
    ! evaluates interpolation weights serially, then sharded across OpenMP
    ! threads, asserting bit-for-bit equality.
    !
    ! PRIMAL concurrency only; derivative-product concurrency lands in M6 (#40).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t
    use fortnum_roots, only: root_bisect, root_newton, root_brent
    use fortnum_interp, only: grid_search
    use fortnum_polynomial, only: lagrange_weights, lagrange_deriv_weights
    implicit none

    ! root_fn_t / root_fn_df_t take no ctx, so the problem cannot be varied per
    ! iteration through a parameter; every iteration solves the same fixed
    ! bracketed root. The concurrency property under test is that many threads
    ! calling the solver concurrently each get the identical, race-free result,
    ! bit-for-bit equal to the serial run. The bracket varies per index so the
    ! solver path is not trivially constant, while the root stays well defined.
    integer, parameter :: nprob = 512
    integer :: nfail, j
    real(dp) :: lo(nprob)
    real(dp) :: ser_b(nprob), par_b(nprob)
    real(dp) :: ser_n(nprob), par_n(nprob)
    real(dp) :: ser_r(nprob), par_r(nprob)

    nfail = 0
    ! Lower bracket endpoint in [0, 1); root of x*x-2 is sqrt(2) ~ 1.414, so
    ! every bracket [lo, 10] still contains it.
    do j = 1, nprob
        lo(j) = real(j - 1, dp) / real(nprob, dp)
    end do

    call serial_bisect()
    call parallel_bisect()
    call check_exact(ser_b, par_b, "root_bisect parallel == serial")

    call serial_newton()
    call parallel_newton()
    call check_exact(ser_n, par_n, "root_newton parallel == serial")

    call serial_brent()
    call parallel_brent()
    call check_exact(ser_r, par_r, "root_brent parallel == serial")

    call check_grid_search()
    call check_lagrange()

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_concurrency_roots_interp: all tests passed"

contains

    ! Fixed problem: root of x*x-2 = 0 at sqrt(2). No ctx, no shared state.
    pure function f_sq(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = x*x - 2.0_dp
    end function f_sq

    pure subroutine fdf_sq(x, fx, dfx)
        real(dp), intent(in)  :: x
        real(dp), intent(out) :: fx, dfx
        fx  = x*x - 2.0_dp
        dfx = 2.0_dp*x
    end subroutine fdf_sq

    subroutine serial_bisect()
        type(fortnum_status_t) :: st
        real(dp) :: x
        do j = 1, nprob
            call root_bisect(f_sq, lo(j), 10.0_dp, x, st, xtol=1.0e-12_dp)
            ser_b(j) = x
        end do
    end subroutine serial_bisect

    subroutine parallel_bisect()
        type(fortnum_status_t) :: st
        real(dp) :: x
        !$omp parallel do default(shared) private(st, x, j) schedule(static)
        do j = 1, nprob
            call root_bisect(f_sq, lo(j), 10.0_dp, x, st, xtol=1.0e-12_dp)
            par_b(j) = x
        end do
        !$omp end parallel do
    end subroutine parallel_bisect

    subroutine serial_newton()
        type(fortnum_status_t) :: st
        real(dp) :: x
        do j = 1, nprob
            call root_newton(fdf_sq, lo(j), 10.0_dp, 5.0_dp, x, st, &
                             xtol=1.0e-12_dp)
            ser_n(j) = x
        end do
    end subroutine serial_newton

    subroutine parallel_newton()
        type(fortnum_status_t) :: st
        real(dp) :: x
        !$omp parallel do default(shared) private(st, x, j) schedule(static)
        do j = 1, nprob
            call root_newton(fdf_sq, lo(j), 10.0_dp, 5.0_dp, x, st, &
                             xtol=1.0e-12_dp)
            par_n(j) = x
        end do
        !$omp end parallel do
    end subroutine parallel_newton

    subroutine serial_brent()
        type(fortnum_status_t) :: st
        real(dp) :: x
        do j = 1, nprob
            call root_brent(f_sq, lo(j), 10.0_dp, x, st, xtol=1.0e-12_dp)
            ser_r(j) = x
        end do
    end subroutine serial_brent

    subroutine parallel_brent()
        type(fortnum_status_t) :: st
        real(dp) :: x
        !$omp parallel do default(shared) private(st, x, j) schedule(static)
        do j = 1, nprob
            call root_brent(f_sq, lo(j), 10.0_dp, x, st, xtol=1.0e-12_dp)
            par_r(j) = x
        end do
        !$omp end parallel do
    end subroutine parallel_brent

    subroutine check_grid_search()
        integer, parameter :: ng = 64
        real(dp) :: grid(ng)
        integer  :: ser_idx(nprob), par_idx(nprob), idx, k
        real(dp) :: xi
        do k = 1, ng
            grid(k) = real(k, dp)
        end do
        do j = 1, nprob
            xi = 1.0_dp + 62.0_dp * real(j - 1, dp) / real(nprob, dp)
            call grid_search(grid, 1, ng, xi, idx)
            ser_idx(j) = idx
        end do
        !$omp parallel do default(shared) private(idx, xi, j) schedule(static)
        do j = 1, nprob
            xi = 1.0_dp + 62.0_dp * real(j - 1, dp) / real(nprob, dp)
            call grid_search(grid, 1, ng, xi, idx)
            par_idx(j) = idx
        end do
        !$omp end parallel do
        do k = 1, nprob
            if (ser_idx(k) /= par_idx(k)) then
                nfail = nfail + 1
                write (error_unit, "(a)") "FAIL: grid_search parallel == serial"
                exit
            end if
        end do
    end subroutine check_grid_search

    subroutine check_lagrange()
        integer, parameter :: np = 8
        real(dp) :: xp(np), ser_c(np, nprob), par_c(np, nprob)
        real(dp) :: ser_d(np, nprob), par_d(np, nprob)
        real(dp) :: coef(np), dcoef(np), x
        integer  :: k
        do k = 1, np
            xp(k) = real(k, dp)
        end do
        do j = 1, nprob
            x = 1.0_dp + 7.0_dp * real(j - 1, dp) / real(nprob, dp)
            call lagrange_weights(np, x, xp, coef)
            call lagrange_deriv_weights(np, x, xp, dcoef)
            ser_c(:, j) = coef
            ser_d(:, j) = dcoef
        end do
        !$omp parallel do default(shared) private(coef, dcoef, x, j) schedule(static)
        do j = 1, nprob
            x = 1.0_dp + 7.0_dp * real(j - 1, dp) / real(nprob, dp)
            call lagrange_weights(np, x, xp, coef)
            call lagrange_deriv_weights(np, x, xp, dcoef)
            par_c(:, j) = coef
            par_d(:, j) = dcoef
        end do
        !$omp end parallel do
        if (any(ser_c /= par_c)) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: lagrange_weights parallel == serial"
        end if
        if (any(ser_d /= par_d)) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: lagrange_deriv_weights parallel == serial"
        end if
    end subroutine check_lagrange

    subroutine check_exact(ref, got, name)
        real(dp),     intent(in) :: ref(:), got(:)
        character(*), intent(in) :: name
        integer :: k
        do k = 1, size(ref)
            if (ref(k) /= got(k)) then
                nfail = nfail + 1
                write (error_unit, "(a,a,a,i0)") "FAIL: ", name, " at index ", k
                return
            end if
        end do
    end subroutine check_exact

end program test_concurrency_roots_interp
