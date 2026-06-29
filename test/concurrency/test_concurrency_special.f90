program test_concurrency_special
    ! M5.1 concurrency: special functions (bessel/dawson/gamma).
    ! The special functions are elemental/pure with no caller-owned state, so
    ! they are trivially thread-safe; this test still proves it by evaluating
    ! each on a large grid serially, then re-evaluating the same grid sharded
    ! across OpenMP threads, and asserting bit-for-bit equality.
    !
    ! PRIMAL concurrency only. Derivative-product concurrency (jvp/vjp/grad)
    ! arrives with the derivative entry points in M6 (#40); no derivative path
    ! exists yet, so nothing to exercise here beyond this note.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special, only: bessel_in, bessel_in_array, bessel_kn, &
        dawson, gamma_lower, gamma_reg_p
    implicit none

    integer, parameter :: ngrid = 4096
    integer :: nfail
    real(dp) :: x(ngrid)
    real(dp) :: ser_bessel(ngrid), ser_dawson(ngrid)
    real(dp) :: ser_gl(ngrid), ser_gp(ngrid)
    real(dp) :: par_bessel(ngrid), par_dawson(ngrid)
    real(dp) :: par_gl(ngrid), par_gp(ngrid)
    integer :: i

    nfail = 0

    ! Grid in (0, 8); gamma_lower/reg_p require a>0, x>=0.
    do i = 1, ngrid
        x(i) = 8.0_dp * real(i, dp) / real(ngrid + 1, dp)
    end do

    ! Serial reference.
    do i = 1, ngrid
        ser_bessel(i) = bessel_in(3, x(i))
        ser_dawson(i) = dawson(x(i))
        ser_gl(i)     = gamma_lower(2.5_dp, x(i))
        ser_gp(i)     = gamma_reg_p(2.5_dp, x(i))
    end do

    ! Parallel: each iteration writes its own index, no shared mutable state.
    !$omp parallel do default(shared) private(i) schedule(static)
    do i = 1, ngrid
        par_bessel(i) = bessel_in(3, x(i))
        par_dawson(i) = dawson(x(i))
        par_gl(i)     = gamma_lower(2.5_dp, x(i))
        par_gp(i)     = gamma_reg_p(2.5_dp, x(i))
    end do
    !$omp end parallel do

    call check_exact(ser_bessel, par_bessel, "bessel_in parallel == serial")
    call check_exact(ser_dawson, par_dawson, "dawson parallel == serial")
    call check_exact(ser_gl, par_gl, "gamma_lower parallel == serial")
    call check_exact(ser_gp, par_gp, "gamma_reg_p parallel == serial")

    call check_bessel_array_threadsafe()
    call check_kn_threadsafe()

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_concurrency_special: all tests passed"

contains

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

    ! bessel_in_array writes a caller-owned values(0:nmax). Each thread must own
    ! its own buffer; a shared buffer would race. Reduce per-thread sums and
    ! compare to the serial sum bit-for-bit.
    subroutine check_bessel_array_threadsafe()
        integer, parameter :: nmax = 6
        real(dp) :: ser_sum, par_sum, vals(0:nmax)
        integer  :: k
        ser_sum = 0.0_dp
        do i = 1, ngrid
            call bessel_in_array(nmax, x(i), vals)
            do k = 0, nmax
                ser_sum = ser_sum + vals(k)
            end do
        end do
        par_sum = 0.0_dp
        !$omp parallel do default(shared) private(i, k, vals) &
        !$omp   reduction(+:par_sum) schedule(static)
        do i = 1, ngrid
            call bessel_in_array(nmax, x(i), vals)
            do k = 0, nmax
                par_sum = par_sum + vals(k)
            end do
        end do
        !$omp end parallel do
        ! Floating reduction reorders adds, so compare at a tight tolerance.
        if (abs(par_sum - ser_sum) > 1.0e-10_dp * abs(ser_sum)) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: bessel_in_array reduction mismatch"
        end if
    end subroutine check_bessel_array_threadsafe

    subroutine check_kn_threadsafe()
        real(dp) :: ser(ngrid), par(ngrid)
        do i = 1, ngrid
            ser(i) = bessel_kn(2, x(i))
        end do
        !$omp parallel do default(shared) private(i) schedule(static)
        do i = 1, ngrid
            par(i) = bessel_kn(2, x(i))
        end do
        !$omp end parallel do
        call check_exact(ser, par, "bessel_kn parallel == serial")
    end subroutine check_kn_threadsafe

end program test_concurrency_special
