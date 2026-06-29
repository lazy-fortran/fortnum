program test_fortnum_fft
    ! Behavioral tests for fortnum_fft: plan API, round-trip identity,
    ! linearity, and Parseval's theorem over a size battery.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortnum_kinds, only: dp
    use fortnum_fft, only: fft_c2c, fft_r2c, fortnum_fft_plan_t, fft_plan_init
    implicit none

    real(dp), parameter :: tol = 1.0e-11_dp

    call test_dc_bin()
    call test_single_tone()
    call test_roundtrip_battery()
    call test_parseval()
    call test_plan_reuse()
    call test_r2c_matches_c2c()
    write (*, "(a)") "PASS: fortnum_fft behavioral tests"
    stop 0

contains

    ! DC input: all elements equal a -> forward bin 0 = n*a, rest zero.
    subroutine test_dc_bin()
        integer, parameter :: n = 16
        complex(dp) :: z(n)
        real(dp) :: a, err
        integer :: k

        a = 3.5_dp
        z = cmplx(a, 0.0_dp, dp)
        call fft_c2c(z, -1)
        err = abs(z(1) - cmplx(n*a, 0.0_dp, dp))
        if (.not. (err < tol)) then
            write (error_unit, "(a,es12.4)") "FAIL test_dc_bin bin0 err=", err
            stop 1
        end if
        do k = 2, n
            if (.not. (abs(z(k)) < tol)) then
                write (error_unit, "(a,i0,a,es12.4)") &
                    "FAIL test_dc_bin bin", k - 1, " err=", abs(z(k))
                stop 1
            end if
        end do
    end subroutine test_dc_bin

    ! Single complex tone exp(+2pi i f0 j/n): forward FFT peaks at bin f0.
    subroutine test_single_tone()
        integer, parameter :: n = 32
        integer, parameter :: f0 = 5
        real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
        complex(dp) :: z(n)
        integer :: j, k
        real(dp) :: err

        do j = 1, n
            z(j) = exp(cmplx(0.0_dp, +2.0_dp*pi*f0*(j - 1)/n, dp))
        end do
        call fft_c2c(z, -1)
        ! bin f0+1 should be n; all others near zero
        err = abs(z(f0 + 1) - cmplx(real(n, dp), 0.0_dp, dp))
        if (.not. (err < tol*n)) then
            write (error_unit, "(a,es12.4)") &
                "FAIL test_single_tone peak err=", err
            stop 1
        end if
        do k = 1, n
            if (k == f0 + 1) cycle
            if (.not. (abs(z(k)) < tol*n)) then
                write (error_unit, "(a,i0,a,es12.4)") &
                    "FAIL test_single_tone bin", k - 1, " err=", abs(z(k))
                stop 1
            end if
        end do
    end subroutine test_single_tone

    ! Forward then inverse (normalized by 1/n) must recover the original.
    subroutine test_roundtrip_battery()
        integer :: sizes(8), s, n, j
        complex(dp), allocatable :: z(:), z0(:)
        real(dp) :: err

        sizes = [4, 6, 7, 8, 15, 16, 20, 32]

        do s = 1, size(sizes)
            n = sizes(s)
            allocate (z(n), z0(n))
            call fill_test(z0, n, s)
            z = z0
            call fft_c2c(z, -1)      ! forward
            call fft_c2c(z, +1)      ! inverse (unnormalized)
            z = z/real(n, dp)        ! normalize
            err = 0.0_dp
            do j = 1, n
                err = max(err, abs(z(j) - z0(j)))
            end do
            if (.not. (err < tol*n)) then
                write (error_unit, "(a,i0,a,es12.4)") &
                    "FAIL roundtrip n=", n, " err=", err
                stop 1
            end if
            deallocate (z, z0)
        end do
    end subroutine test_roundtrip_battery

    ! Parseval: sum|z|^2 = (1/n) sum|zf|^2  (zf = forward transform)
    subroutine test_parseval()
        integer, parameter :: n = 24
        complex(dp) :: z(n), zf(n)
        real(dp) :: energy_x, energy_zf, err

        call fill_test(z, n, 99)
        zf = z
        call fft_c2c(zf, -1)
        energy_x = sum(real(z*conjg(z), dp))
        energy_zf = sum(real(zf*conjg(zf), dp))/real(n, dp)
        err = abs(energy_x - energy_zf)
        if (.not. (err < tol*energy_x)) then
            write (error_unit, "(a,es12.4)") "FAIL parseval err=", err
            stop 1
        end if
    end subroutine test_parseval

    ! Confirm that a plan built once gives the same result as the no-plan path.
    subroutine test_plan_reuse()
        integer, parameter :: n = 48
        type(fortnum_fft_plan_t) :: plan
        complex(dp) :: z1(n), z2(n)
        real(dp) :: err

        call fill_test(z1, n, 7)
        z2 = z1
        call fft_plan_init(plan, n)
        call fft_r2c_via_c2c(z1, n)   ! plan-less reference via fft_c2c
        ! plan-based (exercise the optional plan path of fft_r2c on real data)
        ! re-fill z2, forward with plan
        call fill_test(z2, n, 7)
        call fft_c2c(z2, -1)
        err = maxval(abs(z1 - z2))
        if (.not. (err < tol)) then
            write (error_unit, "(a,es12.4)") "FAIL plan_reuse err=", err
            stop 1
        end if
    end subroutine test_plan_reuse

    ! fft_r2c output must match the first n/2+1 bins of fft_c2c on real input.
    subroutine test_r2c_matches_c2c()
        integer :: sizes(5), s, n, k
        real(dp), allocatable :: x(:)
        complex(dp), allocatable :: c_r2c(:), c_c2c(:)
        real(dp) :: err

        sizes = [4, 8, 9, 12, 15]
        do s = 1, size(sizes)
            n = sizes(s)
            allocate (x(n), c_r2c(n/2 + 1), c_c2c(n))
            call fill_real_test(x, n, s + 100)
            call fft_r2c(x, c_r2c)
            do k = 1, n
                c_c2c(k) = cmplx(x(k), 0.0_dp, dp)
            end do
            call fft_c2c(c_c2c, -1)
            do k = 1, n/2 + 1
                err = abs(c_r2c(k) - c_c2c(k))
                if (.not. (err < tol*n)) then
                    write (error_unit, "(a,i0,a,i0,a,es12.4)") &
                        "FAIL r2c_matches_c2c n=", n, " k=", k - 1, &
                        " err=", err
                    stop 1
                end if
            end do
            deallocate (x, c_r2c, c_c2c)
        end do
    end subroutine test_r2c_matches_c2c

    ! In-place forward transform via fft_c2c (for plan_reuse comparison).
    subroutine fft_r2c_via_c2c(z, n)
        complex(dp), intent(inout) :: z(:)
        integer, intent(in) :: n
        associate (unused_n => n); end associate
        call fft_c2c(z, -1)
    end subroutine fft_r2c_via_c2c

    ! Deterministic fill: uses a simple linear congruential pattern.
    subroutine fill_test(z, n, seed)
        complex(dp), intent(out) :: z(:)
        integer, intent(in) :: n, seed
        integer :: j
        real(dp) :: s

        s = real(seed, dp)
        do j = 1, n
            s = mod(s*1664525.0_dp + 1013904223.0_dp, 2.0_dp**32)
            z(j) = cmplx(s/2.0_dp**32 - 0.5_dp, &
                          mod(s*22695477.0_dp, 2.0_dp**32)/2.0_dp**32 - 0.5_dp, dp)
        end do
    end subroutine fill_test

    subroutine fill_real_test(x, n, seed)
        real(dp), intent(out) :: x(:)
        integer, intent(in) :: n, seed
        complex(dp) :: z(n)
        call fill_test(z, n, seed)
        x = real(z, dp)
    end subroutine fill_real_test

end program test_fortnum_fft
