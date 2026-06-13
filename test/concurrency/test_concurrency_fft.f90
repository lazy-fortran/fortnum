program test_concurrency_fft
    ! M5.1 concurrency: fft (per-thread plan).
    ! A fortnum_fft_plan_t holds read-only twiddle/chirp tables once built; the
    ! transform routines take the plan intent(in) and write a caller-owned
    ! output array. Thread safety rests on each thread owning its own plan and
    ! its own input/output buffers. This test transforms a batch of distinct
    ! signals serially, then transforms the same batch with one private plan per
    ! thread, and asserts bit-for-bit equality of every spectrum.
    !
    ! PRIMAL concurrency only; derivative-product fft concurrency lands in M6
    ! (#40).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_fft, only: fft_r2c, fft_c2c, fortnum_fft_plan_t, fft_plan_init
    implicit none

    integer, parameter :: n = 96            ! mixed-radix length
    integer, parameter :: nbatch = 512
    integer :: nfail, j, k
    complex(dp) :: ser_spec(n/2 + 1, nbatch), par_spec(n/2 + 1, nbatch)
    complex(dp) :: ser_c2c(n, nbatch), par_c2c(n, nbatch)
    real(dp) :: signals(n, nbatch)

    nfail = 0

    do j = 1, nbatch
        do k = 1, n
            signals(k, j) = sin(0.13_dp*k*j) + 0.5_dp*cos(0.07_dp*k)
        end do
    end do

    call serial_r2c()
    call parallel_r2c()
    call check_spec(ser_spec, par_spec, "fft_r2c parallel == serial")

    call serial_c2c()
    call parallel_c2c()
    call check_c2c(ser_c2c, par_c2c, "fft_c2c parallel == serial")

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_concurrency_fft: all tests passed"

contains

    subroutine serial_r2c()
        type(fortnum_fft_plan_t) :: plan
        complex(dp) :: c(n/2 + 1)
        call fft_plan_init(plan, n)
        do j = 1, nbatch
            call fft_r2c(signals(:, j), c, plan)
            ser_spec(:, j) = c
        end do
    end subroutine serial_r2c

    subroutine parallel_r2c()
        type(fortnum_fft_plan_t) :: plan
        complex(dp) :: c(n/2 + 1)
        ! Each thread builds and owns one private plan for the whole region;
        ! no thread reads another's plan or output column.
        !$omp parallel default(shared) private(plan, c, j)
        call fft_plan_init(plan, n)
        !$omp do schedule(static)
        do j = 1, nbatch
            call fft_r2c(signals(:, j), c, plan)
            par_spec(:, j) = c
        end do
        !$omp end do
        !$omp end parallel
    end subroutine parallel_r2c

    subroutine serial_c2c()
        complex(dp) :: z(n)
        do j = 1, nbatch
            z = cmplx(signals(:, j), 0.0_dp, dp)
            call fft_c2c(z, -1)
            ser_c2c(:, j) = z
        end do
    end subroutine serial_c2c

    subroutine parallel_c2c()
        complex(dp) :: z(n)
        ! fft_c2c builds its plan internally on the call stack each call, so the
        ! only shared state is the output array, written by disjoint columns.
        !$omp parallel do default(shared) private(z, j) schedule(static)
        do j = 1, nbatch
            z = cmplx(signals(:, j), 0.0_dp, dp)
            call fft_c2c(z, -1)
            par_c2c(:, j) = z
        end do
        !$omp end parallel do
    end subroutine parallel_c2c

    subroutine check_spec(ref, got, name)
        complex(dp),  intent(in) :: ref(:,:), got(:,:)
        character(*), intent(in) :: name
        integer :: a, b
        do b = 1, size(ref, 2)
            do a = 1, size(ref, 1)
                if (ref(a, b) /= got(a, b)) then
                    nfail = nfail + 1
                    write (error_unit, "(a,a,a,i0,a,i0)") "FAIL: ", name, &
                        " at (", a, ",", b
                    return
                end if
            end do
        end do
    end subroutine check_spec

    subroutine check_c2c(ref, got, name)
        complex(dp),  intent(in) :: ref(:,:), got(:,:)
        character(*), intent(in) :: name
        integer :: a, b
        do b = 1, size(ref, 2)
            do a = 1, size(ref, 1)
                if (ref(a, b) /= got(a, b)) then
                    nfail = nfail + 1
                    write (error_unit, "(a,a)") "FAIL: ", name
                    return
                end if
            end do
        end do
    end subroutine check_c2c

end program test_concurrency_fft
