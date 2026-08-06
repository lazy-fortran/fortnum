program fortnum_bench_toeplitz
    use, intrinsic :: iso_fortran_env, only: int64
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_toeplitz, only: toeplitz_operator_t
    implicit none

    write (*, '(a)') &
        "case n rhs repeats toeplitz_seconds dense_seconds speedup checksum"
    call run_case(128, 4, 400, .false.)
    call run_case(256, 4, 240, .false.)
    call run_case(512, 4, 120, .true.)
    call run_case(1024, 4, 60, .false.)
    call run_case(2048, 4, 24, .false.)
    call run_case(4096, 4, 12, .false.)

contains

    subroutine run_case(n, n_rhs, repetitions, include_dense)
        integer, intent(in) :: n, n_rhs, repetitions
        logical, intent(in) :: include_dense

        type(toeplitz_operator_t) :: toeplitz
        type(fortnum_status_t) :: status
        real(dp), allocatable :: column(:), input(:, :), output(:, :)
        real(dp), allocatable :: dense(:, :)
        real(dp) :: toeplitz_seconds, dense_seconds, speedup, checksum
        integer(int64) :: start_count, stop_count, count_rate
        integer :: i, column_index

        allocate(column(n), input(n, n_rhs), output(n, n_rhs))
        do i = 1, n
            column(i) = exp(-0.015_dp*real(i - 1, dp))
        end do
        do column_index = 1, n_rhs
            do i = 1, n
                input(i, column_index) = sin(0.013_dp*real(i, dp)) + &
                    0.1_dp*real(column_index, dp)
            end do
        end do

        call toeplitz%initialize(column, status)
        if (status%code /= FORTNUM_OK) error stop "Toeplitz setup failed"
        call toeplitz%matmat(input, output)
        call system_clock(count_rate=count_rate)
        call system_clock(start_count)
        do i = 1, repetitions
            call toeplitz%matmat(input, output)
        end do
        call system_clock(stop_count)
        toeplitz_seconds = real(stop_count - start_count, dp)/ &
            real(count_rate, dp)/real(repetitions, dp)
        checksum = sum(output)
        dense_seconds = -1.0_dp
        speedup = -1.0_dp

        if (include_dense) then
            allocate(dense(n, n))
            call assemble_dense(column, dense)
            output = matmul(dense, input)
            call system_clock(start_count)
            do i = 1, repetitions
                output = matmul(dense, input)
            end do
            call system_clock(stop_count)
            dense_seconds = real(stop_count - start_count, dp)/ &
                real(count_rate, dp)/real(repetitions, dp)
            speedup = dense_seconds/toeplitz_seconds
        end if

        write (*, '(a,1x,3(i0,1x),4(es14.6,1x))') "toeplitz", n, n_rhs, &
            repetitions, toeplitz_seconds, dense_seconds, speedup, checksum
    end subroutine run_case

    subroutine assemble_dense(column, dense)
        real(dp), intent(in) :: column(:)
        real(dp), intent(out) :: dense(:, :)
        integer :: i, j

        do i = 1, size(column)
            do j = 1, size(column)
                dense(i, j) = column(abs(i - j) + 1)
            end do
        end do
    end subroutine assemble_dense

end program fortnum_bench_toeplitz
