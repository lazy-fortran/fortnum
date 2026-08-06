program fortnum_bench_tensor_product_device
    use, intrinsic :: iso_fortran_env, only: int64
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use fortnum_tensor_product, only: tensor_factor_t, &
        tensor_product_operator_t
    implicit none

    write (*, '(a)') &
        "case n rhs repeats host_seconds device_transfer_seconds "// &
        "device_resident_seconds host_to_resident transfer_to_resident checksum"
    call run_case(8, 8, 8, 4, 200, .true.)
    call run_case(16, 16, 16, 4, 60, .false.)

contains

    subroutine run_case(n1, n2, n3, n_rhs, repetitions, include_dense)
        integer, intent(in) :: n1, n2, n3, n_rhs, repetitions
        logical, intent(in) :: include_dense

        type(tensor_factor_t) :: factors(3)
        type(tensor_product_operator_t) :: tensor_operator
        type(fortnum_status_t) :: status
        real(dp), allocatable :: input(:, :), output(:, :), dense(:, :)
        real(dp), allocatable :: expected(:, :)
        real(dp) :: host_seconds, transfer_seconds, resident_seconds
        real(dp) :: host_to_resident, transfer_to_resident, checksum
        integer(int64) :: start_count, stop_count, count_rate
        integer :: column, i, n

        n = n1*n2*n3
        allocate( &
            factors(1)%values(n1, n1), factors(2)%values(n2, n2), &
            factors(3)%values(n3, n3), input(n, n_rhs), output(n, n_rhs))
        call fill_factor(factors(1)%values)
        call fill_factor(factors(2)%values)
        call fill_factor(factors(3)%values)
        do column = 1, n_rhs
            do i = 1, n
                input(i, column) = sin(0.013_dp*real(i, dp)) + &
                    0.1_dp*real(column, dp)
            end do
        end do

        call tensor_operator%initialize(factors, status)
        if (status%code /= FORTNUM_OK) error stop "tensor setup failed"
        call tensor_operator%enter_data(status, n_rhs)
        if (status%code /= FORTNUM_OK) error stop "device setup failed"

        call tensor_operator%matmat(input, output, status)
        if (status%code /= FORTNUM_OK) error stop "host warm-up failed"
        call system_clock(count_rate=count_rate)
        call system_clock(start_count)
        do i = 1, repetitions
            call tensor_operator%matmat(input, output, status)
            if (status%code /= FORTNUM_OK) error stop "host timing failed"
        end do
        call system_clock(stop_count)
        host_seconds = real(stop_count - start_count, dp)/ &
            real(count_rate, dp)/real(repetitions, dp)

        !$acc data copyin(input) copyout(output)
        call tensor_operator%matmat_device(input, output, status)
        !$acc end data
        if (status%code /= FORTNUM_OK) error stop "device warm-up failed"
        call system_clock(start_count)
        do i = 1, repetitions
            !$acc data copyin(input) copyout(output)
            call tensor_operator%matmat_device(input, output, status)
            !$acc end data
            if (status%code /= FORTNUM_OK) error stop "device transfer timing failed"
        end do
        call system_clock(stop_count)
        transfer_seconds = real(stop_count - start_count, dp)/ &
            real(count_rate, dp)/real(repetitions, dp)

        !$acc data copyin(input) copyout(output)
        call tensor_operator%matmat_device(input, output, status)
        call system_clock(start_count)
        do i = 1, repetitions
            call tensor_operator%matmat_device(input, output, status)
            if (status%code /= FORTNUM_OK) error stop "device resident timing failed"
        end do
        call system_clock(stop_count)
        resident_seconds = real(stop_count - start_count, dp)/ &
            real(count_rate, dp)/real(repetitions, dp)
        !$acc end data

        checksum = sum(output)
        host_to_resident = host_seconds/resident_seconds
        transfer_to_resident = transfer_seconds/resident_seconds

        if (include_dense) then
            allocate(dense(n, n), expected(n, n_rhs))
            call assemble_dense(factors, n1, n2, n3, dense)
            expected = matmul(dense, input)
            if (maxval(abs(output - expected)) > 2.0e-12_dp) then
                error stop "device benchmark dense oracle failed"
            end if
        end if
        call tensor_operator%exit_data(status)
        if (status%code /= FORTNUM_OK) error stop "device teardown failed"

        write (*, '(a,1x,4(i0,1x),6(es14.6,1x))') "tensor_device", n, &
            n_rhs, repetitions, n1, host_seconds, transfer_seconds, &
            resident_seconds, host_to_resident, transfer_to_resident, checksum
    end subroutine run_case

    subroutine fill_factor(factor)
        real(dp), intent(out) :: factor(:, :)
        integer :: i

        factor = 0.0_dp
        do i = 1, size(factor, 1)
            factor(i, i) = 2.0_dp
        end do
        do i = 2, size(factor, 1)
            factor(i, i - 1) = -0.25_dp
        end do
        do i = 1, size(factor, 1) - 1
            factor(i, i + 1) = -0.25_dp
        end do
    end subroutine fill_factor

    subroutine assemble_dense(factors, n1, n2, n3, dense)
        type(tensor_factor_t), intent(in) :: factors(:)
        integer, intent(in) :: n1, n2, n3
        real(dp), intent(out) :: dense(:, :)
        integer :: i1, i2, i3, j1, j2, j3, row, column

        dense = 0.0_dp
        do i3 = 1, n3
            do i2 = 1, n2
                do i1 = 1, n1
                    row = i1 + n1*(i2 - 1) + n1*n2*(i3 - 1)
                    do j3 = 1, n3
                        do j2 = 1, n2
                            do j1 = 1, n1
                                column = j1 + n1*(j2 - 1) + &
                                    n1*n2*(j3 - 1)
                                dense(row, column) = &
                                    factors(3)%values(i3, j3)* &
                                    factors(2)%values(i2, j2)* &
                                    factors(1)%values(i1, j1)
                            end do
                        end do
                    end do
                end do
            end do
        end do
    end subroutine assemble_dense

end program fortnum_bench_tensor_product_device
