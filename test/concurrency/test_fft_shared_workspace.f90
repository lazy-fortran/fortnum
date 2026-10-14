program test_fft_shared_workspace
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_fft, only: fortnum_fft_plan_t, fft_c2c_plan_init, fft_c2c, &
        fortnum_fft_workspace_t, fft_workspace_init
    implicit none
    type(fortnum_fft_plan_t) :: plan
    type(fortnum_fft_workspace_t) :: scratch(4)
    complex(dp), allocatable :: z(:, :), reference(:, :), expected(:, :)
    integer, parameter :: sizes(5) = [8, 9, 17, 25, 33]
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    real(dp) :: angle, error
    integer :: size_index, n, j, thread, repeat, failures
    character(16) :: command
    call get_command_argument(1, command)
    do size_index = 1, size(sizes)
        n = sizes(size_index)
        allocate (z(n, 4), reference(n, 4), expected(n, 4))
        call fft_c2c_plan_init(plan, n)
        do thread = 1, 4
            call fft_workspace_init(scratch(thread), plan)
            do j = 1, n
                angle = 2*acos(-1.0_dp)*(j - 1)/n
                reference(j, thread) = exp(iu*thread*angle) &
                    + (0.3_dp, -0.1_dp)*exp(-iu*2*angle)
            end do
            expected(:, thread) = 0
            expected(thread + 1, thread) = n
            expected(n - 1, thread) = n*(0.3_dp, -0.1_dp)
        end do
        failures = 0
        !$omp parallel do private(repeat,error) reduction(+:failures) num_threads(4)
        do thread = 1, 4
            do repeat = 1, 64
                z(:, thread) = reference(:, thread)
                if (command == 'legacy') then
                    call fft_c2c(z(:, thread), -1, plan)
                else
                    call fft_c2c(z(:, thread), -1, plan, scratch(thread))
                end if
                error = maxval(abs(z(:, thread) - expected(:, thread)))
                if (error > 2.0e-12_dp*n) failures = failures + 1
                if (command == 'legacy') then
                    call fft_c2c(z(:, thread), 1, plan)
                else
                    call fft_c2c(z(:, thread), 1, plan, scratch(thread))
                end if
                error = maxval(abs(z(:, thread)/n - reference(:, thread)))
                if (error > 2.0e-12_dp) failures = failures + 1
            end do
        end do
        !$omp end parallel do
        if (failures /= 0) error stop 'shared plan/private scratch spectral oracle'
        deallocate (z, reference, expected)
    end do
    print *, 'PASS concurrent mixed-radix/Bluestein spectral and inverse oracles'
end program
