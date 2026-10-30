program bench_gram_factor
    use fortnum_kinds, only: dp
    use fortnum_gram_factor, only: gram_factor_workspace_t
    use fortnum_status, only: fortnum_status_t, status_ok
    implicit none
    type(gram_factor_workspace_t) :: cached
    type(fortnum_status_t) :: status
    real(dp) :: w(7, 16), b(7, 7), started, finished
    integer :: i, j, sample, iteration
    integer, parameter :: repetitions = 10000
    do j = 1, 16
        do i = 1, 7
            w(i, j) = sin(real(3*i+j, dp)) + real(mod(i+j, 3), dp)
        end do
    end do
    call cached%factor(w, b, status)
    if (.not. status_ok(status)) stop 1
    print *, 'n=7 m=16; CPU seconds per factor, 10000 calls/sample'
    do sample = 1, 5
        call cpu_time(started)
        do iteration = 1, repetitions
            call cached%factor(w, b, status)
            if (.not. status_ok(status)) stop 1
        end do
        call cpu_time(finished)
        print *, 'cached', sample, (finished-started)/repetitions
        call cpu_time(started)
        do iteration = 1, repetitions
            block
                type(gram_factor_workspace_t) :: fresh
                call fresh%factor(w, b, status)
                if (.not. status_ok(status)) stop 1
            end block
        end do
        call cpu_time(finished)
        print *, 'fresh ', sample, (finished-started)/repetitions
    end do
end program
