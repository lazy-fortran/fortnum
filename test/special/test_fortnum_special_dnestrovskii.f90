program test_dnestrovskii
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special_dnestrovskii, only: dnestrovskii
    implicit none
    complex(dp) :: f, df, plus, minus, unused
    real(dp) :: row(9), err, derr, fd, h, maximum, dmaximum
    character(1024) :: path
    integer :: unit, ios, status, q, j, count
    call get_command_argument(1, path)
    open (newunit=unit, file=trim(path), status='old')
    maximum = 0.0_dp
    dmaximum = 0.0_dp
    count = 0
    do
        read (unit, *, iostat=ios) row
        if (ios /= 0) exit
        do j = 1, 2
            q = 2*j + 3
            call dnestrovskii(q, row(1), f, df, status)
            if (status /= 0) error stop 'valid input rejected'
            err = abs(f - cmplx(row(4*j - 2), row(4*j - 1), dp))
            derr = abs(df - cmplx(row(4*j), row(4*j + 1), dp))
            maximum = max(maximum, err)
            dmaximum = max(dmaximum, derr)
            if (err > 3.0e-14_dp) error stop '90-digit defining-function value oracle'
            if (derr > 3.0e-13_dp) error stop '90-digit derivative oracle'
            if (abs(row(1)) > 0.01_dp) then
                h = 1.0e-5_dp*max(1.0_dp, abs(row(1)))
                call dnestrovskii(q, row(1) + h, plus, unused, status)
                call dnestrovskii(q, row(1) - h, minus, unused, status)
                fd = abs((plus - minus)/(2*h) - df)
                if (fd > 3.0e-9_dp) error stop 'finite difference derivative'
            end if
            count = count + 1
        end do
    end do
    close (unit)
    call dnestrovskii(3, -1.0_dp, f, df, status)
    if (status == 0) error stop 'unsupported order accepted'
    print '(A,I0,A,ES12.4,A,ES12.4)', 'oracle samples=', count, &
        ' max absolute value error=', maximum, ' derivative error=', dmaximum
end program
