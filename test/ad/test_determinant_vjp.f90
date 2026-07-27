program test_determinant_vjp
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_linalg, only: det2, det3, det2_jvp, det3_jvp, det2_vjp, det3_vjp
    implicit none

    character(32) :: action, candidate
    integer :: matrix_size
    integer(int64) :: iteration_count

    call get_environment_variable("FORTNUM_DET_VJP_ACTION", action)
    call get_environment_variable("FORTNUM_DET_VJP_CANDIDATE", candidate)
    call read_integer_env("FORTNUM_DET_VJP_SIZE", matrix_size, 3)
    call read_int64_env("FORTNUM_DET_VJP_ITERATIONS", iteration_count, 100000_int64)

    if (trim(action) == "--benchmark") then
        call benchmark_candidate(trim(candidate), matrix_size, iteration_count)
    else
        call validate_products()
    end if

contains

    subroutine validate_products()
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: a2(2, 2), v2(2, 2), bar2(2, 2), fd2(2, 2)
        real(dp) :: a3(3, 3), v3(3, 3), bar3(3, 3), fd3(3, 3)
        real(dp) :: u, jv

        u = -0.7_dp
        a2 = reshape([2.0_dp, -0.5_dp, 1.25_dp, 3.0_dp], [2, 2])
        v2 = reshape([0.2_dp, 0.3_dp, -0.4_dp, 0.1_dp], [2, 2])
        call det2_vjp(a2, u, bar2)
        call central_gradient2(a2, u, h, fd2)
        if (maxval(abs(bar2 - fd2)) > 2.0e-10_dp) then
            error stop "det2 VJP finite-difference mismatch"
        end if
        call det2_jvp(a2, v2, jv)
        if (abs(sum(bar2*v2) - u*jv) > 2.0e-14_dp) then
            error stop "det2 adjoint identity mismatch"
        end if

        a3 = reshape([2.0_dp, -0.5_dp, 0.3_dp, 1.25_dp, 3.0_dp, -0.8_dp, &
            0.4_dp, 0.7_dp, 1.8_dp], [3, 3])
        v3 = reshape([0.2_dp, 0.3_dp, -0.4_dp, 0.1_dp, -0.2_dp, 0.5_dp, &
            -0.3_dp, 0.6_dp, 0.25_dp], [3, 3])
        call det3_vjp(a3, u, bar3)
        call central_gradient3(a3, u, h, fd3)
        if (maxval(abs(bar3 - fd3)) > 2.0e-9_dp) then
            error stop "det3 VJP finite-difference mismatch"
        end if
        call det3_jvp(a3, v3, jv)
        if (abs(sum(bar3*v3) - u*jv) > 2.0e-14_dp) then
            error stop "det3 adjoint identity mismatch"
        end if
    end subroutine validate_products

    subroutine benchmark_candidate(name, n, iterations)
        character(*), intent(in) :: name
        integer, intent(in) :: n
        integer(int64), intent(in) :: iterations
        real(dp) :: a2(2, 2), bar2(2, 2)
        real(dp) :: a3(3, 3), bar3(3, 3)
        real(dp) :: value, u, sink, elapsed_ns
        integer(int64) :: start, finish, clock_rate, iteration

        if (n /= 2 .and. n /= 3) error stop "matrix size must be 2 or 3"
        if (name /= "analytical" .and. name /= "diagnostic" .and. &
            name /= "primal") then
            error stop "candidate must be analytical, diagnostic, or primal"
        end if

        a2 = reshape([2.0_dp, -0.5_dp, 1.25_dp, 3.0_dp], [2, 2])
        a3 = reshape([2.0_dp, -0.5_dp, 0.3_dp, 1.25_dp, 3.0_dp, -0.8_dp, &
            0.4_dp, 0.7_dp, 1.8_dp], [3, 3])
        sink = 0.0_dp
        call system_clock(start, clock_rate)
        do iteration = 1, iterations
            u = 0.5_dp + 1.0e-6_dp*real(iand(iteration, 1023_int64), dp)
            if (n == 2) then
                a2(1, 1) = 2.0_dp + 1.0e-12_dp*real(iand(iteration, 1023_int64), dp)
                value = det2(a2)
                if (name == "analytical") then
                    call det2_vjp(a2, u, bar2)
                    sink = sink + sum(bar2)
                else if (name == "diagnostic") then
                    call central_gradient2(a2, u, 1.0e-5_dp, bar2)
                    sink = sink + sum(bar2)
                end if
            else
                a3(1, 1) = 2.0_dp + 1.0e-12_dp*real(iand(iteration, 1023_int64), dp)
                value = det3(a3)
                if (name == "analytical") then
                    call det3_vjp(a3, u, bar3)
                    sink = sink + sum(bar3)
                else if (name == "diagnostic") then
                    call central_gradient3(a3, u, 1.0e-5_dp, bar3)
                    sink = sink + sum(bar3)
                end if
            end if
            sink = sink + value
        end do
        call system_clock(finish)
        elapsed_ns = real(finish - start, dp)*1.0e9_dp / &
            (real(clock_rate, dp)*real(iterations, dp))
        if (sink /= sink) error stop "benchmark produced NaN"
        write (*, "(a,a,a,i0,a,i0,a,f0.6,a,es12.4)") &
            "candidate=", name, " size=", n, " iterations=", iterations, &
            " ns_per_workload=", elapsed_ns, " sink=", sink
    end subroutine benchmark_candidate

    pure subroutine central_gradient2(a, u, h, bar)
        real(dp), intent(in) :: a(2, 2), u, h
        real(dp), intent(out) :: bar(2, 2)
        real(dp) :: plus(2, 2), minus(2, 2)
        integer :: i, j

        do j = 1, 2
            do i = 1, 2
                plus = a
                minus = a
                plus(i, j) = plus(i, j) + h
                minus(i, j) = minus(i, j) - h
                bar(i, j) = u*(det2(plus) - det2(minus))/(2.0_dp*h)
            end do
        end do
    end subroutine central_gradient2

    pure subroutine central_gradient3(a, u, h, bar)
        real(dp), intent(in) :: a(3, 3), u, h
        real(dp), intent(out) :: bar(3, 3)
        real(dp) :: plus(3, 3), minus(3, 3)
        integer :: i, j

        do j = 1, 3
            do i = 1, 3
                plus = a
                minus = a
                plus(i, j) = plus(i, j) + h
                minus(i, j) = minus(i, j) - h
                bar(i, j) = u*(det3(plus) - det3(minus))/(2.0_dp*h)
            end do
        end do
    end subroutine central_gradient3

    subroutine read_integer_env(name, value, default_value)
        character(*), intent(in) :: name
        integer, intent(out) :: value
        integer, intent(in) :: default_value
        character(32) :: text
        integer :: status, ios

        value = default_value
        call get_environment_variable(name, text, status=status)
        if (status /= 0 .or. len_trim(text) == 0) return
        read (text, *, iostat=ios) value
        if (ios /= 0) error stop "invalid integer environment value"
    end subroutine read_integer_env

    subroutine read_int64_env(name, value, default_value)
        character(*), intent(in) :: name
        integer(int64), intent(out) :: value
        integer(int64), intent(in) :: default_value
        character(32) :: text
        integer :: status, ios

        value = default_value
        call get_environment_variable(name, text, status=status)
        if (status /= 0 .or. len_trim(text) == 0) return
        read (text, *, iostat=ios) value
        if (ios /= 0) error stop "invalid int64 environment value"
    end subroutine read_int64_env

end program test_determinant_vjp
