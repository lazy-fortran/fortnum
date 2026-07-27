program test_inverse_vjp
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_linalg, only: inv2, inv3, inv2_jvp, inv3_jvp, inv2_vjp, inv3_vjp, &
        LINALG_OK, LINALG_SINGULAR
    implicit none

    character(32) :: action, candidate
    integer :: matrix_size
    integer(int64) :: iteration_count

    call get_environment_variable("FORTNUM_INV_VJP_ACTION", action)
    call get_environment_variable("FORTNUM_INV_VJP_CANDIDATE", candidate)
    call read_integer_env("FORTNUM_INV_VJP_SIZE", matrix_size, 3)
    call read_int64_env("FORTNUM_INV_VJP_ITERATIONS", iteration_count, 100000_int64)

    if (trim(action) == "--benchmark") then
        call benchmark_candidate(trim(candidate), matrix_size, iteration_count)
    else
        call validate_products()
    end if

contains

    subroutine validate_products()
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: a2(2, 2), u2(2, 2), va2(2, 2), ainv2(2, 2)
        real(dp) :: abar2(2, 2), fd2_h(2, 2), fd2_half(2, 2), vainv2(2, 2)
        real(dp) :: a3(3, 3), u3(3, 3), va3(3, 3), ainv3(3, 3)
        real(dp) :: abar3(3, 3), fd3_h(3, 3), fd3_half(3, 3), vainv3(3, 3)
        integer :: info

        a2 = reshape([2.0_dp, -0.5_dp, 1.25_dp, 3.0_dp], [2, 2])
        u2 = reshape([0.4_dp, -0.2_dp, 0.3_dp, 0.7_dp], [2, 2])
        va2 = reshape([0.2_dp, 0.3_dp, -0.4_dp, 0.1_dp], [2, 2])
        call inv2_vjp(a2, u2, ainv2, abar2, info)
        if (info /= LINALG_OK) error stop "inv2 VJP rejected regular matrix"
        call component_difference2(a2, u2, h, fd2_h)
        call component_difference2(a2, u2, 0.5_dp*h, fd2_half)
        if (maxval(abs(abar2 - fd2_half)) > 2.0e-10_dp) then
            error stop "inv2 VJP finite-difference mismatch"
        end if
        if (maxval(abs(fd2_h - fd2_half)) > 2.0e-10_dp) then
            error stop "inv2 finite-difference refinement mismatch"
        end if
        call inv2_jvp(a2, va2, ainv2, vainv2, info)
        if (abs(sum(u2*vainv2) - sum(abar2*va2)) > 2.0e-14_dp) then
            error stop "inv2 adjoint identity mismatch"
        end if

        a3 = reshape([2.0_dp, -0.5_dp, 0.3_dp, 1.25_dp, 3.0_dp, -0.8_dp, &
            0.4_dp, 0.7_dp, 1.8_dp], [3, 3])
        u3 = reshape([0.4_dp, -0.2_dp, 0.3_dp, 0.7_dp, 0.1_dp, -0.6_dp, &
            0.2_dp, 0.8_dp, -0.5_dp], [3, 3])
        va3 = reshape([0.2_dp, 0.3_dp, -0.4_dp, 0.1_dp, -0.2_dp, 0.5_dp, &
            -0.3_dp, 0.6_dp, 0.25_dp], [3, 3])
        call inv3_vjp(a3, u3, ainv3, abar3, info)
        if (info /= LINALG_OK) error stop "inv3 VJP rejected regular matrix"
        call component_difference3(a3, u3, h, fd3_h)
        call component_difference3(a3, u3, 0.5_dp*h, fd3_half)
        if (maxval(abs(abar3 - fd3_half)) > 2.0e-9_dp) then
            error stop "inv3 VJP finite-difference mismatch"
        end if
        if (maxval(abs(fd3_h - fd3_half)) > 2.0e-9_dp) then
            error stop "inv3 finite-difference refinement mismatch"
        end if
        call inv3_jvp(a3, va3, ainv3, vainv3, info)
        if (abs(sum(u3*vainv3) - sum(abar3*va3)) > 2.0e-14_dp) then
            error stop "inv3 adjoint identity mismatch"
        end if

        a2 = reshape([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], [2, 2])
        call inv2_vjp(a2, u2, ainv2, abar2, info)
        if (info /= LINALG_SINGULAR) error stop "inv2 VJP accepted singular matrix"
        if (maxval(abs(ainv2)) /= 0.0_dp) error stop "singular inv2 value not zero"
        if (maxval(abs(abar2)) /= 0.0_dp) error stop "singular inv2 VJP not zero"
    end subroutine validate_products

    subroutine benchmark_candidate(name, n, iterations)
        character(*), intent(in) :: name
        integer, intent(in) :: n
        integer(int64), intent(in) :: iterations
        real(dp) :: a2(2, 2), u2(2, 2), ainv2(2, 2), abar2(2, 2)
        real(dp) :: a3(3, 3), u3(3, 3), ainv3(3, 3), abar3(3, 3)
        real(dp) :: sink, elapsed_ns
        integer(int64) :: start, finish, clock_rate, iteration
        integer :: info

        if (n /= 2 .and. n /= 3) error stop "matrix size must be 2 or 3"
        if (name /= "analytical" .and. name /= "diagnostic" .and. &
            name /= "primal") then
            error stop "candidate must be analytical, diagnostic, or primal"
        end if

        a2 = reshape([2.0_dp, -0.5_dp, 1.25_dp, 3.0_dp], [2, 2])
        u2 = reshape([0.4_dp, -0.2_dp, 0.3_dp, 0.7_dp], [2, 2])
        a3 = reshape([2.0_dp, -0.5_dp, 0.3_dp, 1.25_dp, 3.0_dp, -0.8_dp, &
            0.4_dp, 0.7_dp, 1.8_dp], [3, 3])
        u3 = reshape([0.4_dp, -0.2_dp, 0.3_dp, 0.7_dp, 0.1_dp, -0.6_dp, &
            0.2_dp, 0.8_dp, -0.5_dp], [3, 3])
        sink = 0.0_dp
        call system_clock(start, clock_rate)
        do iteration = 1, iterations
            if (n == 2) then
                a2(1, 1) = 2.0_dp + 1.0e-12_dp*real(iand(iteration, 1023_int64), dp)
                if (name == "analytical") then
                    call inv2_vjp(a2, u2, ainv2, abar2, info)
                    sink = sink + sum(abar2)
                else
                    call inv2(a2, ainv2, info)
                    if (name == "diagnostic") then
                        call component_difference2(a2, u2, 1.0e-5_dp, abar2)
                        sink = sink + sum(abar2)
                    end if
                end if
                sink = sink + sum(ainv2)
            else
                a3(1, 1) = 2.0_dp + 1.0e-12_dp*real(iand(iteration, 1023_int64), dp)
                if (name == "analytical") then
                    call inv3_vjp(a3, u3, ainv3, abar3, info)
                    sink = sink + sum(abar3)
                else
                    call inv3(a3, ainv3, info)
                    if (name == "diagnostic") then
                        call component_difference3(a3, u3, 1.0e-5_dp, abar3)
                        sink = sink + sum(abar3)
                    end if
                end if
                sink = sink + sum(ainv3)
            end if
            if (info /= LINALG_OK) error stop "benchmark inverse failed"
        end do
        call system_clock(finish)
        elapsed_ns = real(finish - start, dp)*1.0e9_dp / &
            (real(clock_rate, dp)*real(iterations, dp))
        if (sink /= sink) error stop "benchmark produced NaN"
        write (*, "(a,a,a,i0,a,i0,a,f0.6,a,es12.4)") &
            "candidate=", name, " size=", n, " iterations=", iterations, &
            " ns_per_workload=", elapsed_ns, " sink=", sink
    end subroutine benchmark_candidate

    pure subroutine component_difference2(a, u, h, abar)
        real(dp), intent(in) :: a(2, 2), u(2, 2), h
        real(dp), intent(out) :: abar(2, 2)
        real(dp) :: plus(2, 2), minus(2, 2), plus_inv(2, 2), minus_inv(2, 2)
        integer :: i, j, info

        do j = 1, 2
            do i = 1, 2
                plus = a
                minus = a
                plus(i, j) = plus(i, j) + h
                minus(i, j) = minus(i, j) - h
                call inv2(plus, plus_inv, info)
                if (info /= LINALG_OK) error stop "FD plus inverse failed"
                call inv2(minus, minus_inv, info)
                if (info /= LINALG_OK) error stop "FD minus inverse failed"
                abar(i, j) = sum(u*(plus_inv - minus_inv))/(2.0_dp*h)
            end do
        end do
    end subroutine component_difference2

    pure subroutine component_difference3(a, u, h, abar)
        real(dp), intent(in) :: a(3, 3), u(3, 3), h
        real(dp), intent(out) :: abar(3, 3)
        real(dp) :: plus(3, 3), minus(3, 3), plus_inv(3, 3), minus_inv(3, 3)
        integer :: i, j, info

        do j = 1, 3
            do i = 1, 3
                plus = a
                minus = a
                plus(i, j) = plus(i, j) + h
                minus(i, j) = minus(i, j) - h
                call inv3(plus, plus_inv, info)
                if (info /= LINALG_OK) error stop "FD plus inverse failed"
                call inv3(minus, minus_inv, info)
                if (info /= LINALG_OK) error stop "FD minus inverse failed"
                abar(i, j) = sum(u*(plus_inv - minus_inv))/(2.0_dp*h)
            end do
        end do
    end subroutine component_difference3

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

end program test_inverse_vjp
