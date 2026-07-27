program test_determinant_jvp
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_linalg, only: det2, det3, det2_jvp, det3_jvp
    implicit none

    integer, parameter :: MAX_DIRECTIONS = 64
    character(32) :: action, candidate
    integer :: matrix_size, direction_count
    integer(int64) :: iteration_count

    call get_environment_variable("FORTNUM_DET_ACTION", action)
    call get_environment_variable("FORTNUM_DET_CANDIDATE", candidate)
    call read_integer_env("FORTNUM_DET_SIZE", matrix_size, 3)
    call read_integer_env("FORTNUM_DET_DIRECTIONS", direction_count, 1)
    call read_int64_env("FORTNUM_DET_ITERATIONS", iteration_count, 100000_int64)

    if (trim(action) == "--benchmark") then
        call benchmark_candidate(trim(candidate), matrix_size, direction_count, &
            iteration_count)
    else
        call validate_products()
    end if

contains

    subroutine validate_products()
        real(dp), parameter :: h = 1.0e-5_dp
        real(dp) :: a2(2, 2), v2(2, 2), a3(3, 3), v3(3, 3)
        real(dp) :: jv, fd_h, fd_half

        a2 = reshape([2.0_dp, -0.5_dp, 1.25_dp, 3.0_dp], [2, 2])
        v2 = reshape([0.2_dp, 0.3_dp, -0.4_dp, 0.1_dp], [2, 2])
        call det2_jvp(a2, v2, jv)
        fd_h = central_difference2(a2, v2, h)
        fd_half = central_difference2(a2, v2, 0.5_dp*h)
        if (abs(jv - fd_half) > 2.0e-10_dp) then
            error stop "det2 JVP finite-difference mismatch"
        end if
        if (abs(fd_h - fd_half) > 2.0e-10_dp) then
            error stop "det2 finite-difference refinement mismatch"
        end if

        a3 = reshape([2.0_dp, -0.5_dp, 0.3_dp, 1.25_dp, 3.0_dp, -0.8_dp, &
            0.4_dp, 0.7_dp, 1.8_dp], [3, 3])
        v3 = reshape([0.2_dp, 0.3_dp, -0.4_dp, 0.1_dp, -0.2_dp, 0.5_dp, &
            -0.3_dp, 0.6_dp, 0.25_dp], [3, 3])
        call det3_jvp(a3, v3, jv)
        fd_h = central_difference3(a3, v3, h)
        fd_half = central_difference3(a3, v3, 0.5_dp*h)
        if (abs(jv - fd_half) > 2.0e-9_dp) then
            error stop "det3 JVP finite-difference mismatch"
        end if
        if (abs(fd_h - fd_half) > 2.0e-9_dp) then
            error stop "det3 finite-difference refinement mismatch"
        end if
    end subroutine validate_products

    subroutine benchmark_candidate(name, n, ndirection, iterations)
        character(*), intent(in) :: name
        integer, intent(in) :: n, ndirection
        integer(int64), intent(in) :: iterations
        real(dp) :: a2(2, 2), directions2(2, 2, MAX_DIRECTIONS)
        real(dp) :: a3(3, 3), directions3(3, 3, MAX_DIRECTIONS)
        real(dp) :: value, jv, sink, elapsed_ns
        integer(int64) :: start, finish, clock_rate, iteration
        integer :: direction

        if (n /= 2 .and. n /= 3) error stop "matrix size must be 2 or 3"
        if (ndirection < 1 .or. ndirection > MAX_DIRECTIONS) then
            error stop "direction count must be between 1 and 64"
        end if
        if (name /= "analytical" .and. name /= "diagnostic" .and. &
            name /= "primal") then
            error stop "candidate must be analytical, diagnostic, or primal"
        end if

        call initialize_inputs(a2, directions2, a3, directions3)
        sink = 0.0_dp
        call system_clock(start, clock_rate)
        do iteration = 1, iterations
            if (n == 2) then
                a2(1, 1) = 2.0_dp + 1.0e-12_dp*real(iand(iteration, 1023_int64), dp)
                value = det2(a2)
                if (name == "analytical") then
                    do direction = 1, ndirection
                        call det2_jvp(a2, directions2(:, :, direction), jv)
                        sink = sink + jv
                    end do
                else if (name == "diagnostic") then
                    do direction = 1, ndirection
                        jv = central_difference2(a2, &
                            directions2(:, :, direction), 1.0e-5_dp)
                        sink = sink + jv
                    end do
                end if
            else
                a3(1, 1) = 2.0_dp + 1.0e-12_dp*real(iand(iteration, 1023_int64), dp)
                value = det3(a3)
                if (name == "analytical") then
                    do direction = 1, ndirection
                        call det3_jvp(a3, directions3(:, :, direction), jv)
                        sink = sink + jv
                    end do
                else if (name == "diagnostic") then
                    do direction = 1, ndirection
                        jv = central_difference3(a3, &
                            directions3(:, :, direction), 1.0e-5_dp)
                        sink = sink + jv
                    end do
                end if
            end if
            sink = sink + value
        end do
        call system_clock(finish)
        elapsed_ns = real(finish - start, dp)*1.0e9_dp / &
            (real(clock_rate, dp)*real(iterations, dp))
        if (sink /= sink) error stop "benchmark produced NaN"
        write (*, "(a,a,a,i0,a,i0,a,i0,a,f0.6,a,es12.4)") &
            "candidate=", name, " size=", n, " directions=", ndirection, &
            " iterations=", iterations, " ns_per_workload=", elapsed_ns, &
            " sink=", sink
    end subroutine benchmark_candidate

    subroutine initialize_inputs(a2, directions2, a3, directions3)
        real(dp), intent(out) :: a2(2, 2), directions2(2, 2, MAX_DIRECTIONS)
        real(dp), intent(out) :: a3(3, 3), directions3(3, 3, MAX_DIRECTIONS)
        integer :: direction, i, j

        a2 = reshape([2.0_dp, -0.5_dp, 1.25_dp, 3.0_dp], [2, 2])
        a3 = reshape([2.0_dp, -0.5_dp, 0.3_dp, 1.25_dp, 3.0_dp, -0.8_dp, &
            0.4_dp, 0.7_dp, 1.8_dp], [3, 3])
        do direction = 1, MAX_DIRECTIONS
            do j = 1, 2
                do i = 1, 2
                    directions2(i, j, direction) = &
                        0.01_dp*real(i + 2*j + direction, dp)
                end do
            end do
            do j = 1, 3
                do i = 1, 3
                    directions3(i, j, direction) = &
                        0.01_dp*real(i + 3*j + direction, dp)
                end do
            end do
        end do
    end subroutine initialize_inputs

    pure function central_difference2(a, v, h) result(jv)
        real(dp), intent(in) :: a(2, 2), v(2, 2), h
        real(dp) :: jv, plus(2, 2), minus(2, 2)

        plus = a + h*v
        minus = a - h*v
        jv = (det2(plus) - det2(minus))/(2.0_dp*h)
    end function central_difference2

    pure function central_difference3(a, v, h) result(jv)
        real(dp), intent(in) :: a(3, 3), v(3, 3), h
        real(dp) :: jv, plus(3, 3), minus(3, 3)

        plus = a + h*v
        minus = a - h*v
        jv = (det3(plus) - det3(minus))/(2.0_dp*h)
    end function central_difference3

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

end program test_determinant_jvp
