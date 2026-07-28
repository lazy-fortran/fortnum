module fft_jvp_tournament_kernel
    use, intrinsic :: iso_c_binding, only: c_double
    use fortnum_fft8_kernel, only: fft_c2c8
    implicit none
    private

    integer, parameter, public :: fft_size = 8
    integer, parameter, public :: real_size = 2*fft_size
    public :: fft8_primal, analytical_jvp, analytical_vjp

contains

    subroutine fft8_primal(x, y) bind(c, name="fortnum_fft8_primal")
        real(c_double), intent(in) :: x(real_size)
        real(c_double), intent(out) :: y(real_size)
        complex(c_double) :: z(fft_size)

        call pack_complex(x, z)
        call fft_c2c8(z(1), z(2), z(3), z(4), z(5), z(6), z(7), z(8), -1)
        call unpack_complex(z, y)
    end subroutine fft8_primal

    subroutine analytical_jvp(direction, product)
        real(c_double), intent(in) :: direction(real_size)
        real(c_double), intent(out) :: product(real_size)
        complex(c_double) :: tangent(fft_size)

        call pack_complex(direction, tangent)
        call fft_c2c8(tangent(1), tangent(2), tangent(3), tangent(4), &
            tangent(5), tangent(6), tangent(7), tangent(8), -1)
        call unpack_complex(tangent, product)
    end subroutine analytical_jvp

    subroutine analytical_vjp(cotangent, product)
        real(c_double), intent(in) :: cotangent(real_size)
        real(c_double), intent(out) :: product(real_size)
        complex(c_double) :: adjoint(fft_size)

        call pack_complex(cotangent, adjoint)
        call fft_c2c8(adjoint(1), adjoint(2), adjoint(3), adjoint(4), &
            adjoint(5), adjoint(6), adjoint(7), adjoint(8), 1)
        call unpack_complex(adjoint, product)
    end subroutine analytical_vjp

    pure subroutine pack_complex(input, output)
        real(c_double), intent(in) :: input(real_size)
        complex(c_double), intent(out) :: output(fft_size)
        integer :: i

        do i = 1, fft_size
            output(i) = cmplx(input(2*i - 1), input(2*i), c_double)
        end do
    end subroutine pack_complex

    pure subroutine unpack_complex(input, output)
        complex(c_double), intent(in) :: input(fft_size)
        real(c_double), intent(out) :: output(real_size)
        integer :: i

        do i = 1, fft_size
            output(2*i - 1) = real(input(i), c_double)
            output(2*i) = aimag(input(i))
        end do
    end subroutine unpack_complex

end module fft_jvp_tournament_kernel

program enzyme_fft_jvp_tournament
    use, intrinsic :: iso_c_binding, only: c_double
    use fft_jvp_tournament_kernel, only: fft8_primal, analytical_jvp, &
        analytical_vjp, fft_size, real_size
    use fortnum_generated_enzyme_fft8, only: fortnum_enzyme_fft8_jvp
    use fortnum_enzyme_fixture_support, only: collect_fixture_samples, &
        fixture_peak_rss_bytes, fixture_sample_count, fixture_timer_t, &
        median_mad, read_fixture_environment, read_fixture_integer, &
        write_fixture_result
    implicit none

    real(c_double) :: samples(fixture_sample_count), sink
    character(32) :: action, candidate, product_kind
    integer :: directions, iterations
    logical :: valid

    call read_fixture_environment("FORTNUM_FFT_ACTION", "validate", action)
    if (trim(action) == "validate") then
        call validate_candidates()
        write (*, "(a)") "PASS"
        stop
    end if
    call read_fixture_environment("FORTNUM_FFT_CANDIDATE", "analytical", candidate)
    call read_fixture_environment("FORTNUM_FFT_PRODUCT", "jvp", product_kind)
    call read_fixture_integer("FORTNUM_FFT_DIRECTIONS", 16, directions, valid)
    if (.not. valid) error stop "invalid direction count"
    call read_fixture_integer("FORTNUM_FFT_ITERATIONS", 10000, iterations, valid)
    if (.not. valid) error stop "invalid iteration count"
    if (trim(candidate) /= "analytical" .and. trim(candidate) /= "autodiff") &
        error stop "candidate must be analytical or autodiff"
    if (trim(product_kind) /= "jvp" .and. trim(product_kind) /= "vjp") &
        error stop "product must be jvp or vjp"
    if (directions < 1 .or. directions > 16) &
        error stop "directions must be 1..16"
    if (iterations < 1) error stop "iterations must be positive"

    select case (trim(action))
    case ("benchmark", "--benchmark")
        call collect_fixture_samples(measure_candidate, samples)
        call report_samples()
    case ("peak-rss", "--peak-rss")
        sink = measure_candidate()
        write (*, "(i0)") fixture_peak_rss_bytes()
    case default
        error stop "action must be validate, benchmark, or peak-rss"
    end select
    if (sink /= sink) error stop "FFT benchmark produced NaN"

contains

    subroutine validate_candidates()
        real(c_double), parameter :: tolerance = 4.0e-13_c_double
        real(c_double) :: x(real_size), direction(real_size)
        real(c_double) :: primal(real_size), enzyme_primal(real_size)
        real(c_double) :: expected(real_size), analytical(real_size)
        real(c_double) :: autodiff(real_size), scale, lhs, rhs
        integer :: i

        do i = 1, real_size
            x(i) = sin(0.31_c_double*real(i, c_double))
            direction(i) = cos(0.47_c_double*real(i, c_double))
        end do
        call fft8_primal(x, primal)
        call direct_dft(x, expected)
        scale = max(1.0_c_double, maxval(abs(expected)))
        if (maxval(abs(primal - expected)) > tolerance*scale) &
            error stop "production FFT disagrees with direct DFT"

        call analytical_jvp(direction, analytical)
        call direct_dft(direction, expected)
        if (maxval(abs(analytical - expected)) > tolerance*scale) &
            error stop "analytical FFT JVP disagrees with direct DFT"

        call fortnum_enzyme_fft8_jvp(x, direction, enzyme_primal, autodiff)
        if (maxval(abs(enzyme_primal - primal)) > tolerance*scale) &
            error stop "Enzyme FFT primal disagrees with production"
        if (maxval(abs(autodiff - expected)) > tolerance*scale) &
            error stop "autodiff FFT JVP disagrees with direct DFT"

        call analytical_vjp(direction, analytical)
        call direct_dft_adjoint(direction, expected)
        scale = max(1.0_c_double, maxval(abs(expected)))
        if (maxval(abs(analytical - expected)) > tolerance*scale) &
            error stop "analytical FFT VJP disagrees with direct adjoint DFT"

        lhs = dot_product(direction, primal)
        rhs = dot_product(analytical, x)
        if (abs(lhs - rhs) > tolerance*max(1.0_c_double, abs(lhs), abs(rhs))) &
            error stop "FFT adjoint dot identity failed"

        call autodiff_vjp(x, direction, autodiff)
        if (maxval(abs(autodiff - expected)) > tolerance*scale) &
            error stop "autodiff FFT VJP disagrees with direct adjoint DFT"
    end subroutine validate_candidates

    function measure_candidate() result(elapsed_ns)
        type(fixture_timer_t) :: timer
        real(c_double) :: x(real_size), direction(real_size)
        real(c_double) :: primal(real_size), product(real_size)
        real(c_double) :: elapsed_ns, local_sink
        integer :: direction_index, i, iteration

        local_sink = 0.0_c_double
        call timer%start()
        do iteration = 1, iterations
            do i = 1, real_size
                x(i) = sin(0.013_c_double*real(i + iteration, c_double))
            end do
            do direction_index = 1, directions
                do i = 1, real_size
                    direction(i) = cos(0.017_c_double* &
                        real(i + 3*direction_index, c_double))
                end do
                if (trim(product_kind) == "jvp") then
                    if (trim(candidate) == "analytical") then
                        call analytical_jvp(direction, product)
                    else
                        call fortnum_enzyme_fft8_jvp(x, direction, primal, product)
                    end if
                else
                    if (trim(candidate) == "analytical") then
                        call analytical_vjp(direction, product)
                    else
                        call autodiff_vjp(x, direction, product)
                    end if
                end if
                local_sink = local_sink + product( &
                    1 + mod(direction_index - 1, real_size))
            end do
        end do
        elapsed_ns = timer%elapsed_ns()/real(iterations, c_double)
        sink = local_sink
    end function measure_candidate

    subroutine report_samples()
        real(c_double) :: median, mad
        character(96) :: name

        call median_mad(samples, median, mad)
        write (name, "('fft8_',a,'_',a,'_d',i0)") &
            trim(candidate), trim(product_kind), directions
        call write_fixture_result(trim(name), iterations, median, mad)
    end subroutine report_samples

    pure subroutine direct_dft(input, output)
        real(c_double), intent(in) :: input(real_size)
        real(c_double), intent(out) :: output(real_size)
        real(c_double), parameter :: pi = &
            3.141592653589793238462643383279502884_c_double
        complex(c_double) :: value, z(fft_size)
        real(c_double) :: angle
        integer :: frequency, i, sample

        do i = 1, fft_size
            z(i) = cmplx(input(2*i - 1), input(2*i), c_double)
        end do
        do frequency = 0, fft_size - 1
            value = cmplx(0.0_c_double, 0.0_c_double, c_double)
            do sample = 0, fft_size - 1
                angle = -2.0_c_double*pi* &
                    real(frequency*sample, c_double)/real(fft_size, c_double)
                value = value + z(sample + 1)* &
                    cmplx(cos(angle), sin(angle), c_double)
            end do
            output(2*frequency + 1) = real(value, c_double)
            output(2*frequency + 2) = aimag(value)
        end do
    end subroutine direct_dft

    subroutine autodiff_vjp(x, cotangent, product)
        real(c_double), intent(in) :: x(real_size), cotangent(real_size)
        real(c_double), intent(out) :: product(real_size)
        real(c_double) :: basis(real_size), primal(real_size)
        real(c_double) :: column(real_size)
        integer :: input_index

        do input_index = 1, real_size
            basis = 0.0_c_double
            basis(input_index) = 1.0_c_double
            call fortnum_enzyme_fft8_jvp(x, basis, primal, column)
            product(input_index) = dot_product(cotangent, column)
        end do
    end subroutine autodiff_vjp

    pure subroutine direct_dft_adjoint(input, output)
        real(c_double), intent(in) :: input(real_size)
        real(c_double), intent(out) :: output(real_size)
        real(c_double), parameter :: pi = &
            3.141592653589793238462643383279502884_c_double
        complex(c_double) :: value, z(fft_size)
        real(c_double) :: angle
        integer :: frequency, i, sample

        do i = 1, fft_size
            z(i) = cmplx(input(2*i - 1), input(2*i), c_double)
        end do
        do sample = 0, fft_size - 1
            value = cmplx(0.0_c_double, 0.0_c_double, c_double)
            do frequency = 0, fft_size - 1
                angle = 2.0_c_double*pi* &
                    real(frequency*sample, c_double)/real(fft_size, c_double)
                value = value + z(frequency + 1)* &
                    cmplx(cos(angle), sin(angle), c_double)
            end do
            output(2*sample + 1) = real(value, c_double)
            output(2*sample + 2) = aimag(value)
        end do
    end subroutine direct_dft_adjoint

end program enzyme_fft_jvp_tournament
