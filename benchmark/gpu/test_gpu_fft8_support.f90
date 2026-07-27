module test_gpu_fft8_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_gpu_fft8_wrapper, only: fft8_jvp_batch, fft8_vjp_batch
    implicit none
    private

    public :: run_fft8_test

contains

    subroutine run_fft8_test()
        integer, parameter :: batch_size = 4096
        real(dp), parameter :: tolerance = 2.0e-13_dp
        complex(dp) :: tangents(batch_size, 8), cotangents(batch_size, 8)
        complex(dp) :: jvps(batch_size, 8), vjps(batch_size, 8)
        complex(dp) :: expected(8)
        real(dp) :: lhs, rhs, scale
        integer :: batch_index, entry

        do entry = 1, 8
            do batch_index = 1, batch_size
                tangents(batch_index, entry) = cmplx( &
                    sin(real(3*batch_index + entry, dp)/19.0_dp), &
                    cos(real(batch_index + 5*entry, dp)/23.0_dp), dp)
                cotangents(batch_index, entry) = cmplx( &
                    cos(real(2*batch_index + entry, dp)/17.0_dp), &
                    sin(real(batch_index + 7*entry, dp)/29.0_dp), dp)
            end do
        end do

        !$acc data copyin(tangents, cotangents) create(jvps, vjps)
        !$omp target data map(to: tangents, cotangents) map(alloc: jvps, vjps)
        call fft8_jvp_batch(batch_size, tangents, jvps)
        call fft8_vjp_batch(batch_size, cotangents, vjps)
        !$omp target update from(jvps, vjps)
        !$omp end target data
        !$acc update self(jvps, vjps)
        !$acc end data

        do batch_index = 1, batch_size
            call direct_dft(tangents(batch_index, :), -1, expected)
            scale = max(1.0_dp, maxval(abs(expected)))
            if (maxval(abs(expected - jvps(batch_index, :))) > &
                tolerance*scale) then
                error stop "GPU FFT8 JVP disagrees with direct DFT"
            end if
            call direct_dft(cotangents(batch_index, :), 1, expected)
            scale = max(1.0_dp, maxval(abs(expected)))
            if (maxval(abs(expected - vjps(batch_index, :))) > &
                tolerance*scale) then
                error stop "GPU FFT8 VJP disagrees with direct adjoint DFT"
            end if
        end do

        lhs = real(sum(conjg(cotangents)*jvps), dp)
        rhs = real(sum(conjg(vjps)*tangents), dp)
        scale = max(1.0_dp, abs(lhs), abs(rhs))
        if (abs(lhs - rhs) > tolerance*scale) then
            error stop "GPU FFT8 products violate complex adjoint identity"
        end if
    end subroutine run_fft8_test

    pure subroutine direct_dft(input, sign, output)
        complex(dp), intent(in) :: input(8)
        integer, intent(in) :: sign
        complex(dp), intent(out) :: output(8)
        real(dp), parameter :: pi = 3.1415926535897932384626433832795_dp
        real(dp) :: angle
        integer :: frequency, sample

        do frequency = 0, 7
            output(frequency + 1) = cmplx(0.0_dp, 0.0_dp, dp)
            do sample = 0, 7
                angle = real(sign, dp)*2.0_dp*pi* &
                    real(frequency*sample, dp)/8.0_dp
                output(frequency + 1) = output(frequency + 1) + &
                    input(sample + 1)*cmplx(cos(angle), sin(angle), dp)
            end do
        end do
    end subroutine direct_dft

end module test_gpu_fft8_support
