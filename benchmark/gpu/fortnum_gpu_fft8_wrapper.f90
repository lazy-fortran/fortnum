module fortnum_gpu_fft8_wrapper
    use fortnum_kinds, only: dp
    implicit none
    private

    public :: fft8_jvp_batch, fft8_vjp_batch

contains

    subroutine fft8_jvp_batch(batch_size, tangents, products)
        integer, intent(in) :: batch_size
        complex(dp), intent(in) :: tangents(batch_size, 8)
        complex(dp), intent(out) :: products(batch_size, 8)
        integer :: i

        !$acc parallel loop present(tangents, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: tangents) map(from: products)
        do i = 1, batch_size
            call fft8_transform( &
                tangents(i, 1), tangents(i, 2), tangents(i, 3), &
                tangents(i, 4), tangents(i, 5), tangents(i, 6), &
                tangents(i, 7), tangents(i, 8), -1, &
                products(i, 1), products(i, 2), products(i, 3), &
                products(i, 4), products(i, 5), products(i, 6), &
                products(i, 7), products(i, 8))
        end do
    end subroutine fft8_jvp_batch

    subroutine fft8_vjp_batch(batch_size, cotangents, products)
        integer, intent(in) :: batch_size
        complex(dp), intent(in) :: cotangents(batch_size, 8)
        complex(dp), intent(out) :: products(batch_size, 8)
        integer :: i

        !$acc parallel loop present(cotangents, products)
        !$omp target teams distribute parallel do &
        !$omp& map(to: cotangents) map(from: products)
        do i = 1, batch_size
            call fft8_transform( &
                cotangents(i, 1), cotangents(i, 2), cotangents(i, 3), &
                cotangents(i, 4), cotangents(i, 5), cotangents(i, 6), &
                cotangents(i, 7), cotangents(i, 8), 1, &
                products(i, 1), products(i, 2), products(i, 3), &
                products(i, 4), products(i, 5), products(i, 6), &
                products(i, 7), products(i, 8))
        end do
    end subroutine fft8_vjp_batch

    pure subroutine fft8_transform( &
            z1, z2, z3, z4, z5, z6, z7, z8, sign, &
            w1, w2, w3, w4, w5, w6, w7, w8)
        !$omp declare target
        !$acc routine seq
        complex(dp), intent(in) :: z1, z2, z3, z4, z5, z6, z7, z8
        integer, intent(in) :: sign
        complex(dp), intent(out) :: w1, w2, w3, w4, w5, w6, w7, w8
        real(dp), parameter :: root_half = 0.70710678118654752440084436210485_dp
        complex(dp) :: values(8), even, odd, twiddle
        integer :: block, half, index, offset, stage_size

        values(1) = z1
        values(2) = z5
        values(3) = z3
        values(4) = z7
        values(5) = z2
        values(6) = z6
        values(7) = z4
        values(8) = z8
        stage_size = 2
        do while (stage_size <= 8)
            half = stage_size/2
            do block = 1, 8, stage_size
                do offset = 0, half - 1
                    index = block + offset
                    select case (stage_size)
                    case (2)
                        twiddle = cmplx(1.0_dp, 0.0_dp, dp)
                    case (4)
                        if (offset == 0) then
                            twiddle = cmplx(1.0_dp, 0.0_dp, dp)
                        else
                            twiddle = cmplx(0.0_dp, real(sign, dp), dp)
                        end if
                    case default
                        select case (offset)
                        case (0)
                            twiddle = cmplx(1.0_dp, 0.0_dp, dp)
                        case (1)
                            twiddle = cmplx( &
                                root_half, real(sign, dp)*root_half, dp)
                        case (2)
                            twiddle = cmplx(0.0_dp, real(sign, dp), dp)
                        case default
                            twiddle = cmplx( &
                                -root_half, real(sign, dp)*root_half, dp)
                        end select
                    end select
                    even = values(index)
                    odd = twiddle*values(index + half)
                    values(index) = even + odd
                    values(index + half) = even - odd
                end do
            end do
            stage_size = 2*stage_size
        end do
        w1 = values(1)
        w2 = values(2)
        w3 = values(3)
        w4 = values(4)
        w5 = values(5)
        w6 = values(6)
        w7 = values(7)
        w8 = values(8)
    end subroutine fft8_transform

end module fortnum_gpu_fft8_wrapper
