module fortnum_fft8_kernel
    use fortnum_kinds, only: dp
    implicit none
    private

    public :: fft_c2c8

contains

    pure subroutine fft_c2c8(z1, z2, z3, z4, z5, z6, z7, z8, sign)
        !$omp declare target
        !$acc routine seq
        complex(dp), intent(inout) :: z1, z2, z3, z4, z5, z6, z7, z8
        integer, intent(in) :: sign
        real(dp), parameter :: root_half = &
            0.707106781186547524400844362104849039_dp
        complex(dp) :: even, odd, twiddle, values(8)
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
        z1 = values(1)
        z2 = values(2)
        z3 = values(3)
        z4 = values(4)
        z5 = values(5)
        z6 = values(6)
        z7 = values(7)
        z8 = values(8)
    end subroutine fft_c2c8

end module fortnum_fft8_kernel
