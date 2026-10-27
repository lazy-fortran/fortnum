module fortnum_special_dnestrovskii
    !! F_q(z) on the upper lip of the real axis, q=5/2 or 7/2.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_dnestro_re_5_1_m, only: fortnum_dnestro_re_5_1
    use fortnum_dnestro_re_5_2_m, only: fortnum_dnestro_re_5_2
    use fortnum_dnestro_re_5_3_m, only: fortnum_dnestro_re_5_3
    use fortnum_dnestro_re_5_4_m, only: fortnum_dnestro_re_5_4
    use fortnum_dnestro_re_5_5_m, only: fortnum_dnestro_re_5_5
    use fortnum_dnestro_re_5_6_m, only: fortnum_dnestro_re_5_6
    use fortnum_dnestro_re_5_7_m, only: fortnum_dnestro_re_5_7
    use fortnum_dnestro_re_5_8_m, only: fortnum_dnestro_re_5_8
    use fortnum_dnestro_tail_5_m, only: fortnum_dnestro_tail_5
    use fortnum_dnestro_im_5_m, only: fortnum_dnestro_im_5
    use fortnum_dnestro_re_7_1_m, only: fortnum_dnestro_re_7_1
    use fortnum_dnestro_re_7_2_m, only: fortnum_dnestro_re_7_2
    use fortnum_dnestro_re_7_3_m, only: fortnum_dnestro_re_7_3
    use fortnum_dnestro_re_7_4_m, only: fortnum_dnestro_re_7_4
    use fortnum_dnestro_re_7_5_m, only: fortnum_dnestro_re_7_5
    use fortnum_dnestro_re_7_6_m, only: fortnum_dnestro_re_7_6
    use fortnum_dnestro_re_7_7_m, only: fortnum_dnestro_re_7_7
    use fortnum_dnestro_re_7_8_m, only: fortnum_dnestro_re_7_8
    use fortnum_dnestro_tail_7_m, only: fortnum_dnestro_tail_7
    use fortnum_dnestro_im_7_m, only: fortnum_dnestro_im_7
    use fortnum_dnestro_pos_5_1_m, only: fortnum_dnestro_pos_5_1
    use fortnum_dnestro_pos_5_2_m, only: fortnum_dnestro_pos_5_2
    use fortnum_dnestro_pos_5_3_m, only: fortnum_dnestro_pos_5_3
    use fortnum_dnestro_pos_5_4_m, only: fortnum_dnestro_pos_5_4
    use fortnum_dnestro_pos_5_5_m, only: fortnum_dnestro_pos_5_5
    use fortnum_dnestro_pos_5_6_m, only: fortnum_dnestro_pos_5_6
    use fortnum_dnestro_pos_5_7_m, only: fortnum_dnestro_pos_5_7
    use fortnum_dnestro_pos_5_8_m, only: fortnum_dnestro_pos_5_8
    use fortnum_dnestro_pos_tail_5_m, only: fortnum_dnestro_pos_tail_5
    use fortnum_dnestro_pos_7_1_m, only: fortnum_dnestro_pos_7_1
    use fortnum_dnestro_pos_7_2_m, only: fortnum_dnestro_pos_7_2
    use fortnum_dnestro_pos_7_3_m, only: fortnum_dnestro_pos_7_3
    use fortnum_dnestro_pos_7_4_m, only: fortnum_dnestro_pos_7_4
    use fortnum_dnestro_pos_7_5_m, only: fortnum_dnestro_pos_7_5
    use fortnum_dnestro_pos_7_6_m, only: fortnum_dnestro_pos_7_6
    use fortnum_dnestro_pos_7_7_m, only: fortnum_dnestro_pos_7_7
    use fortnum_dnestro_pos_7_8_m, only: fortnum_dnestro_pos_7_8
    use fortnum_dnestro_pos_tail_7_m, only: fortnum_dnestro_pos_tail_7
    implicit none
    private
    public :: dnestrovskii
contains
    pure subroutine dnestrovskii(q_twice, z, value, derivative, status)
        integer, intent(in) :: q_twice
        real(dp), intent(in) :: z
        complex(dp), intent(out) :: value, derivative
        integer, intent(out) :: status
        real(dp), parameter :: edges(9) = &
                       [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp, 32.0_dp, 64.0_dp, 128.0_dp]
        real(dp) :: x, t, re, im, dr, di, scale
        integer :: j
        value = cmplx(0.0_dp, 0.0_dp, dp)
        derivative = value
        status = 1
        if (q_twice /= 5 .and. q_twice /= 7) return
        if (.not. ieee_is_finite(z)) return
        if (z > 0.0_dp) then
            call positive_axis(q_twice, z, re, dr)
            value = cmplx(re, 0.0_dp, dp)
            derivative = cmplx(dr, 0.0_dp, dp)
            status = 0
            return
        end if
        x = -z
        if (x <= edges(9)) then
            j = 1
            do while (j < 8)
                if (x <= edges(j + 1)) exit
                j = j + 1
            end do
            scale = 2.0_dp/(edges(j + 1) - edges(j))
            t = (x - edges(j))*scale - 1.0_dp
            select case (q_twice)
            case (5)
                select case (j)
                case (1); call fortnum_dnestro_re_5_1(t, re, dr)
                case (2); call fortnum_dnestro_re_5_2(t, re, dr)
                case (3); call fortnum_dnestro_re_5_3(t, re, dr)
                case (4); call fortnum_dnestro_re_5_4(t, re, dr)
                case (5); call fortnum_dnestro_re_5_5(t, re, dr)
                case (6); call fortnum_dnestro_re_5_6(t, re, dr)
                case (7); call fortnum_dnestro_re_5_7(t, re, dr)
                case (8); call fortnum_dnestro_re_5_8(t, re, dr)
                end select
            case (7)
                select case (j)
                case (1); call fortnum_dnestro_re_7_1(t, re, dr)
                case (2); call fortnum_dnestro_re_7_2(t, re, dr)
                case (3); call fortnum_dnestro_re_7_3(t, re, dr)
                case (4); call fortnum_dnestro_re_7_4(t, re, dr)
                case (5); call fortnum_dnestro_re_7_5(t, re, dr)
                case (6); call fortnum_dnestro_re_7_6(t, re, dr)
                case (7); call fortnum_dnestro_re_7_7(t, re, dr)
                case (8); call fortnum_dnestro_re_7_8(t, re, dr)
                end select
            end select
            dr = -scale*dr
        else
            t = 1.0_dp/x
            if (q_twice == 5) then
                call fortnum_dnestro_tail_5(t, re, dr)
            else
                call fortnum_dnestro_tail_7(t, re, dr)
            end if
            dr = dr*t*t
        end if
        ! The exponential tail is negligible before a power can overflow.
        im = 0.0_dp
        di = 0.0_dp
        if (x < 700.0_dp) then
            if (q_twice == 5) then
                call fortnum_dnestro_im_5(x, im, di)
            else
                call fortnum_dnestro_im_7(x, im, di)
            end if
            di = -di
        end if
        value = cmplx(re, im, dp)
        derivative = cmplx(dr, di, dp)
        status = 0
    end subroutine
    pure subroutine positive_axis(q, z, value, derivative)
        integer, intent(in) :: q
        real(dp), intent(in) :: z
        real(dp), intent(out) :: value, derivative
        real(dp), parameter :: edges(9) = &
                       [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp, 32.0_dp, 64.0_dp, 128.0_dp]
        real(dp) :: t, scale
        integer :: j
        if (z <= 1.0_dp) then
            if (q == 5) then
                call fortnum_dnestro_pos_5_1(z, value, derivative)
            else
                call fortnum_dnestro_pos_7_1(z, value, derivative)
            end if
        else if (z > 128.0_dp) then
            t = 1.0_dp/z
            if (q == 5) then
                call fortnum_dnestro_pos_tail_5(t, value, derivative)
            else
                call fortnum_dnestro_pos_tail_7(t, value, derivative)
            end if
            derivative = -derivative*t*t
        else
            j = 2
            do while (j < 8)
                if (z <= edges(j + 1)) exit
                j = j + 1
            end do
            scale = 2.0_dp/(edges(j + 1) - edges(j))
            t = (z - edges(j))*scale - 1.0_dp
            select case (q)
            case (5)
                select case (j)
                case (2); call fortnum_dnestro_pos_5_2(t, value, derivative)
                case (3); call fortnum_dnestro_pos_5_3(t, value, derivative)
                case (4); call fortnum_dnestro_pos_5_4(t, value, derivative)
                case (5); call fortnum_dnestro_pos_5_5(t, value, derivative)
                case (6); call fortnum_dnestro_pos_5_6(t, value, derivative)
                case (7); call fortnum_dnestro_pos_5_7(t, value, derivative)
                case (8); call fortnum_dnestro_pos_5_8(t, value, derivative)
                end select
            case (7)
                select case (j)
                case (2); call fortnum_dnestro_pos_7_2(t, value, derivative)
                case (3); call fortnum_dnestro_pos_7_3(t, value, derivative)
                case (4); call fortnum_dnestro_pos_7_4(t, value, derivative)
                case (5); call fortnum_dnestro_pos_7_5(t, value, derivative)
                case (6); call fortnum_dnestro_pos_7_6(t, value, derivative)
                case (7); call fortnum_dnestro_pos_7_7(t, value, derivative)
                case (8); call fortnum_dnestro_pos_7_8(t, value, derivative)
                end select
            end select
            derivative = derivative*scale
        end if
    end subroutine positive_axis
end module
