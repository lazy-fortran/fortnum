program test_generated_hyperg_region_array
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_generated_hyperg_asymptotic_outer_jvp, only: &
        fortnum_hyperg_asymptotic_outer_jvp
    use fortnum_generated_hyperg_asymptotic_outer_vjp, only: &
        fortnum_hyperg_asymptotic_outer_vjp
    implicit none

    real(dp), parameter :: gamma_25 = 1.3293403881791370205_dp
    real(dp), parameter :: x(4) = [64.0_dp, 80.0_dp, 120.0_dp, 160.0_dp]
    real(dp), parameter :: seed(4) = [0.25_dp, -0.5_dp, 1.25_dp, -1.5_dp]
    real(dp) :: actual(4), expected(4), m(4), scalar

    m = gamma_25*exp(x)*x**(-1.5_dp)
    expected = 2.0_dp*seed*(1.0_dp - 1.5_dp/x)/(1.0_dp + 1.0_dp/(m*m))

    call fortnum_hyperg_asymptotic_outer_jvp(x, seed, actual)
    call require_close("elemental array JVP", actual, expected)

    call fortnum_hyperg_asymptotic_outer_vjp(x, seed, actual)
    call require_close("elemental array VJP", actual, expected)

    call fortnum_hyperg_asymptotic_outer_jvp(x(2), seed(2), scalar)
    if (abs(scalar - expected(2)) > 64.0_dp*epsilon(1.0_dp)*abs(expected(2))) &
        error stop "elemental scalar JVP disagrees with analytical oracle"

contains

    subroutine require_close(label, got, want)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got(:), want(:)
        real(dp) :: scale

        scale = max(1.0_dp, maxval(abs(want)))
        if (maxval(abs(got - want)) > 64.0_dp*epsilon(1.0_dp)*scale) then
            write (*, "(a,1x,es12.4)") trim(label), maxval(abs(got - want))
            error stop "generated array product disagrees with analytical oracle"
        end if
    end subroutine require_close

end program test_generated_hyperg_region_array
