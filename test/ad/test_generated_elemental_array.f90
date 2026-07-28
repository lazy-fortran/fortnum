program test_generated_elemental_array
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_generated_dawson_identity_jvp, only: &
        fortnum_dawson_identity_jvp_kernel
    implicit none

    real(dp), parameter :: x(4) = [-1.5_dp, -0.25_dp, 0.5_dp, 2.0_dp]
    real(dp), parameter :: d(4) = [0.2_dp, -0.7_dp, 1.1_dp, -1.4_dp]
    real(dp), parameter :: tx(4) = [0.3_dp, 0.8_dp, -0.4_dp, 0.6_dp]
    real(dp), parameter :: td(4) = [-0.9_dp, 0.1_dp, 0.5_dp, -0.2_dp]
    real(dp) :: actual(4), expected(4), scalar

    call fortnum_dawson_identity_jvp_kernel(x, d, tx, td, actual)
    expected = d*tx + td*(x - 2.0_dp*d*exp(-d*d))
    if (maxval(abs(actual - expected)) > 32.0_dp*epsilon(1.0_dp)) then
        error stop "elemental array result disagrees with analytical oracle"
    end if

    call fortnum_dawson_identity_jvp_kernel(x(3), d(3), tx(3), td(3), scalar)
    if (abs(scalar - expected(3)) > 32.0_dp*epsilon(1.0_dp)) then
        error stop "elemental scalar result disagrees with analytical oracle"
    end if
end program test_generated_elemental_array
