program test_stable_sqrt_difference
    use, intrinsic :: iso_fortran_env, only: real64, real128
    use fortnum_generated_sqrt1pm1_raw, only: fortnum_sqrt1pm1_raw_kernel
    use fortnum_generated_sqrt1pm1_stable, only: fortnum_sqrt1pm1_stable_kernel
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: qp = real128
    real(dp) :: x, raw, stable
    real(qp) :: reference, raw_error, stable_error
    integer :: nfail

    nfail = 0
    x = nearest(1.0_dp, 1.0_dp) - 1.0_dp
    call fortnum_sqrt1pm1_raw_kernel(x, raw)
    call fortnum_sqrt1pm1_stable_kernel(x, stable)
    reference = sqrt(1.0_qp + real(x, qp)) - 1.0_qp
    raw_error = abs(real(raw, qp) - reference)
    stable_error = abs(real(stable, qp) - reference)
    call check("stable boundary beats raw", stable_error < raw_error)
    call check("stable boundary stays nonzero", &
        stable > 0.0_dp .and. raw == 0.0_dp)

    x = 0.75_dp
    call fortnum_sqrt1pm1_raw_kernel(x, raw)
    call fortnum_sqrt1pm1_stable_kernel(x, stable)
    reference = sqrt(1.0_qp + real(x, qp)) - 1.0_qp
    call check("stable ordinary-region accuracy", &
        abs(real(stable, qp) - reference) <= &
        2.0_qp*real(epsilon(1.0_dp), qp)*abs(reference))
    call check("raw and stable ordinary-region agreement", &
        abs(raw - stable) <= 2.0_dp*epsilon(1.0_dp))

    if (nfail /= 0) error stop 1
    print *, "PASS stable square-root difference"

contains

    subroutine check(label, condition)
        character(*), intent(in) :: label
        logical, intent(in) :: condition
        if (.not. condition) then
            nfail = nfail + 1
            print *, "FAIL ", label
        end if
    end subroutine check

end program test_stable_sqrt_difference
