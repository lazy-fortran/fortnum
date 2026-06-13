program test_fortnum_special_dawson
    ! Behavioral unit tests for the Dawson integral module.
    ! Each test exercises a distinct property; no external reference needed.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special_dawson, only: dawson
    implicit none

    logical :: ok

    ok = .true.
    call test_odd(ok)
    call test_zero(ok)
    call test_known_values(ok)
    call test_elemental(ok)

    if (.not. ok) stop 1
    write (*, "(a)") "dawson unit tests passed"
    stop 0

contains

    subroutine test_odd(ok)
        ! F is an odd function: F(-x) = -F(x).
        logical, intent(inout) :: ok
        real(dp), parameter :: xs(*) = [0.3_dp, 1.5_dp, 7.0_dp, 12.0_dp]
        integer :: i
        real(dp) :: err
        do i = 1, size(xs)
            err = abs(dawson(-xs(i)) + dawson(xs(i)))
            if (err > 1.0e-15_dp) then
                write (error_unit, "(a,es12.4,a,es12.4)") &
                    "FAIL test_odd at x=", xs(i), " err=", err
                ok = .false.
            end if
        end do
    end subroutine test_odd

    subroutine test_zero(ok)
        ! F(0) = 0 exactly.
        logical, intent(inout) :: ok
        if (dawson(0.0_dp) /= 0.0_dp) then
            write (error_unit, "(a)") "FAIL test_zero: dawson(0) /= 0"
            ok = .false.
        end if
    end subroutine test_zero

    subroutine test_known_values(ok)
        ! Spot-check a few values against the scipy reference to cover the
        ! three algorithm branches: series (|x|<1), Rybicki (1<=|x|<10),
        ! and asymptotic (|x|>=10).
        logical, intent(inout) :: ok
        real(dp), parameter :: ref_x(*) = [0.5_dp, 1.0_dp, 5.0_dp, 10.0_dp]
        real(dp), parameter :: ref_f(*) = [ &
            0.42443638350202226_dp, &   ! scipy.special.dawsn(0.5)
            0.5380795069127684_dp, &    ! scipy.special.dawsn(1.0)
            0.10213407442427686_dp, &   ! scipy.special.dawsn(5.0)
            0.050253847187598542_dp]    ! scipy.special.dawsn(10.0)
        real(dp), parameter :: rtol = 1.0e-13_dp
        integer :: i
        real(dp) :: err
        do i = 1, size(ref_x)
            err = abs(dawson(ref_x(i)) - ref_f(i)) / abs(ref_f(i))
            if (err > rtol) then
                write (error_unit, "(a,es12.4,a,es12.4,a,es12.4)") &
                    "FAIL test_known_values at x=", ref_x(i), &
                    " got=", dawson(ref_x(i)), " relerr=", err
                ok = .false.
            end if
        end do
    end subroutine test_known_values

    subroutine test_elemental(ok)
        ! elemental: dawson applied to an array must equal element-wise calls.
        logical, intent(inout) :: ok
        real(dp) :: xs(5), ys(5)
        integer :: i
        xs = [0.1_dp, 0.9_dp, 2.0_dp, 8.5_dp, 11.0_dp]
        ys = dawson(xs)
        do i = 1, size(xs)
            if (ys(i) /= dawson(xs(i))) then
                write (error_unit, "(a,i0)") "FAIL test_elemental at index ", i
                ok = .false.
            end if
        end do
    end subroutine test_elemental

end program test_fortnum_special_dawson
