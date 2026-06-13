program test_fortnum_special_bessel
    ! Behavioral tests for fortnum_special_bessel.
    ! Exercises: symmetry, recurrence identity, small-x limits, array fill.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special_bessel, only: bessel_in, bessel_in_array, bessel_kn

    implicit none

    integer  :: nfail
    real(dp) :: tol

    nfail = 0
    tol   = 1.0e-13_dp

    call test_in_symmetry(nfail, tol)
    call test_in_zero_argument(nfail, tol)
    call test_in_recurrence(nfail, tol)
    call test_in_array_matches_scalar(nfail, tol)
    call test_kn_symmetry(nfail, tol)
    call test_kn_recurrence(nfail, tol)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine check(label, got, expected, tol, nfail)
        character(*), intent(in)    :: label
        real(dp),     intent(in)    :: got, expected, tol
        integer,      intent(inout) :: nfail

        real(dp) :: scale, err

        scale = max(abs(expected), 1.0e-280_dp)
        err   = abs(got - expected)/scale
        if (err > tol) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es13.6,a,es24.16,a,es24.16)") &
                "FAIL [", label, "] relerr=", err, &
                " got=", got, " expected=", expected
        end if
    end subroutine check

    subroutine test_in_symmetry(nfail, tol)
        ! I_{-n}(x) = I_n(x); I_n(-x) = (-1)^n I_n(x)
        integer,  intent(inout) :: nfail
        real(dp), intent(in)    :: tol

        integer  :: n
        real(dp) :: x, pos, neg

        x = 3.5_dp
        do n = 0, 5
            pos = bessel_in( n, x)
            neg = bessel_in(-n, x)
            call check("I_{-n}=I_n", neg, pos, tol, nfail)
            pos = bessel_in(n,  x)
            neg = bessel_in(n, -x)
            if (mod(n, 2) == 0) then
                call check("I_n(-x)=I_n(x) even", neg, pos, tol, nfail)
            else
                call check("I_n(-x)=-I_n(x) odd", neg, -pos, tol, nfail)
            end if
        end do
    end subroutine test_in_symmetry

    subroutine test_in_zero_argument(nfail, tol)
        ! I_0(0) = 1; I_n(0) = 0 for n > 0
        integer,  intent(inout) :: nfail
        real(dp), intent(in)    :: tol

        integer :: n

        call check("I_0(0)=1", bessel_in(0, 0.0_dp), 1.0_dp, tol, nfail)
        do n = 1, 5
            call check("I_n(0)=0", bessel_in(n, 0.0_dp), 0.0_dp, tol, nfail)
        end do
    end subroutine test_in_zero_argument

    subroutine test_in_recurrence(nfail, tol)
        ! DLMF 10.29.1: I_{n-1} - I_{n+1} = (2n/x) I_n
        integer,  intent(inout) :: nfail
        real(dp), intent(in)    :: tol

        integer  :: n
        real(dp) :: x

        x = 5.0_dp
        do n = 1, 4
            call check("recurrence", &
                bessel_in(n - 1, x) - bessel_in(n + 1, x), &
                (2.0_dp*real(n, dp)/x)*bessel_in(n, x), tol, nfail)
        end do
    end subroutine test_in_recurrence

    subroutine test_in_array_matches_scalar(nfail, tol)
        ! bessel_in_array must agree with bessel_in for each order
        integer,  intent(inout) :: nfail
        real(dp), intent(in)    :: tol

        integer,  parameter :: nmax = 10
        real(dp) :: x, arr(0:nmax)
        integer  :: n
        character(len=32) :: label

        x = 7.3_dp
        call bessel_in_array(nmax, x, arr)
        do n = 0, nmax
            write (label, "(a,i0)") "array_vs_scalar_n=", n
            call check(trim(label), arr(n), bessel_in(n, x), tol, nfail)
        end do
    end subroutine test_in_array_matches_scalar

    subroutine test_kn_symmetry(nfail, tol)
        ! K_{-n}(x) = K_n(x)
        integer,  intent(inout) :: nfail
        real(dp), intent(in)    :: tol

        integer  :: n
        real(dp) :: x

        x = 2.5_dp
        do n = 0, 5
            call check("K_{-n}=K_n", bessel_kn(-n, x), bessel_kn(n, x), tol, nfail)
        end do
    end subroutine test_kn_symmetry

    subroutine test_kn_recurrence(nfail, tol)
        ! DLMF 10.29.1: K_{n+1}(x) = K_{n-1}(x) + (2n/x) K_n(x)
        integer,  intent(inout) :: nfail
        real(dp), intent(in)    :: tol

        integer  :: n
        real(dp) :: x

        x = 4.0_dp
        do n = 1, 4
            call check("K recurrence", &
                bessel_kn(n + 1, x), &
                bessel_kn(n - 1, x) + (2.0_dp*real(n, dp)/x)*bessel_kn(n, x), &
                tol, nfail)
        end do
    end subroutine test_kn_recurrence

end program test_fortnum_special_bessel
