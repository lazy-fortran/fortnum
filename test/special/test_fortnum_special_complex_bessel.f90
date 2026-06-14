program test_fortnum_special_complex_bessel
    ! Behavioral tests for fortnum_special_complex_bessel: known reference
    ! values, the AMOS-style order-sequence convention KiLCA uses, the I(i z)
    ! relation to J, and the scaling-flag contract.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_special_complex_bessel, only: bessel_j_complex, &
        bessel_i_complex, bessel_i_complex_array, bessel_k_complex, &
        bessel_k_complex_array

    implicit none

    real(dp),    parameter :: tol = 1.0e-12_dp
    complex(dp), parameter :: ii  = (0.0_dp, 1.0_dp)

    integer :: nfail
    nfail = 0

    call test_known_values(nfail)
    call test_sequence_convention(nfail)
    call test_i_of_iz_is_j(nfail)
    call test_scaling_flag(nfail)
    call test_negative_order(nfail)
    call test_k_domain_guard(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine check(label, got, want, nfail)
        character(*), intent(in)    :: label
        complex(dp),  intent(in)    :: got, want
        integer,      intent(inout) :: nfail
        real(dp) :: e
        e = abs(got - want)/max(abs(want), 1.0e-300_dp)
        if (e > tol) then
            write (error_unit, "(a,a,a,2es24.16,a,2es24.16,a,es12.4)") &
                "FAIL ", label, " got=", got, " want=", want, " relerr=", e
            nfail = nfail + 1
        end if
    end subroutine check

    subroutine test_known_values(nfail)
        integer, intent(inout) :: nfail
        complex(dp) :: r
        type(fortnum_status_t) :: s
        ! J_0(1+i), I_0(1+i), K_0(1+i) from scipy.special (machine reference).
        call bessel_j_complex(0, (1.0_dp, 1.0_dp), r, s)
        call check("J0(1+i)", r, (0.9376084768060293_dp, -0.4965299476091221_dp), nfail)
        call bessel_i_complex(0, (1.0_dp, 1.0_dp), .false., r, s)
        call check("I0(1+i)", r, (0.9376084768060292_dp, 0.4965299476091221_dp), nfail)
        call bessel_k_complex(0, (1.0_dp, 1.0_dp), .false., r, s)
        call check("K0(1+i)", r, (0.08019772694651779_dp, -0.35727745928533017_dp), nfail)
    end subroutine test_known_values

    subroutine test_sequence_convention(nfail)
        ! KiLCA hom_medium calls zbesi(...,FNU,KODE,3,...) for orders
        ! FNU..FNU+2.  Sequence base order0=1 -> I_1, I_2, I_3.
        integer, intent(inout) :: nfail
        complex(dp) :: seq(3), single
        type(fortnum_status_t) :: s
        integer :: k
        call bessel_i_complex_array(1, 3, (2.0_dp, 0.5_dp), .false., seq, s)
        if (.not. status_ok(s)) then
            write (error_unit, "(a)") "FAIL sequence status"; nfail = nfail + 1
        end if
        do k = 1, 3
            call bessel_i_complex(k, (2.0_dp, 0.5_dp), .false., single, s)
            call check("seq vs single I", seq(k), single, nfail)
        end do
    end subroutine test_sequence_convention

    subroutine test_i_of_iz_is_j(nfail)
        ! KiLCA besseli builds I_n(z) = i^{-n} J_n(i z).  Verify the relation
        ! holds for our independent J and I implementations.
        integer, intent(inout) :: nfail
        complex(dp) :: zarg, jval, ival, recovered
        type(fortnum_status_t) :: s
        integer :: n
        zarg = (1.5_dp, -0.7_dp)
        do n = 0, 3
            call bessel_j_complex(n, ii*zarg, jval, s)
            call bessel_i_complex(n, zarg, .false., ival, s)
            recovered = ii**(-n)*jval
            call check("I_n = i^-n J_n(iz)", recovered, ival, nfail)
        end do
    end subroutine test_i_of_iz_is_j

    subroutine test_scaling_flag(nfail)
        ! KODE=2: scaled I = e^{-Re z} I; scaled K = e^{z} K.  Check against the
        ! unscaled values times the documented factor at moderate z.
        integer, intent(inout) :: nfail
        complex(dp) :: z, iun, isc, kun, ksc
        type(fortnum_status_t) :: s
        z = (3.0_dp, 1.0_dp)
        call bessel_i_complex(2, z, .false., iun, s)
        call bessel_i_complex(2, z, .true.,  isc, s)
        call check("I scaling", isc, exp(-real(z, dp))*iun, nfail)
        call bessel_k_complex(2, z, .false., kun, s)
        call bessel_k_complex(2, z, .true.,  ksc, s)
        call check("K scaling", ksc, exp(z)*kun, nfail)
    end subroutine test_scaling_flag

    subroutine test_negative_order(nfail)
        integer, intent(inout) :: nfail
        complex(dp) :: z, jm, jp, im, ip, km, kp
        type(fortnum_status_t) :: s
        z = (2.5_dp, 0.6_dp)
        ! J_{-n} = (-1)^n J_n
        call bessel_j_complex(-3, z, jm, s)
        call bessel_j_complex(3, z, jp, s)
        call check("J_-3 = -J_3", jm, -jp, nfail)
        ! I_{-n} = I_n
        call bessel_i_complex(-2, z, .false., im, s)
        call bessel_i_complex(2, z, .false., ip, s)
        call check("I_-2 = I_2", im, ip, nfail)
        ! K_{-n} = K_n
        call bessel_k_complex(-2, z, .false., km, s)
        call bessel_k_complex(2, z, .false., kp, s)
        call check("K_-2 = K_2", km, kp, nfail)
    end subroutine test_negative_order

    subroutine test_k_domain_guard(nfail)
        ! K requires Re z > 0; a non-positive real part must flag a status error.
        integer, intent(inout) :: nfail
        complex(dp) :: r
        type(fortnum_status_t) :: s
        call bessel_k_complex(0, (-1.0_dp, 1.0_dp), .false., r, s)
        if (status_ok(s)) then
            write (error_unit, "(a)") "FAIL K domain guard not raised"
            nfail = nfail + 1
        end if
    end subroutine test_k_domain_guard

end program test_fortnum_special_complex_bessel
