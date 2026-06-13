program test_enzyme_vs_analytic
    ! Enzyme-vs-analytic JVP compare. Compiles in every build; runs the
    ! comparison only when fortnum is built with the Enzyme pass available,
    ! and otherwise prints a skip line and exits 0. Until real generated
    ! products exist (#40) the "Enzyme" side is a stand-in analytic JVP, so the
    ! plumbing is exercised whenever Enzyme is present.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_ad_test_utils, only: enzyme_available, check_enzyme_vs_analytic
    implicit none

    integer :: nfail
    real(dp) :: x(3), v(3)

    if (.not. enzyme_available()) then
        write (*, '(a)') "SKIP enzyme_vs_analytic: Enzyme not available"
        stop 0
    end if

    nfail = 0
    x = [1.0_dp, -2.0_dp, 0.5_dp]
    v = [0.3_dp, 1.1_dp, -0.7_dp]
    if (.not. check_enzyme_vs_analytic("linear", jvp_enzyme, jvp_analytic, &
            x, v, tol=1.0e-12_dp)) nfail = nfail + 1

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    ! Stand-in for the generated entry point until #40 wires a real one.
    subroutine jvp_enzyme(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        jv(1) = 2.0_dp*v(1) - v(2) + 3.0_dp*v(3)
        jv(2) = v(1) + 4.0_dp*v(2)
        jv(3) = -v(1) + 0.5_dp*v(2) + 2.0_dp*v(3)
    end subroutine jvp_enzyme

    subroutine jvp_analytic(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        jv(1) = 2.0_dp*v(1) - v(2) + 3.0_dp*v(3)
        jv(2) = v(1) + 4.0_dp*v(2)
        jv(3) = -v(1) + 0.5_dp*v(2) + 2.0_dp*v(3)
    end subroutine jvp_analytic

end program test_enzyme_vs_analytic
