program test_fortnum_special_reexport
    ! Verify that fortnum_special re-exports all public names from
    ! the constituent special modules (bessel, dawson, gamma).
    !
    ! Tests one routine from each module to ensure the re-export
    ! names resolve and produce correct results.

    use fortnum_kinds, only: dp
    use fortnum_special, only: &
        bessel_in, bessel_in_array, bessel_kn, &
        dawson, &
        gamma_lower, gamma_reg_p

    implicit none

    integer :: nfail
    real(dp) :: val, tol

    nfail = 0
    tol = 1.0e-14_dp

    ! Test bessel_in re-export: I_0(1.0)
    val = bessel_in(0, 1.0_dp)
    if (abs(val - 1.2660658777520082_dp) > tol * abs(val) + 1.0e-15_dp) then
        print *, "FAIL: bessel_in(0, 1.0)"
        nfail = nfail + 1
    else
        print *, "PASS: bessel_in(0, 1.0) =", val
    end if

    ! Test bessel_kn re-export: K_0(1.0)
    val = bessel_kn(0, 1.0_dp)
    if (abs(val - 0.4210244382407084_dp) > tol * abs(val) + 1.0e-15_dp) then
        print *, "FAIL: bessel_kn(0, 1.0)"
        nfail = nfail + 1
    else
        print *, "PASS: bessel_kn(0, 1.0) =", val
    end if

    ! Test bessel_in_array re-export: I_0..I_2(1.0)
    block
        real(dp) :: vals(0:2)
        call bessel_in_array(2, 1.0_dp, vals)
        if (abs(vals(0) - 1.2660658777520082_dp) > tol * abs(vals(0)) + 1.0e-15_dp .or. &
            abs(vals(1) - 0.56515910399248503_dp) > tol * abs(vals(1)) + 1.0e-15_dp .or. &
            abs(vals(2) - 0.13574766976703831_dp) > tol * abs(vals(2)) + 1.0e-15_dp) then
            print *, "FAIL: bessel_in_array(2, 1.0)"
            nfail = nfail + 1
        else
            print *, "PASS: bessel_in_array(2, 1.0)"
        end if
    end block

    ! Test dawson re-export: F(1.0)
    val = dawson(1.0_dp)
    if (abs(val - 0.53807950691276851_dp) > tol * abs(val) + 1.0e-15_dp) then
        print *, "FAIL: dawson(1.0)"
        nfail = nfail + 1
    else
        print *, "PASS: dawson(1.0) =", val
    end if

    ! Test gamma_lower re-export: gamma_lower(2.0, 1.0)
    val = gamma_lower(2.0_dp, 1.0_dp)
    if (abs(val - 0.26424111765711528_dp) > tol * abs(val) + 1.0e-15_dp) then
        print *, "FAIL: gamma_lower(2.0, 1.0)"
        nfail = nfail + 1
    else
        print *, "PASS: gamma_lower(2.0, 1.0) =", val
    end if

    ! Test gamma_reg_p re-export: P(2.0, 1.0)
    val = gamma_reg_p(2.0_dp, 1.0_dp)
    if (abs(val - 0.26424111765711528_dp) > tol * abs(val) + 1.0e-15_dp) then
        print *, "FAIL: gamma_reg_p(2.0, 1.0)"
        nfail = nfail + 1
    else
        print *, "PASS: gamma_reg_p(2.0, 1.0) =", val
    end if

    if (nfail == 0) then
        print *, ""
        print *, "All tests passed!"
    else
        print *, ""
        print *, nfail, "test(s) failed"
        stop 1
    end if

end program test_fortnum_special_reexport
