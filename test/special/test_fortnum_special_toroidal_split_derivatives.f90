program test_fortnum_special_toroidal_split_derivatives
    ! dp_n/dx and dr_n/dx of the logarithmic split Q = L p + r.  Oracles:
    ! mpmath legenp/legenq (type=3, 60 digits, exact double x, derivative by
    ! a relative step of 1e-20 at working precision 60); centred finite
    ! differences of p and r away from the cut; continuity across the switch
    ! between the Taylor series and the recurrence-based formulas.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special_toroidal, only: toroidal_q_log_split, &
        toroidal_q_log_split_batch
    implicit none

    integer, parameter :: nx = 7, nns = 5, nmax = 10
    real(dp), parameter :: x_ref(nx) = [ &
        1.0000000000000002220446049250313_dp, 1.000000000001_dp, &
        1.000001_dp, 1.001_dp, 1.05_dp, 1.5_dp, 3.0_dp]
    integer, parameter :: n_ref(nns) = [0, 1, 2, 5, 10]
    real(dp), parameter :: dp_ref(5, 7) = reshape([ &
        -0.12499999999999998439_dp, &
        0.37499999999999997398_dp, &
        1.8750000000000001821_dp, &
        12.375000000000015628_dp, &
        49.875000000000270632_dp, &
        -0.12499999999992968125_dp, &
        0.37499999999988280208_dp, &
        1.8750000000008203854_dp, &
        12.37500000007038907_dp, &
        49.875000001218928666_dp, &
        -0.12499992968753662686_dp, &
        0.37499988281255127915_dp, &
        1.875000820312346124_dp, &
        12.375070382922467413_dp, &
        49.876218829834468237_dp, &
        -0.1249297241024111992_dp, &
        0.37488286374551023017_dp, &
        1.8758201587442563941_dp, &
        12.445492844074693658_dp, &
        51.103377230619979402_dp, &
        -0.12157364892069505512_dp, &
        0.36926586531982249378_dp, &
        1.9156375811053397048_dp, &
        16.176419563667198196_dp, &
        139.43289363032492382_dp, &
        -0.09713482677963942655_dp, &
        0.32680092529541988718_dp, &
        2.2522100321114376027_dp, &
        82.778949928521087558_dp, &
        15304.39472320492156_dp, &
        -0.056655869990395824031_dp, &
        0.24734581086584912105_dp, &
        3.1381173403613769247_dp, &
        1157.2399841367271281_dp, &
        11512101.600260730419_dp &
        ], [5, 7])
    real(dp), parameter :: dr_ref(5, 7) = reshape([ &
        -0.048286795139993639384_dp, &
        0.39486038541995165742_dp, &
        -0.27569807290021218473_dp, &
        -14.455321566855656032_dp, &
        -93.525807872561434727_dp, &
        -0.048286795139970883409_dp, &
        0.39486038541991761577_dp, &
        -0.27569807289921233809_dp, &
        -14.455321566898025684_dp, &
        -93.525807874218630384_dp, &
        -0.04828677969741966398_dp, &
        0.39486034405734931365_dp, &
        -0.27569708023692856065_dp, &
        -14.455363949480699316_dp, &
        -93.527465305730273592_dp, &
        -0.048271358167129546184_dp, &
        0.39481903356786359305_dp, &
        -0.27470547603892818689_dp, &
        -14.49771905630720724_dp, &
        -95.192776960480565956_dp, &
        -0.04752840208321740935_dp, &
        0.39281874479886345062_dp, &
        -0.22622946267681755263_dp, &
        -16.608799791777872392_dp, &
        -203.49038228221484853_dp, &
        -0.041738370392008242364_dp, &
        0.37649405050447046008_dp, &
        0.20550657890004823971_dp, &
        -35.876776939013005836_dp, &
        -8771.7881525305794531_dp, &
        -0.030105651395713704549_dp, &
        0.33857395442430554285_dp, &
        1.5148490913764170033_dp, &
        183.24023386786453094_dp, &
        859775.41803914496647_dp &
        ], [5, 7])

    integer :: nfail = 0

    call test_mpmath_reference()
    call test_finite_differences()
    call test_switch_continuity()
    call test_batch_and_domain()

    if (nfail > 0) then
        write (error_unit, "(a,i0,a)") "FAILED: ", nfail, " checks"
        error stop 1
    end if
    print "(a)", "toroidal split derivative tests passed"

contains

    subroutine expect(label, value, reference, tol)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: value, reference, tol

        if (.not. (abs(value - reference) <= tol*abs(reference))) then
            nfail = nfail + 1
            write (error_unit, "(a,a,3es24.15)") "FAIL: ", label, value, &
                reference, abs(value - reference)/abs(reference)
        end if
    end subroutine expect

    subroutine test_mpmath_reference()
        ! The smallest x - 1 = 2.2e-16 (t ~ 1e-16) is where difference
        ! formulas lose every digit; x = 1.05 and above use the recurrence
        ! formulas for nmax = 10.
        integer :: i, j
        real(dp) :: p(0:nmax), r(0:nmax), dpdx(0:nmax), drdx(0:nmax)
        character(len=48) :: label

        do i = 1, nx
            call toroidal_q_log_split(nmax, x_ref(i), p, r, dpdx, drdx)
            do j = 1, nns
                write (label, "(a,i0,a,es9.2)") "dp/dx n=", n_ref(j), &
                    " x-1=", x_ref(i) - 1.0_dp
                call expect(trim(label), dpdx(n_ref(j)), dp_ref(j, i), &
                    1.0e-12_dp)
                write (label, "(a,i0,a,es9.2)") "dr/dx n=", n_ref(j), &
                    " x-1=", x_ref(i) - 1.0_dp
                call expect(trim(label), drdx(n_ref(j)), dr_ref(j, i), &
                    1.0e-12_dp)
            end do
        end do
    end subroutine test_mpmath_reference

    subroutine test_finite_differences()
        integer :: i, n
        real(dp) :: x, h, pa(0:nmax), ra(0:nmax), pb(0:nmax), rb(0:nmax)
        real(dp) :: pc(0:nmax), rc(0:nmax), pd(0:nmax), rd(0:nmax)
        real(dp) :: p(0:nmax), r(0:nmax), dpdx(0:nmax), drdx(0:nmax)
        real(dp), parameter :: xs(5) = [1.001_dp, 1.05_dp, 1.5_dp, 3.0_dp, &
            8.0_dp]

        do i = 1, size(xs)
            x = xs(i)
            h = 1.0e-3_dp*(x - 1.0_dp)
            call toroidal_q_log_split(nmax, x - 2.0_dp*h, pa, ra)
            call toroidal_q_log_split(nmax, x - h, pb, rb)
            call toroidal_q_log_split(nmax, x + h, pc, rc)
            call toroidal_q_log_split(nmax, x + 2.0_dp*h, pd, rd)
            call toroidal_q_log_split(nmax, x, p, r, dpdx, drdx)
            do n = 0, nmax
                call expect("dp/dx vs finite difference", dpdx(n), &
                    (pa(n) - 8.0_dp*pb(n) + 8.0_dp*pc(n) - pd(n))/(12.0_dp*h), &
                    1.0e-8_dp)
                call expect("dr/dx vs finite difference", drdx(n), &
                    (ra(n) - 8.0_dp*rb(n) + 8.0_dp*rc(n) - rd(n))/(12.0_dp*h), &
                    1.0e-8_dp)
            end do
        end do
    end subroutine test_finite_differences

    subroutine test_switch_continuity()
        ! Series for t < min(0.2, 2/(nmax+1)^2), recurrence formulas above:
        ! both sides must agree to the first-order change over the gap.
        integer, parameter :: ns(4) = [0, 1, 3, 10]
        integer :: i, n, m
        real(dp) :: xs, d, p(0:10), r(0:10)
        real(dp) :: dpa(0:10), dra(0:10), dpb(0:10), drb(0:10)

        do i = 1, size(ns)
            m = ns(i)
            xs = 1.0_dp + 2.0_dp*min(0.2_dp, 2.0_dp/real(m + 1, dp)**2)
            d = 1.0e-9_dp*(xs - 1.0_dp)
            call toroidal_q_log_split(m, xs - d, p(0:m), r(0:m), dpa(0:m), &
                dra(0:m))
            call toroidal_q_log_split(m, xs + d, p(0:m), r(0:m), dpb(0:m), &
                drb(0:m))
            do n = 0, m
                call expect("dp/dx switch continuity", dpb(n), dpa(n), &
                    1.0e-7_dp)
                call expect("dr/dx switch continuity", drb(n), dra(n), &
                    1.0e-7_dp)
            end do
        end do
    end subroutine test_switch_continuity

    subroutine test_batch_and_domain()
        integer, parameter :: np = 100
        real(dp) :: x(np), p(0:6, np), r(0:6, np), dpdx(0:6, np), drdx(0:6, np)
        real(dp) :: p1(0:6), r1(0:6), dp1(0:6), dr1(0:6), p0(0:6), r0(0:6)
        real(dp) :: dpo(0:6, np), nanv(0:2), nanr(0:2), nand(0:2)
        integer :: i

        do i = 1, np
            x(i) = 1.0_dp + 10.0_dp**(-15.0_dp + 15.0_dp*real(i - 1, dp)/ &
                real(np - 1, dp))
        end do
        call toroidal_q_log_split_batch(6, x, p, r, dpdx, drdx)
        call toroidal_q_log_split_batch(6, x, p, r, dpdx=dpo)
        do i = 1, np
            call toroidal_q_log_split(6, x(i), p1, r1, dp1, dr1)
            call toroidal_q_log_split(6, x(i), p0, r0)
            if (any(dpdx(:, i) /= dp1) .or. any(drdx(:, i) /= dr1) .or. &
                any(dpo(:, i) /= dp1) .or. any(p(:, i) /= p0) .or. &
                any(r(:, i) /= r0)) then
                nfail = nfail + 1
                write (error_unit, "(a,i0)") "FAIL: split derivative batch ", i
            end if
        end do
        call toroidal_q_log_split(2, 1.0_dp, nanv, nanr, nand, nand)
        if (.not. all(ieee_is_nan(nand))) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: x = 1 derivatives must be NaN"
        end if
    end subroutine test_batch_and_domain

end program test_fortnum_special_toroidal_split_derivatives
