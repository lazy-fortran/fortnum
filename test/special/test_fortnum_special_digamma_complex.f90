program test_fortnum_special_digamma_complex
    ! Complex digamma/trigamma against closed forms that share no formula
    ! with the shifted Bernoulli series:
    !   psi(1) = -gamma, psi(1/2) = -gamma - 2 ln 2 (DLMF 5.4.12, 5.4.13),
    !   Im psi(iy), Im psi(1/2+iy), Im psi(1+iy) (DLMF 5.4.17-5.4.19),
    !   Re psi(1+iy) = -gamma + y^2 sum 1/(n(n^2+y^2)) (DLMF 5.15.1 summed
    !   in quad precision with an Euler-Maclaurin tail),
    !   psi'(1) = pi^2/6, psi'(1/2) = pi^2/2 (DLMF 5.15.2, 5.15.3),
    !   Re psi'(1+iy) = d/dy Im psi(1+iy) = 1/(2y^2) - (pi^2/2) csch^2(pi y),
    ! plus the recurrences psi(z+1) = psi(z) + 1/z and
    ! psi'(z+1) = psi'(z) - 1/z^2 across both half planes.
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_special, only: digamma_complex, trigamma_complex
    implicit none

    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp), parameter :: euler = 0.577215664901532860606512090082402431_dp
    real(dp), parameter :: tol = 4.0e-15_dp
    real(dp), parameter :: ys(8) = [1.0e-3_dp, 0.1_dp, 0.7_dp, 1.0_dp, &
        3.3_dp, 10.0_dp, 47.0_dp, 400.0_dp]
    complex(dp), parameter :: zs(8) = [(0.3_dp, 0.2_dp), (2.5_dp, -1.0_dp), &
        (-3.7_dp, 0.4_dp), (-0.5_dp, -6.0_dp), (12.0_dp, 30.0_dp), &
        (-20.3_dp, -0.01_dp), (0.0_dp, 1.0e-4_dp), (1.0e3_dp, -2.0e2_dp)]
    integer :: nfail, i
    real(dp) :: y
    complex(dp) :: z, iu

    nfail = 0
    iu = (0.0_dp, 1.0_dp)
    call check("psi(1)", digamma_complex((1.0_dp, 0.0_dp)), &
        cmplx(-euler, 0.0_dp, dp))
    call check("psi(1/2)", digamma_complex((0.5_dp, 0.0_dp)), &
        cmplx(-euler - 2.0_dp*log(2.0_dp), 0.0_dp, dp))
    call check("psi(-1/2)", digamma_complex((-0.5_dp, 0.0_dp)), &
        cmplx(2.0_dp - euler - 2.0_dp*log(2.0_dp), 0.0_dp, dp), 1.0_dp)
    call check("psi'(1)", trigamma_complex((1.0_dp, 0.0_dp)), &
        cmplx(pi*pi/6.0_dp, 0.0_dp, dp))
    call check("psi'(1/2)", trigamma_complex((0.5_dp, 0.0_dp)), &
        cmplx(pi*pi/2.0_dp, 0.0_dp, dp))

    do i = 1, size(ys)
        y = ys(i)
        call check_real("Im psi(iy)", aimag(digamma_complex(iu*y)), &
            0.5_dp/y + 0.5_dp*pi/tanh(pi*y))
        call check_real("Im psi(1/2+iy)", &
            aimag(digamma_complex(cmplx(0.5_dp, y, dp))), 0.5_dp*pi*tanh(pi*y))
        call check_real("Im psi(1+iy)", &
            aimag(digamma_complex(cmplx(1.0_dp, y, dp))), &
            -0.5_dp/y + 0.5_dp*pi/tanh(pi*y), abs(0.5_dp*pi/tanh(pi*y)))
        call check_real("Re psi(1+iy)", &
            real(digamma_complex(cmplx(1.0_dp, y, dp))), re_psi_one(y), &
            abs(digamma_complex(cmplx(1.0_dp, y, dp))))
        call check_real("Re psi'(1+iy)", &
            real(trigamma_complex(cmplx(1.0_dp, y, dp))), &
            re_trigamma_one(y), 0.5_dp/(y*y))
    end do

    do i = 1, size(zs)
        z = zs(i)
        call check("psi recurrence", digamma_complex(z + 1.0_dp), &
            digamma_complex(z) + 1.0_dp/z, &
            2.5_dp*(abs(digamma_complex(z)) + abs(1.0_dp/z)))
        call check("psi' recurrence", trigamma_complex(z + 1.0_dp), &
            trigamma_complex(z) - 1.0_dp/(z*z), &
            2.5_dp*(abs(trigamma_complex(z)) + abs(1.0_dp/(z*z))))
    end do

    if (ieee_is_finite(real(digamma_complex((0.0_dp, 0.0_dp))))) &
        call fail("pole at 0")
    if (ieee_is_finite(real(digamma_complex((-3.0_dp, 0.0_dp))))) &
        call fail("pole at -3")
    if (ieee_is_finite(real(trigamma_complex((-2.0_dp, 0.0_dp))))) &
        call fail("trigamma pole at -2")

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) FAILED"
        error stop 1
    end if
    write (*, "(a)") "PASS"

contains

    function re_psi_one(yd) result(value)
        ! -gamma + sum_{n>=1} y^2/(n(n^2+y^2)) in quad precision.
        real(dp), intent(in) :: yd
        real(dp) :: value
        integer, parameter :: nmax = 200000
        real(qp) :: y2, total, fn
        integer :: n

        y2 = real(yd, qp)**2
        total = 0.0_qp
        do n = nmax, 1, -1
            total = total + y2/(real(n, qp)*(real(n, qp)**2 + y2))
        end do
        fn = y2/(real(nmax, qp)*(real(nmax, qp)**2 + y2))
        total = total + 0.5_qp*log(1.0_qp + y2/real(nmax, qp)**2) - 0.5_qp*fn
        value = real(total, dp) - euler
    end function re_psi_one

    function re_trigamma_one(yd) result(value)
        real(dp), intent(in) :: yd
        real(dp) :: value
        real(qp) :: yq, s

        yq = real(yd, qp)
        s = sinh(acos(-1.0_qp)*yq)
        value = real(0.5_qp/(yq*yq) - 0.5_qp*acos(-1.0_qp)**2/(s*s), dp)
    end function re_trigamma_one

    subroutine check(label, got, expected, scale)
        character(*), intent(in) :: label
        complex(dp), intent(in) :: got, expected
        real(dp), intent(in), optional :: scale
        real(dp) :: s

        s = abs(expected)
        if (present(scale)) s = scale
        if (.not. abs(got - expected) <= tol*max(s, 1.0e-300_dp)) then
            nfail = nfail + 1
            write (error_unit, "(a,4es24.16)") label, got, expected
        end if
    end subroutine check

    subroutine check_real(label, got, expected, scale)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected
        real(dp), intent(in), optional :: scale
        real(dp) :: s

        s = abs(expected)
        if (present(scale)) s = scale
        if (.not. abs(got - expected) <= tol*s) then
            nfail = nfail + 1
            write (error_unit, "(a,' y=',es10.3,2es24.16)") label, y, got, &
                expected
        end if
    end subroutine check_real

    subroutine fail(label)
        character(*), intent(in) :: label
        nfail = nfail + 1
        write (error_unit, "(a)") label
    end subroutine fail

end program test_fortnum_special_digamma_complex
