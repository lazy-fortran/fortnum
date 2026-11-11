program test_fortnum_matrix_exponential
    use fortnum_kinds, only: dp
    use fortnum_matrix_exponential, only: matrix_exponential
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)
    complex(dp) :: a(3, 3), e(3, 3), expected(3, 3), s(3, 3), sinv(3, 3), lam(3)
    complex(dp) :: c4(4, 4), e4(4, 4), f4(4, 4)
    real(dp) :: r(2, 2), er(2, 2), rr(3, 3), err(3, 3), theta
    integer :: info, i, j, failures

    failures = 0

    ! Diagonal: exp acts on each entry.
    lam = [(-3.0_dp, 0.7_dp), (0.2_dp, -40.0_dp), (1.5_dp, 0.0_dp)]
    a = 0
    expected = 0
    do i = 1, 3
        a(i, i) = lam(i)
        expected(i, i) = exp(lam(i))
    end do
    call matrix_exponential(a, e, info)
    call check('diagonal', info, relative(e, expected), 1.0e-13_dp)

    ! Jordan block: exp(l I + N) = e^l (I + N + N^2/2), N nilpotent.
    a = 0
    do i = 1, 3
        a(i, i) = (0.3_dp, 2.0_dp)
    end do
    a(1, 2) = 5; a(2, 3) = -7
    expected = 0
    do i = 1, 3
        expected(i, i) = 1
    end do
    expected(1, 2) = 5; expected(2, 3) = -7; expected(1, 3) = -35.0_dp/2
    expected = exp((0.3_dp, 2.0_dp))*expected
    call matrix_exponential(a, e, info)
    call check('Jordan block', info, relative(e, expected), 1.0e-13_dp)

    ! Stiff similarity with eigenvalues -1e6, -1 + 3i, 2i: S diag(e^l) S^-1,
    ! S unit upper triangular so that S^-1 is exact.
    s = reshape([(1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
                 (2.0_dp, 1.0_dp), (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
                 (-1.0_dp, 0.5_dp), (3.0_dp, 0.0_dp), (1.0_dp, 0.0_dp)], [3, 3])
    sinv = 0
    do i = 1, 3
        sinv(i, i) = 1
    end do
    sinv(1, 2) = -s(1, 2); sinv(2, 3) = -s(2, 3)
    sinv(1, 3) = s(1, 2)*s(2, 3) - s(1, 3)
    lam = [(-1.0e6_dp, 0.0_dp), (-1.0_dp, 3.0_dp), (0.0_dp, 2.0_dp)]
    a = 0; expected = 0
    do i = 1, 3
        a(i, i) = lam(i); expected(i, i) = exp(lam(i))
    end do
    a = matmul(s, matmul(a, sinv))
    expected = matmul(s, matmul(expected, sinv))
    call check('S^-1 S exact', 0, relative(matmul(s, sinv), identity3()), 1.0e-15_dp)
    call matrix_exponential(a, e, info)
    ! Twenty squarings of a non-normal matrix amplify rounding (overscaling).
    call check('stiff similarity', info, relative(e, expected), 1.0e-9_dp)

    ! Semigroup: exp(A) exp(-A) = I for a dense complex matrix of norm ~ 30.
    do j = 1, 4
        do i = 1, 4
            c4(i, j) = cmplx(3*sin(1.7_dp*i + 0.3_dp*j), 2*cos(0.9_dp*i*j), dp)
        end do
    end do
    call matrix_exponential(c4, e4, info)
    call matrix_exponential(-c4, f4, info)
    c4 = matmul(e4, f4)
    do i = 1, 4
        c4(i, i) = c4(i, i) - 1
    end do
    call check('exp(A) exp(-A) = I', info, maxval(abs(c4)), 1.0e-9_dp)

    ! Real: rotation by a large angle needs squaring.
    theta = 50
    r = reshape([0.0_dp, theta, -theta, 0.0_dp], [2, 2])
    call matrix_exponential(r, er, info)
    call check('real rotation', info, maxval(abs(er - reshape([cos(theta), sin(theta), &
        -sin(theta), cos(theta)], [2, 2]))), 1.0e-12_dp)

    ! Real and complex paths agree on a real matrix.
    rr = reshape([-2.0_dp, 1.0_dp, 0.5_dp, 4.0_dp, -1.0_dp, 0.0_dp, 0.3_dp, 7.0_dp, &
                  -9.0_dp], [3, 3])
    call matrix_exponential(rr, err, info)
    call matrix_exponential(cmplx(rr, 0, dp), e, info)
    call check('real equals complex', info, &
               maxval(abs(err - e%re))/maxval(abs(err)) + maxval(abs(e%im)), 1.0e-14_dp)

    ! Failure status.
    a = 0
    a(2, 2) = cmplx(ieee_value(0.0_dp, ieee_quiet_nan), 0, dp)
    call matrix_exponential(a, e, info)
    call expect_failure('NaN input', info)
    a = 0
    a(1, 1) = 1.0e300_dp; a(2, 2) = (0.0_dp, 1.0e300_dp)
    call matrix_exponential(a*1.0e10_dp, e, info)
    call expect_failure('infinite input', info)
    call matrix_exponential(c4(:, 1:3), e4(:, 1:3), info)
    call expect_failure('non-square input', info)

    if (failures > 0) error stop 'test_fortnum_matrix_exponential failed'
    print '(a)', 'test_fortnum_matrix_exponential passed'
contains
    function identity3() result(id)
        complex(dp) :: id(3, 3)
        integer :: k
        id = 0
        do k = 1, 3
            id(k, k) = 1
        end do
    end function

    real(dp) function relative(x, y)
        complex(dp), intent(in) :: x(:, :), y(:, :)
        relative = maxval(abs(x - y))/maxval(abs(y))
    end function

    subroutine check(name, status, value, tolerance)
        character(*), intent(in) :: name
        integer, intent(in) :: status
        real(dp), intent(in) :: value, tolerance
        if (status == 0 .and. value <= tolerance) then
            write (*, '(a,a,es10.3)') 'ok   ', name//': ', value
        else
            write (*, '(a,a,es10.3,a,i0)') 'FAIL ', name//': ', value, ' info=', status
            failures = failures + 1
        end if
    end subroutine

    subroutine expect_failure(name, status)
        character(*), intent(in) :: name
        integer, intent(in) :: status
        if (status /= 0) then
            print '(a)', 'ok   '//name//' rejected'
        else
            print '(a)', 'FAIL '//name//' accepted'
            failures = failures + 1
        end if
    end subroutine
end program test_fortnum_matrix_exponential
