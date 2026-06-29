program test_fortnum_linalg
    ! Behavioral tests for fortnum_linalg: closed-form 2x2/3x3 determinant and
    ! inverse round-trips A*inv(A) = I, the near-singular jacobian_ok3 predicate
    ! (true / false / NaN), the in-place partial-pivot LU solve against known
    ! systems at several sizes including a pivot-required case, an inv3 cross
    ! check, and the singular-status contract for inv2/inv3/lu_solve.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_linalg, only: det2, det3, inv2, inv3, jacobian_ok3, lu_solve, &
        LINALG_OK, LINALG_SINGULAR, LINALG_MAX_N
    implicit none

    integer :: nfail
    nfail = 0

    call test_det(nfail)
    call test_inv2_roundtrip(nfail)
    call test_inv3_roundtrip(nfail)
    call test_jacobian_ok3(nfail)
    call test_lu_known(nfail)
    call test_lu_pivot(nfail)
    call test_lu_vs_inv3(nfail)
    call test_singular(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine check(label, got, expected, atol, nfail)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected, atol
        integer, intent(inout) :: nfail
        if (abs(got - expected) > atol) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,es12.4)") "FAIL [", label, "] err=", &
                abs(got - expected)
        end if
    end subroutine check

    subroutine check_true(label, cond, nfail)
        character(*), intent(in) :: label
        logical, intent(in) :: cond
        integer, intent(inout) :: nfail
        if (.not. cond) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a)") "FAIL [", label, "] expected .true."
        end if
    end subroutine check_true

    subroutine check_int(label, got, expected, nfail)
        character(*), intent(in) :: label
        integer, intent(in) :: got, expected
        integer, intent(inout) :: nfail
        if (got /= expected) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,i0,a,i0)") "FAIL [", label, "] got=", &
                got, " expected=", expected
        end if
    end subroutine check_int

    subroutine test_det(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: a2(2, 2), a3(3, 3)
        a2 = reshape([2.0_dp, 1.0_dp, 3.0_dp, 4.0_dp], [2, 2]) ! det = 8-3 = 5
        call check("det2", det2(a2), 5.0_dp, 1.0e-13_dp, nfail)
        ! Upper triangular: det = product of diagonal = 2*3*4 = 24.
        a3 = reshape([2.0_dp, 0.0_dp, 0.0_dp, 5.0_dp, 3.0_dp, 0.0_dp, &
            7.0_dp, 9.0_dp, 4.0_dp], [3, 3])
        call check("det3_triangular", det3(a3), 24.0_dp, 1.0e-12_dp, nfail)
    end subroutine test_det

    subroutine test_inv2_roundtrip(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: a(2, 2), ainv(2, 2), prod(2, 2)
        integer :: info
        a = reshape([4.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], [2, 2])
        call inv2(a, ainv, info)
        call check_int("inv2_ok", info, LINALG_OK, nfail)
        prod = matmul(a, ainv)
        call check("inv2_rt", maxval(abs(prod - eye(2))), 0.0_dp, 1.0e-12_dp, nfail)
    end subroutine test_inv2_roundtrip

    subroutine test_inv3_roundtrip(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: a(3, 3), ainv(3, 3), prod(3, 3)
        integer :: info
        a = reshape([4.0_dp, 1.0_dp, 0.5_dp, 0.2_dp, 3.0_dp, 1.1_dp, &
            0.7_dp, 0.3_dp, 5.0_dp], [3, 3])
        call inv3(a, ainv, info)
        call check_int("inv3_ok", info, LINALG_OK, nfail)
        prod = matmul(a, ainv)
        call check("inv3_rt", maxval(abs(prod - eye(3))), 0.0_dp, 1.0e-12_dp, nfail)
    end subroutine test_inv3_roundtrip

    subroutine test_jacobian_ok3(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: a(3, 3), bad(3, 3), nanm(3, 3)
        real(dp) :: nan
        a = reshape([4.0_dp, 1.0_dp, 0.5_dp, 0.2_dp, 3.0_dp, 1.1_dp, &
            0.7_dp, 0.3_dp, 5.0_dp], [3, 3])
        call check_true("jac_ok_true", jacobian_ok3(a), nfail)
        ! Row scaled by 1e-10 pushes |det| below the 1e-8*scale^3 floor.
        bad = a
        bad(2, :) = a(2, :)*1.0e-10_dp
        call check_true("jac_ok_false",.not. jacobian_ok3(bad), nfail)
        nan = 0.0_dp
        nan = nan/nan
        nanm = a
        nanm(1, 1) = nan
        call check_true("jac_ok_nan",.not. jacobian_ok3(nanm), nfail)
    end subroutine test_jacobian_ok3

    subroutine test_lu_known(nfail)
        integer, intent(inout) :: nfail
        call lu_case(2, nfail)
        call lu_case(3, nfail)
        call lu_case(5, nfail)
        call lu_case(LINALG_MAX_N, nfail)
    end subroutine test_lu_known

    ! Solve A x = b with A = I + small structured perturbation, x_true ramp.
    subroutine lu_case(n, nfail)
        integer, intent(in) :: n
        integer, intent(inout) :: nfail
        real(dp) :: a(n, n), b(n), xtrue(n)
        integer :: i, j, info
        character(len=24) :: lbl
        do i = 1, n
            xtrue(i) = real(i, dp)
            do j = 1, n
                if (i == j) then
                    a(i, j) = 2.0_dp + real(i, dp)
                else
                    a(i, j) = 0.3_dp/real(i + j, dp)
                end if
            end do
        end do
        b = matmul(a, xtrue)
        call lu_solve(n, a, b, info)
        write (lbl, "(a,i0)") "lu_n", n
        call check_int(trim(lbl), info, LINALG_OK, nfail)
        call check(trim(lbl)//"_x", maxval(abs(b - xtrue)), 0.0_dp, 1.0e-10_dp, nfail)
    end subroutine lu_case

    ! A = [[0,1],[1,0]] requires a pivot swap; x_true = [3, 5].
    subroutine test_lu_pivot(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: a(2, 2), b(2), xtrue(2)
        integer :: info
        a = reshape([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
        xtrue = [3.0_dp, 5.0_dp]
        b = matmul(a, xtrue)
        call lu_solve(2, a, b, info)
        call check_int("lu_pivot_ok", info, LINALG_OK, nfail)
        call check("lu_pivot_x", maxval(abs(b - xtrue)), 0.0_dp, 1.0e-12_dp, nfail)
    end subroutine test_lu_pivot

    ! lu_solve(3,A,b) must match matmul(inv3(A), b).
    subroutine test_lu_vs_inv3(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: a(3, 3), acopy(3, 3), ainv(3, 3), b(3), bcopy(3), xinv(3)
        integer :: info
        a = reshape([4.0_dp, 1.0_dp, 0.5_dp, 0.2_dp, 3.0_dp, 1.1_dp, &
            0.7_dp, 0.3_dp, 5.0_dp], [3, 3])
        b = [1.0_dp, -2.0_dp, 0.5_dp]
        acopy = a
        bcopy = b
        call inv3(a, ainv, info)
        xinv = matmul(ainv, b)
        call lu_solve(3, acopy, bcopy, info)
        call check("lu_vs_inv3", maxval(abs(bcopy - xinv)), 0.0_dp, 1.0e-11_dp, nfail)
    end subroutine test_lu_vs_inv3

    subroutine test_singular(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: a2(2, 2), ainv2(2, 2), a3(3, 3), ainv3(3, 3)
        real(dp) :: as(3, 3), bs(3)
        integer :: info
        ! Two identical rows -> rank deficient.
        a2 = reshape([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], [2, 2])
        call inv2(a2, ainv2, info)
        call check_int("inv2_sing", info, LINALG_SINGULAR, nfail)
        call check("inv2_sing_zero", maxval(abs(ainv2)), 0.0_dp, 0.0_dp, nfail)
        a3 = reshape([1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 6.0_dp, 3.0_dp, &
            0.5_dp, 1.0_dp, 0.5_dp], [3, 3]) ! columns dependent
        call inv3(a3, ainv3, info)
        call check_int("inv3_sing", info, LINALG_SINGULAR, nfail)
        call check("inv3_sing_zero", maxval(abs(ainv3)), 0.0_dp, 0.0_dp, nfail)
        ! lu_solve on a singular A: column 2 = 2*column 1.
        as = reshape([1.0_dp, 0.0_dp, 2.0_dp, 2.0_dp, 0.0_dp, 4.0_dp, &
            0.0_dp, 1.0_dp, 1.0_dp], [3, 3])
        bs = [1.0_dp, 1.0_dp, 1.0_dp]
        call lu_solve(3, as, bs, info)
        call check_true("lu_sing", info > 0, nfail)
    end subroutine test_singular

    pure function eye(n) result(m)
        integer, intent(in) :: n
        real(dp) :: m(n, n)
        integer :: i, j
        do j = 1, n
            do i = 1, n
                if (i == j) then
                    m(i, j) = 1.0_dp
                else
                    m(i, j) = 0.0_dp
                end if
            end do
        end do
    end function eye

end program test_fortnum_linalg
