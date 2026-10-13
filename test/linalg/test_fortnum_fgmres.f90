program test_fortnum_fgmres
    !! FGMRES: exact-inverse, block-Jacobi, variable (inner GMRES) and
    !! stagnation behaviour against dense LU oracles.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_krylov, only: complex_fgmres_operator, real_fgmres_operator, &
        real_gmres_operator, KRYLOV_OK, KRYLOV_STAGNATION
    use fortnum_linalg, only: dense_solve, LINALG_OK
    implicit none

    integer, parameter :: n = 120, nb = 10
    real(dp) :: a(n, n), b(n), x(n), x_plain(n), x_ref(n), res
    real(dp) :: spd(6, 6), spd_inv(6, 6), rb(6), rx(6)
    complex(dp) :: cm(5, 5), cinv(5, 5), cb(5), cx(5), cref(5)
    real(dp) :: rot(2, 2)
    integer :: info, it_plain, it_jac, it_var, i, j, it
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp)

    ! (1a) SPD, exact inverse: one iteration
    spd = 0.0_dp
    do i = 1, 6
        spd(i, i) = 4.0_dp + i
        if (i < 6) spd(i, i + 1) = -1.0_dp
        if (i < 6) spd(i + 1, i) = -1.0_dp
    end do
    call invert(6, spd, spd_inv)
    rb = [(real(i, dp), i=1, 6)]
    rx = 0.0_dp
    call real_fgmres_operator(spd_mv, spd_pc, rb, rx, 1.0e-12_dp, 20, 10, &
        info, it, res)
    call require(info == KRYLOV_OK .and. it == 1, "SPD exact inverse: 1 iteration")
    call require(maxval(abs(matmul(spd, rx) - rb)) < 1.0e-11_dp, "SPD solution")

    ! (1b) nonsymmetric complex, exact inverse: one iteration
    cm = czero
    do i = 1, 5
        cm(i, i) = cmplx(2.5_dp + 0.2_dp*i, 0.1_dp, dp)
        if (i < 5) cm(i, i + 1) = cmplx(-0.7_dp, 0.3_dp, dp)
        if (i < 5) cm(i + 1, i) = cmplx(0.15_dp, -0.05_dp, dp)
    end do
    cm(1, 5) = cmplx(0.4_dp, -0.2_dp, dp)
    cb = [(cmplx(real(i, dp), (-1.0_dp)**i*0.25_dp, dp), i=1, 5)]
    call dense_solve(cm, cb, cref, info)
    call require(info == LINALG_OK, "complex dense oracle")
    do j = 1, 5
        cx = czero
        cx(j) = (1.0_dp, 0.0_dp)
        call dense_solve(cm, cx, cinv(:, j), info)
    end do
    cx = czero
    call complex_fgmres_operator(c_mv, c_pc, cb, cx, 1.0e-12_dp, 20, 10, &
        info, it, res)
    call require(info == KRYLOV_OK .and. it == 1, "complex exact inverse: 1 iteration")
    call require(maxval(abs(cx - cref)) < 1.0e-11_dp, "complex solution vs LU")
    ! complex, no useful preconditioner (identity): converges, restarted
    cx = czero
    call complex_fgmres_operator(c_mv, c_id, cb, cx, 1.0e-12_dp, 40, 3, &
        info, it, res)
    call require(info == KRYLOV_OK .and. it > 1, "complex restarted identity")
    call require(maxval(abs(cx - cref)) < 1.0e-10_dp, "complex restarted solution")

    ! (2) convection-diffusion, block Jacobi vs none
    call build_convdiff()
    b = 1.0_dp
    call dense_solve(a, b, x_ref, info)
    call require(info == LINALG_OK, "convdiff dense oracle")
    x_plain = 0.0_dp
    call real_fgmres_operator(a_mv, id_pc, b, x_plain, 1.0e-10_dp, n, n, &
        info, it_plain, res)
    call require(info == KRYLOV_OK, "unpreconditioned FGMRES converges")
    x = 0.0_dp
    call real_fgmres_operator(a_mv, jacobi_pc, b, x, 1.0e-10_dp, n, n, &
        info, it_jac, res)
    call require(info == KRYLOV_OK, "block-Jacobi FGMRES converges")
    call require(it_jac < it_plain, "block Jacobi reduces iterations")
    call require(maxval(abs(x - x_ref)) < 1.0e-7_dp*maxval(abs(x_ref)), &
        "block-Jacobi solution matches LU")
    call require(maxval(abs(x_plain - x_ref)) < 1.0e-7_dp*maxval(abs(x_ref)), &
        "unpreconditioned solution matches LU")
    print '(a,i0,a,i0)', "iterations plain/jacobi: ", it_plain, "/", it_jac

    ! (3) variable preconditioner: 3 steps of inner GMRES (nonlinear in input)
    x = 0.0_dp
    call real_fgmres_operator(a_mv, inner_pc, b, x, 1.0e-10_dp, 400, 20, &
        info, it_var, res)
    call require(info == KRYLOV_OK, "variable preconditioner converges")
    call require(maxval(abs(x - x_ref)) < 1.0e-7_dp*maxval(abs(x_ref)), &
        "variable-preconditioner solution matches LU")

    ! (4) stagnation: GMRES(1) on a rotation makes no progress
    rot = reshape([0.0_dp, 1.0_dp, -1.0_dp, 0.0_dp], [2, 2])
    rb(1:2) = [1.0_dp, 0.0_dp]
    rx(1:2) = 0.0_dp
    call real_fgmres_operator(rot_mv, rot_id, rb(1:2), rx(1:2), 1.0e-10_dp, &
        50, 1, info, it, res)
    call require(info == KRYLOV_STAGNATION, "stagnation is reported")
    call require(it < 50 .and. res > 0.5_dp, "stagnation exits early, residual honest")

    print '(a)', "PASS: fortnum_fgmres"

contains

    subroutine invert(m, mat, inv)
        integer, intent(in) :: m
        real(dp), intent(in) :: mat(m, m)
        real(dp), intent(out) :: inv(m, m)
        real(dp) :: e(m)
        integer :: k, st
        do k = 1, m
            e = 0.0_dp
            e(k) = 1.0_dp
            call dense_solve(mat, e, inv(:, k), st)
        end do
    end subroutine invert

    subroutine build_convdiff()
        real(dp), parameter :: eps = 0.02_dp, conv = 1.0_dp
        real(dp) :: h
        integer :: k
        h = 1.0_dp/(n + 1)
        a = 0.0_dp
        do k = 1, n
            a(k, k) = 2.0_dp*eps/h**2 + 0.5_dp
            if (k > 1) a(k, k - 1) = -eps/h**2 - conv/(2.0_dp*h)
            if (k < n) a(k, k + 1) = -eps/h**2 + conv/(2.0_dp*h)
        end do
        a = a*h
    end subroutine build_convdiff

    subroutine spd_mv(u, v)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: v(:)
        v = matmul(spd, u)
    end subroutine spd_mv

    subroutine spd_pc(u, v)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: v(:)
        v = matmul(spd_inv, u)
    end subroutine spd_pc

    subroutine c_mv(u, v)
        complex(dp), intent(in) :: u(:)
        complex(dp), intent(out) :: v(:)
        v = matmul(cm, u)
    end subroutine c_mv

    subroutine c_pc(u, v)
        complex(dp), intent(in) :: u(:)
        complex(dp), intent(out) :: v(:)
        v = matmul(cinv, u)
    end subroutine c_pc

    subroutine c_id(u, v)
        complex(dp), intent(in) :: u(:)
        complex(dp), intent(out) :: v(:)
        v = u
    end subroutine c_id

    subroutine a_mv(u, v)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: v(:)
        v = matmul(a, u)
    end subroutine a_mv

    subroutine id_pc(u, v)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: v(:)
        v = u
    end subroutine id_pc

    subroutine jacobi_pc(u, v)
        !! Block-diagonal solve with nb x nb diagonal blocks of a.
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: v(:)
        integer :: k, lo, hi, st
        do k = 1, n/nb
            lo = (k - 1)*nb + 1
            hi = k*nb
            call dense_solve(a(lo:hi, lo:hi), u(lo:hi), v(lo:hi), st)
        end do
    end subroutine jacobi_pc

    subroutine inner_pc(u, v)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: v(:)
        integer :: st, its
        real(dp) :: r
        v = 0.0_dp
        call real_gmres_operator(a_mv, u, v, 1.0e-14_dp, 3, 3, st, its, r)
    end subroutine inner_pc

    subroutine rot_mv(u, v)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: v(:)
        v = matmul(rot, u)
    end subroutine rot_mv

    subroutine rot_id(u, v)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: v(:)
        v = u
    end subroutine rot_id

    subroutine require(condition, description)
        logical, intent(in) :: condition
        character(*), intent(in) :: description
        if (.not. condition) then
            write (error_unit, '(a)') "FAIL: "//description
            error stop 1
        end if
    end subroutine require

end program test_fortnum_fgmres
