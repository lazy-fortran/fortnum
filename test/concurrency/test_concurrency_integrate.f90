program test_concurrency_integrate
    ! M5.1 concurrency: quadrature (gauss_legendre) and adaptive integration
    ! (QAG/QAGS/QAGP/QAGIU). The adaptive drivers carry caller-owned
    ! integrate_workspace_t and integrate_epstab_t; thread safety rests on each
    ! thread owning its own workspace, epstab, and result. This test integrates
    ! a family of problems serially, then with one private workspace per thread,
    ! and asserts bit-for-bit equality of every reported value and error.
    !
    ! PRIMAL concurrency only; derivative-product (frozen-subdivision) integrate
    ! concurrency lands in M6 (#40).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t
    use fortnum_quadrature, only: gauss_legendre
    use fortnum_integrate, only: integrate_qag, integrate_qags, &
        integrate_qagp, integrate_qagiu, &
        integrate_workspace_t, integrate_epstab_t, &
        integrate_result_t
    implicit none

    integer, parameter :: nprob = 256
    integer :: nfail, j
    real(dp) :: widths(nprob)
    real(dp) :: ser_v(nprob), ser_e(nprob), par_v(nprob), par_e(nprob)

    nfail = 0
    do j = 1, nprob
        widths(j) = 1.0e-3_dp + 0.4_dp * real(j - 1, dp) / real(nprob, dp)
    end do

    call check_gauss_legendre()

    call serial_qag()
    call parallel_qag()
    call check_exact(ser_v, par_v, "QAG value parallel == serial")
    call check_exact(ser_e, par_e, "QAG abserr parallel == serial")

    call serial_qags()
    call parallel_qags()
    call check_exact(ser_v, par_v, "QAGS value parallel == serial")
    call check_exact(ser_e, par_e, "QAGS abserr parallel == serial")

    call serial_qagp()
    call parallel_qagp()
    call check_exact(ser_v, par_v, "QAGP value parallel == serial")

    call serial_qagiu()
    call parallel_qagiu()
    call check_exact(ser_v, par_v, "QAGIU value parallel == serial")

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "test_concurrency_integrate: all tests passed"

contains

    function f_lorentz(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx, w
        w = 1.0e-2_dp
        if (present(ctx)) then
            select type (ctx)
                type is (real(dp))
                w = ctx
            end select
        end if
        fx = 1.0_dp/(1.0_dp + ((x - 0.5_dp)/w)**2)
    end function f_lorentz

    function f_recip(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = 1.0_dp/(1.0_dp + x*x) ! integral 0..inf = pi/2
    end function f_recip

    function f_sqrtsing(x, ctx) result(fx)
        real(dp), intent(in) :: x
        class(*), intent(in), optional :: ctx
        real(dp) :: fx
        fx = 1.0_dp/sqrt(abs(x - 0.5_dp)) ! integrable singularity at 0.5
    end function f_sqrtsing

    ! gauss_legendre fills caller-owned x,w; each thread owns its own arrays.
    subroutine check_gauss_legendre()
        integer, parameter :: m = 24
        real(dp) :: serx(m, nprob), serw(m, nprob)
        real(dp) :: parx(m, nprob), parw(m, nprob)
        do j = 1, nprob
            call gauss_legendre(m, serx(:, j), serw(:, j))
        end do
        !$omp parallel do default(shared) private(j) schedule(static)
        do j = 1, nprob
            call gauss_legendre(m, parx(:, j), parw(:, j))
        end do
        !$omp end parallel do
        if (any(serx /= parx) .or. any(serw /= parw)) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: gauss_legendre parallel == serial"
        end if
    end subroutine check_gauss_legendre

    subroutine serial_qag()
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        do j = 1, nprob
            call integrate_qag(f_lorentz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-9_dp, &
                ws, res, st, ctx=widths(j))
            ser_v(j) = res%value
            ser_e(j) = res%abserr
        end do
    end subroutine serial_qag

    subroutine parallel_qag()
        type(integrate_workspace_t) :: ws
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        !$omp parallel do default(shared) private(ws, res, st, j) &
        !$omp   schedule(static)
        do j = 1, nprob
            call integrate_qag(f_lorentz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-9_dp, &
                ws, res, st, ctx=widths(j))
            par_v(j) = res%value
            par_e(j) = res%abserr
        end do
        !$omp end parallel do
    end subroutine parallel_qag

    subroutine serial_qags()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        do j = 1, nprob
            call integrate_qags(f_lorentz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-9_dp, &
                ws, eps, res, st, ctx=widths(j))
            ser_v(j) = res%value
            ser_e(j) = res%abserr
        end do
    end subroutine serial_qags

    subroutine parallel_qags()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        !$omp parallel do default(shared) private(ws, eps, res, st, j) &
        !$omp   schedule(static)
        do j = 1, nprob
            call integrate_qags(f_lorentz, 0.0_dp, 1.0_dp, 0.0_dp, 1.0e-9_dp, &
                ws, eps, res, st, ctx=widths(j))
            par_v(j) = res%value
            par_e(j) = res%abserr
        end do
        !$omp end parallel do
    end subroutine parallel_qags

    subroutine serial_qagp()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        real(dp) :: pts(1)
        pts = [0.5_dp]
        do j = 1, nprob
            call integrate_qagp(f_sqrtsing, 0.0_dp, 1.0_dp, pts, 0.0_dp, &
                1.0e-8_dp, ws, eps, res, st)
            ser_v(j) = res%value
        end do
    end subroutine serial_qagp

    subroutine parallel_qagp()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        real(dp) :: pts(1)
        !$omp parallel do default(shared) private(ws, eps, res, st, pts, j) &
        !$omp   schedule(static)
        do j = 1, nprob
            pts = [0.5_dp]
            call integrate_qagp(f_sqrtsing, 0.0_dp, 1.0_dp, pts, 0.0_dp, &
                1.0e-8_dp, ws, eps, res, st)
            par_v(j) = res%value
        end do
        !$omp end parallel do
    end subroutine parallel_qagp

    subroutine serial_qagiu()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        do j = 1, nprob
            call integrate_qagiu(f_recip, 0.0_dp, 1, 0.0_dp, 1.0e-9_dp, &
                ws, eps, res, st)
            ser_v(j) = res%value
        end do
    end subroutine serial_qagiu

    subroutine parallel_qagiu()
        type(integrate_workspace_t) :: ws
        type(integrate_epstab_t)    :: eps
        type(integrate_result_t)    :: res
        type(fortnum_status_t)      :: st
        !$omp parallel do default(shared) private(ws, eps, res, st, j) &
        !$omp   schedule(static)
        do j = 1, nprob
            call integrate_qagiu(f_recip, 0.0_dp, 1, 0.0_dp, 1.0e-9_dp, &
                ws, eps, res, st)
            par_v(j) = res%value
        end do
        !$omp end parallel do
    end subroutine parallel_qagiu

    subroutine check_exact(ref, got, name)
        real(dp),     intent(in) :: ref(:), got(:)
        character(*), intent(in) :: name
        integer :: k
        do k = 1, size(ref)
            if (ref(k) /= got(k)) then
                nfail = nfail + 1
                write (error_unit, "(a,a,a,i0)") "FAIL: ", name, " at index ", k
                return
            end if
        end do
    end subroutine check_exact

end program test_concurrency_integrate
