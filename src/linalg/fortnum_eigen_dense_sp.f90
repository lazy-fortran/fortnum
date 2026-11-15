module fortnum_eigen_dense_sp
    ! The dense steps of the single precision Krylov-Schur solver
    ! (fortnum_eigen_sp) with the arguments of CGEMM, CGEMV, CGEES, CTREXC,
    ! CTREVC and CHEEV. They are computed in double precision by fortnum's
    ! ZGEMM, ZGEMV, ZGEES, ZTREXC, ZTREVC and ZHEEV and rounded to single,
    ! which keeps the backward errors at single precision level. eigen_sp_gees
    ! supports sort 'N' only, the solver's use.
    use fortnum_kinds, only: sp, dp
    use fortnum_blas, only: fortnum_zgemm
    use fortnum_blas12, only: fortnum_zgemv
    use fortnum_lapack_gees, only: fortnum_zgees, fortnum_ztrexc, fortnum_ztrevc
    use fortnum_lapack_syev, only: fortnum_zheev
    implicit none
    private

    public :: eigen_sp_gemm, eigen_sp_gemv, eigen_sp_gees, eigen_sp_trexc
    public :: eigen_sp_trevc, eigen_sp_heev

contains

    subroutine eigen_sp_gemm(ta, tb, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc)
        character, intent(in) :: ta, tb
        integer, intent(in) :: m, n, k, lda, ldb, ldc
        complex(sp), intent(in) :: alpha, beta, a(lda, *), b(ldb, *)
        complex(sp), intent(inout) :: c(ldc, *)
        complex(dp), allocatable :: ad(:, :), bd(:, :), cd(:, :)
        integer :: ka, kb
        ka = k
        if (.not. any(ta == ['N', 'n'])) ka = m
        kb = n
        if (.not. any(tb == ['N', 'n'])) kb = k
        ad = cmplx(a(1:lda, 1:ka), kind=dp)
        bd = cmplx(b(1:ldb, 1:kb), kind=dp)
        cd = cmplx(c(1:m, 1:n), kind=dp)
        call fortnum_zgemm(ta, tb, m, n, k, cmplx(alpha, kind=dp), ad, lda, bd, ldb, &
            cmplx(beta, kind=dp), cd, max(1, m))
        c(1:m, 1:n) = cmplx(cd, kind=sp)
    end subroutine eigen_sp_gemm

    subroutine eigen_sp_gemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        character, intent(in) :: trans
        integer, intent(in) :: m, n, lda, incx, incy
        complex(sp), intent(in) :: alpha, beta, a(lda, *), x(*)
        complex(sp), intent(inout) :: y(*)
        complex(dp), allocatable :: ad(:, :), xd(:), yd(:)
        integer :: lenx, leny
        lenx = n
        leny = m
        if (.not. any(trans == ['N', 'n'])) then
            lenx = m
            leny = n
        end if
        ad = cmplx(a(1:lda, 1:n), kind=dp)
        xd = cmplx(x(1:1 + (lenx - 1)*abs(incx)), kind=dp)
        yd = cmplx(y(1:1 + (leny - 1)*abs(incy)), kind=dp)
        call fortnum_zgemv(trans, m, n, cmplx(alpha, kind=dp), ad, lda, xd, incx, &
            cmplx(beta, kind=dp), yd, incy)
        y(1:1 + (leny - 1)*abs(incy)) = cmplx(yd, kind=sp)
    end subroutine eigen_sp_gemv

    subroutine eigen_sp_gees(jobvs, sort, select, n, a, lda, sdim, w, vs, ldvs, &
            work, lwork, rwork, bwork, info)
        character, intent(in) :: jobvs, sort
        interface
            logical function select(w)
                import :: sp
                complex(sp), intent(in) :: w
            end function select
        end interface
        integer, intent(in) :: n, lda, ldvs, lwork
        complex(sp), intent(inout) :: a(lda, *)
        integer, intent(out) :: sdim, info
        complex(sp), intent(out) :: w(*), vs(ldvs, *), work(*)
        real(sp), intent(out) :: rwork(*)
        logical, intent(out) :: bwork(*)
        complex(dp), allocatable :: ad(:, :), vd(:, :), wd(:), wk(:)
        real(dp), allocatable :: rw(:)
        if (.not. any(sort == ['N', 'n'])) then
            info = -2
            return
        end if
        ad = cmplx(a(1:lda, 1:n), kind=dp)
        allocate (vd(ldvs, max(1, n)), wd(max(1, n)), wk(max(1, 2*n)), rw(max(1, n)))
        call fortnum_zgees(jobvs, sort, no_select, n, ad, lda, sdim, wd, vd, ldvs, wk, &
            size(wk), rw, bwork, info)
        a(1:lda, 1:n) = cmplx(ad, kind=sp)
        w(1:n) = cmplx(wd(1:n), kind=sp)
        if (any(jobvs == ['V', 'v'])) vs(1:ldvs, 1:n) = cmplx(vd(:, 1:n), kind=sp)
        work(1) = max(1, 2*n)
        rwork(1) = 0
    end subroutine eigen_sp_gees

    logical function no_select(w)
        complex(dp), intent(in) :: w
        no_select = .false.
    end function no_select

    subroutine eigen_sp_trexc(compq, n, t, ldt, q, ldq, ifst, ilst, info)
        character, intent(in) :: compq
        integer, intent(in) :: n, ldt, ldq, ifst, ilst
        integer, intent(out) :: info
        complex(sp), intent(inout) :: t(ldt, *), q(ldq, *)
        complex(dp), allocatable :: td(:, :), qd(:, :)
        logical :: wantq
        wantq = any(compq == ['V', 'v'])
        td = cmplx(t(1:ldt, 1:n), kind=dp)
        if (wantq) then
            qd = cmplx(q(1:ldq, 1:n), kind=dp)
        else
            allocate (qd(ldq, 1))
        end if
        call fortnum_ztrexc(compq, n, td, ldt, qd, ldq, ifst, ilst, info)
        t(1:ldt, 1:n) = cmplx(td, kind=sp)
        if (wantq) q(1:ldq, 1:n) = cmplx(qd, kind=sp)
    end subroutine eigen_sp_trexc

    subroutine eigen_sp_trevc(side, how, select, n, t, ldt, vl, ldvl, vr, ldvr, &
            mm, m, work, rwork, info)
        character, intent(in) :: side, how
        logical, intent(inout) :: select(*)
        integer, intent(in) :: n, ldt, ldvl, ldvr, mm
        integer, intent(out) :: m, info
        complex(sp), intent(in) :: t(ldt, *)
        complex(sp), intent(inout) :: vl(ldvl, *), vr(ldvr, *)
        complex(sp), intent(out) :: work(*)
        real(sp), intent(out) :: rwork(*)
        complex(dp), allocatable :: td(:, :), vld(:, :), vrd(:, :), wk(:)
        real(dp), allocatable :: rw(:)
        logical :: leftv, rightv
        leftv = any(side == ['L', 'l', 'B', 'b'])
        rightv = any(side == ['R', 'r', 'B', 'b'])
        td = cmplx(t(1:ldt, 1:n), kind=dp)
        if (leftv) then
            vld = cmplx(vl(1:ldvl, 1:mm), kind=dp)
        else
            allocate (vld(ldvl, 1))
        end if
        if (rightv) then
            vrd = cmplx(vr(1:ldvr, 1:mm), kind=dp)
        else
            allocate (vrd(ldvr, 1))
        end if
        allocate (wk(max(1, 2*n)), rw(max(1, n)))
        call fortnum_ztrevc(side, how, select, n, td, ldt, vld, ldvl, vrd, ldvr, mm, m, &
            wk, rw, info)
        if (leftv) vl(1:ldvl, 1:mm) = cmplx(vld, kind=sp)
        if (rightv) vr(1:ldvr, 1:mm) = cmplx(vrd, kind=sp)
        work(1) = 0
        rwork(1) = 0
    end subroutine eigen_sp_trevc

    subroutine eigen_sp_heev(job, uplo, n, a, lda, w, work, lwork, rwork, info)
        character, intent(in) :: job, uplo
        integer, intent(in) :: n, lda, lwork
        integer, intent(out) :: info
        complex(sp), intent(inout) :: a(lda, *)
        real(sp), intent(out) :: w(*), rwork(*)
        complex(sp), intent(out) :: work(*)
        complex(dp), allocatable :: ad(:, :), wk(:)
        real(dp), allocatable :: wd(:), rw(:)
        complex(dp) :: wq(1)
        real(dp) :: rq(1)
        ad = cmplx(a(1:lda, 1:n), kind=dp)
        allocate (wd(max(1, n)), rw(max(1, 3*n - 2)))
        call fortnum_zheev(job, uplo, n, ad, lda, wd, wq, -1, rq, info)
        allocate (wk(max(1, int(real(wq(1), dp)))))
        call fortnum_zheev(job, uplo, n, ad, lda, wd, wk, size(wk), rw, info)
        a(1:lda, 1:n) = cmplx(ad, kind=sp)
        w(1:n) = real(wd(1:n), sp)
        work(1) = 0
        rwork(1) = 0
    end subroutine eigen_sp_heev

end module fortnum_eigen_dense_sp
