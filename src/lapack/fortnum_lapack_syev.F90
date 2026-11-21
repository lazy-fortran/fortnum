module fortnum_lapack_syev
    ! Symmetric and Hermitian eigenproblems on fortnum's BLAS, with reference
    ! LAPACK arguments and info codes (prefix fortnum_): the reductions
    ! DSYTRD/ZHETRD, DORGTR/ZUNGTR, DSYGST/ZHEGST, the tridiagonal solvers
    ! DSTEQR/ZSTEQR, DSTERF, DSTEV, and the drivers DSYEV, ZHEEV, DSYEVD,
    ! DSYGV, ZHEGV, DSBEVD, DSBEVX. The band drivers expand the band and use
    ! the dense path; a band reduction is planned (issue 167).
    use fortnum_kinds, only: dp
    use fortnum_blas3, only: fortnum_dtrsm, fortnum_ztrsm, fortnum_dtrmm, fortnum_ztrmm
    use fortnum_lapack_chol, only: fortnum_dpotrf, fortnum_zpotrf
    use fortnum_lapack_qr, only: fortnum_dlarfg, fortnum_zlarfg, fortnum_dorgqr, fortnum_zungqr
    implicit none
    private

    public :: fortnum_dsytrd, fortnum_dorgtr, fortnum_dsygst, fortnum_dsteqr, fortnum_dsterf
    public :: fortnum_dsyev, fortnum_dsyevd, fortnum_dsygv, fortnum_dstev, fortnum_dsbevd, fortnum_dsbevx
    public :: fortnum_zhetrd, fortnum_zungtr, fortnum_zhegst, fortnum_zsteqr, fortnum_zheev, fortnum_zhegv

contains

#define SCALAR real(dp)
#define ZERO 0.0_dp
#define ONE 1.0_dp
#define CONJ(x) (x)
#define TRANS_C 'T'
#define LARFG fortnum_dlarfg
#define TRSM fortnum_dtrsm
#define TRMM fortnum_dtrmm
#define POTRF fortnum_dpotrf
#define SYTRD fortnum_dsytrd
#define ORGTR fortnum_dorgtr
#define ORGQR fortnum_dorgqr
#define STEQR dsteqr_ql
#define SYEV_CORE dsyev_core
#define SYGST fortnum_dsygst
#define SYGV_CORE dsygv_core
#include "fortnum_lapack_syev.inc"
#undef SCALAR
#undef ZERO
#undef ONE
#undef CONJ
#undef TRANS_C
#undef LARFG
#undef TRSM
#undef TRMM
#undef POTRF
#undef SYTRD
#undef ORGTR
#undef ORGQR
#undef STEQR
#undef SYEV_CORE
#undef SYGST
#undef SYGV_CORE
#define SCALAR complex(dp)
#define ZERO (0.0_dp, 0.0_dp)
#define ONE (1.0_dp, 0.0_dp)
#define CONJ(x) conjg(x)
#define TRANS_C 'C'
#define LARFG fortnum_zlarfg
#define TRSM fortnum_ztrsm
#define TRMM fortnum_ztrmm
#define POTRF fortnum_zpotrf
#define SYTRD fortnum_zhetrd
#define ORGTR fortnum_zungtr
#define ORGQR fortnum_zungqr
#define STEQR zsteqr_ql
#define SYEV_CORE zheev_core
#define SYGST fortnum_zhegst
#define SYGV_CORE zhegv_core
#include "fortnum_lapack_syev.inc"

    !> Tridiagonal eigenproblem (DSTEQR): compz N (values), V (z holds the
    !> orthogonal matrix of a reduction) or I (z starts as the identity).
    subroutine fortnum_dsteqr(compz, n, d, e, z, ldz, work, info)
        character, intent(in) :: compz
        integer, intent(in) :: n, ldz
        real(dp), intent(inout) :: d(*), e(*), z(ldz, *), work(*)
        integer, intent(out) :: info
        integer :: j
        info = 0
        if (index('NnVvIi', compz) == 0) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (ldz < 1 .or. (index('VvIi', compz) > 0 .and. ldz < max(1, n))) then
            info = -6
        end if
        if (info /= 0) return
        if (compz == 'I' .or. compz == 'i') then
            z(1:n, 1:n) = 0
            do j = 1, n
                z(j, j) = 1
            end do
        end if
        call dsteqr_ql(index('NnVv', compz) == 0 .or. compz == 'V' .or. compz == 'v', n, d, e, z, ldz, info)
    end subroutine fortnum_dsteqr

    subroutine fortnum_zsteqr(compz, n, d, e, z, ldz, work, info)
        character, intent(in) :: compz
        integer, intent(in) :: n, ldz
        real(dp), intent(inout) :: d(*), e(*), work(*)
        complex(dp), intent(inout) :: z(ldz, *)
        integer, intent(out) :: info
        integer :: j
        info = 0
        if (index('NnVvIi', compz) == 0) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (ldz < 1 .or. (index('VvIi', compz) > 0 .and. ldz < max(1, n))) then
            info = -6
        end if
        if (info /= 0) return
        if (compz == 'I' .or. compz == 'i') then
            z(1:n, 1:n) = 0
            do j = 1, n
                z(j, j) = 1
            end do
        end if
        call zsteqr_ql(index('VvIi', compz) > 0, n, d, e, z, ldz, info)
    end subroutine fortnum_zsteqr

    !> Eigenvalues of a symmetric tridiagonal matrix (DSTERF), ascending.
    subroutine fortnum_dsterf(n, d, e, info)
        integer, intent(in) :: n
        real(dp), intent(inout) :: d(*), e(*)
        integer, intent(out) :: info
        real(dp) :: z(1, 1)
        info = 0
        if (n < 0) then
            info = -1
            return
        end if
        call dsteqr_ql(.false., n, d, e, z, 1, info)
    end subroutine fortnum_dsterf

    !> Eigenvalues and optionally eigenvectors of a symmetric tridiagonal
    !> matrix (DSTEV).
    subroutine fortnum_dstev(jobz, n, d, e, z, ldz, work, info)
        character, intent(in) :: jobz
        integer, intent(in) :: n, ldz
        real(dp), intent(inout) :: d(*), e(*), z(ldz, *), work(*)
        integer, intent(out) :: info
        logical :: wantz
        integer :: j
        info = 0
        wantz = jobz == 'V' .or. jobz == 'v'
        if (.not. (wantz .or. jobz == 'N' .or. jobz == 'n')) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (ldz < 1 .or. (wantz .and. ldz < n)) then
            info = -6
        end if
        if (info /= 0 .or. n == 0) return
        if (wantz) then
            z(1:n, 1:n) = 0
            do j = 1, n
                z(j, j) = 1
            end do
        end if
        call dsteqr_ql(wantz, n, d, e, z, ldz, info)
    end subroutine fortnum_dstev

    logical function bad_job(jobz, uplo, n, lda, info)
        character, intent(in) :: jobz, uplo
        integer, intent(in) :: n, lda
        integer, intent(out) :: info
        info = 0
        if (index('NnVv', jobz) == 0) then
            info = -1
        else if (index('UuLl', uplo) == 0) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (lda < max(1, n)) then
            info = -5
        end if
        bad_job = info /= 0
    end function bad_job

    !> Eigenvalues and optionally eigenvectors of a symmetric matrix (DSYEV).
    subroutine fortnum_dsyev(jobz, uplo, n, a, lda, w, work, lwork, info)
        character, intent(in) :: jobz, uplo
        integer, intent(in) :: n, lda, lwork
        real(dp), intent(inout) :: a(lda, *), work(*)
        real(dp), intent(out) :: w(*)
        integer, intent(out) :: info
        if (bad_job(jobz, uplo, n, lda, info)) return
        if (lwork < max(1, 3*n - 1) .and. lwork /= -1) then
            info = -8
            return
        end if
        work(1) = max(1, 34*n)
        if (lwork == -1) return
        call dsyev_core(jobz, uplo, n, a, lda, w, info)
    end subroutine fortnum_dsyev

    !> DSYEVD: the results of DSYEV with LAPACK's divide-and-conquer
    !> workspace contract.
    subroutine fortnum_dsyevd(jobz, uplo, n, a, lda, w, work, lwork, iwork, liwork, info)
        character, intent(in) :: jobz, uplo
        integer, intent(in) :: n, lda, lwork, liwork
        real(dp), intent(inout) :: a(lda, *), work(*)
        real(dp), intent(out) :: w(*)
        integer, intent(inout) :: iwork(*)
        integer, intent(out) :: info
        integer :: lwmin, liwmin
        logical :: wantz
        if (bad_job(jobz, uplo, n, lda, info)) return
        wantz = jobz == 'V' .or. jobz == 'v'
        if (n <= 1) then
            lwmin = 1
            liwmin = 1
        else if (wantz) then
            lwmin = 1 + 6*n + 2*n**2
            liwmin = 3 + 5*n
        else
            lwmin = 2*n + 1
            liwmin = 1
        end if
        work(1) = lwmin
        iwork(1) = liwmin
        if (lwork < lwmin .and. lwork /= -1) then
            info = -8
        else if (liwork < liwmin .and. liwork /= -1) then
            info = -10
        end if
        if (info /= 0 .or. lwork == -1 .or. liwork == -1) return
        call dsyev_core(jobz, uplo, n, a, lda, w, info)
        work(1) = lwmin
        iwork(1) = liwmin
    end subroutine fortnum_dsyevd

    !> Eigenvalues and optionally eigenvectors of a Hermitian matrix (ZHEEV).
    subroutine fortnum_zheev(jobz, uplo, n, a, lda, w, work, lwork, rwork, info)
        character, intent(in) :: jobz, uplo
        integer, intent(in) :: n, lda, lwork
        complex(dp), intent(inout) :: a(lda, *), work(*)
        real(dp), intent(out) :: w(*)
        real(dp), intent(inout) :: rwork(*)
        integer, intent(out) :: info
        if (bad_job(jobz, uplo, n, lda, info)) return
        if (lwork < max(1, 2*n - 1) .and. lwork /= -1) then
            info = -8
            return
        end if
        work(1) = max(1, 33*n)
        if (lwork == -1) return
        call zheev_core(jobz, uplo, n, a, lda, w, info)
    end subroutine fortnum_zheev

    !> Generalized symmetric-definite eigenproblem (DSYGV).
    subroutine fortnum_dsygv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, info)
        integer, intent(in) :: itype, n, lda, ldb, lwork
        character, intent(in) :: jobz, uplo
        real(dp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
        real(dp), intent(out) :: w(*)
        integer, intent(out) :: info
        info = 0
        if (itype < 1 .or. itype > 3) then
            info = -1
        else if (index('NnVv', jobz) == 0) then
            info = -2
        else if (index('UuLl', uplo) == 0) then
            info = -3
        else if (n < 0) then
            info = -4
        else if (lda < max(1, n)) then
            info = -6
        else if (ldb < max(1, n)) then
            info = -8
        else if (lwork < max(1, 3*n - 1) .and. lwork /= -1) then
            info = -11
        end if
        if (info /= 0) return
        work(1) = max(1, 34*n)
        if (lwork == -1 .or. n == 0) return
        call dsygv_core(itype, jobz, uplo, n, a, lda, b, ldb, w, info)
    end subroutine fortnum_dsygv

    !> Generalized Hermitian-definite eigenproblem (ZHEGV).
    subroutine fortnum_zhegv(itype, jobz, uplo, n, a, lda, b, ldb, w, work, lwork, rwork, info)
        integer, intent(in) :: itype, n, lda, ldb, lwork
        character, intent(in) :: jobz, uplo
        complex(dp), intent(inout) :: a(lda, *), b(ldb, *), work(*)
        real(dp), intent(out) :: w(*)
        real(dp), intent(inout) :: rwork(*)
        integer, intent(out) :: info
        info = 0
        if (itype < 1 .or. itype > 3) then
            info = -1
        else if (index('NnVv', jobz) == 0) then
            info = -2
        else if (index('UuLl', uplo) == 0) then
            info = -3
        else if (n < 0) then
            info = -4
        else if (lda < max(1, n)) then
            info = -6
        else if (ldb < max(1, n)) then
            info = -8
        else if (lwork < max(1, 2*n - 1) .and. lwork /= -1) then
            info = -11
        end if
        if (info /= 0) return
        work(1) = max(1, 33*n)
        if (lwork == -1 .or. n == 0) return
        call zhegv_core(itype, jobz, uplo, n, a, lda, b, ldb, w, info)
    end subroutine fortnum_zhegv

    ! the symmetric band matrix ab (kd off-diagonals, uplo) as a dense matrix
    subroutine expand_band(uplo, n, kd, ab, ldab, a)
        character, intent(in) :: uplo
        integer, intent(in) :: n, kd, ldab
        real(dp), intent(in) :: ab(ldab, *)
        real(dp), intent(out) :: a(n, n)
        integer :: i, j
        a = 0
        do j = 1, n
            if (uplo == 'U' .or. uplo == 'u') then
                do i = max(1, j - kd), j
                    a(i, j) = ab(kd + 1 + i - j, j)
                    a(j, i) = a(i, j)
                end do
            else
                do i = j, min(n, j + kd)
                    a(i, j) = ab(1 + i - j, j)
                    a(j, i) = a(i, j)
                end do
            end if
        end do
    end subroutine expand_band

    !> Symmetric band eigenproblem (DSBEVD), LAPACK's workspace contract.
    subroutine fortnum_dsbevd(jobz, uplo, n, kd, ab, ldab, w, z, ldz, work, lwork, iwork, liwork, info)
        character, intent(in) :: jobz, uplo
        integer, intent(in) :: n, kd, ldab, ldz, lwork, liwork
        real(dp), intent(inout) :: ab(ldab, *), work(*)
        real(dp), intent(out) :: w(*), z(ldz, *)
        integer, intent(inout) :: iwork(*)
        integer, intent(out) :: info
        logical :: wantz
        integer :: lwmin, liwmin
        real(dp), allocatable :: a(:, :)
        info = 0
        wantz = jobz == 'V' .or. jobz == 'v'
        if (n <= 1) then
            lwmin = 1
            liwmin = 1
        else if (wantz) then
            lwmin = 1 + 5*n + 2*n**2
            liwmin = 3 + 5*n
        else
            lwmin = 2*n
            liwmin = 1
        end if
        work(1) = lwmin
        iwork(1) = liwmin
        if (.not. (wantz .or. jobz == 'N' .or. jobz == 'n')) then
            info = -1
        else if (index('UuLl', uplo) == 0) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (kd < 0) then
            info = -4
        else if (ldab < kd + 1) then
            info = -6
        else if (ldz < 1 .or. (wantz .and. ldz < n)) then
            info = -9
        else if (lwork < lwmin .and. lwork /= -1) then
            info = -11
        else if (liwork < liwmin .and. liwork /= -1) then
            info = -13
        end if
        if (info /= 0 .or. lwork == -1 .or. liwork == -1 .or. n == 0) return
        allocate (a(n, n))
        call expand_band(uplo, n, kd, ab, ldab, a)
        call dsyev_core(jobz, 'L', n, a, n, w, info)
        if (wantz) z(1:n, 1:n) = a
        work(1) = lwmin
        iwork(1) = liwmin
    end subroutine fortnum_dsbevd

    !> Selected eigenvalues (range A all, V in (vl, vu], I indices il..iu)
    !> and eigenvectors of a symmetric band matrix (DSBEVX); q is the identity
    !> of the dense path, ifail is zero.
    subroutine fortnum_dsbevx(jobz, range, uplo, n, kd, ab, ldab, q, ldq, vl, vu, il, iu, abstol, m, w, &
            z, ldz, work, iwork, ifail, info)
        character, intent(in) :: jobz, range, uplo
        integer, intent(in) :: n, kd, ldab, ldq, il, iu, ldz
        real(dp), intent(inout) :: ab(ldab, *), work(*)
        real(dp), intent(out) :: q(ldq, *), w(*), z(ldz, *)
        real(dp), intent(in) :: vl, vu, abstol
        integer, intent(out) :: m, ifail(*), info
        integer, intent(inout) :: iwork(*)
        logical :: wantz, alleig, valeig, indeig
        real(dp), allocatable :: a(:, :), ev(:)
        integer :: j, first, last
        info = 0
        m = 0
        wantz = jobz == 'V' .or. jobz == 'v'
        alleig = range == 'A' .or. range == 'a'
        valeig = range == 'V' .or. range == 'v'
        indeig = range == 'I' .or. range == 'i'
        if (.not. (wantz .or. jobz == 'N' .or. jobz == 'n')) then
            info = -1
        else if (.not. (alleig .or. valeig .or. indeig)) then
            info = -2
        else if (index('UuLl', uplo) == 0) then
            info = -3
        else if (n < 0) then
            info = -4
        else if (kd < 0) then
            info = -5
        else if (ldab < kd + 1) then
            info = -7
        else if (wantz .and. ldq < max(1, n)) then
            info = -9
        else if (valeig .and. n > 0 .and. vu <= vl) then
            info = -11
        else if (indeig .and. (il < 1 .or. il > max(1, n))) then
            info = -12
        else if (indeig .and. (iu < min(n, il) .or. iu > n)) then
            info = -13
        else if (ldz < 1 .or. (wantz .and. ldz < n)) then
            info = -18
        end if
        if (info /= 0 .or. n == 0) return
        allocate (a(n, n), ev(n))
        call expand_band(uplo, n, kd, ab, ldab, a)
        call dsyev_core(jobz, 'L', n, a, n, ev, info)
        if (info /= 0) return
        if (alleig) then
            first = 1
            last = n
        else if (indeig) then
            first = il
            last = iu
        else
            first = n + 1
            last = 0
            do j = 1, n
                if (ev(j) > vl .and. ev(j) <= vu) then
                    first = min(first, j)
                    last = j
                end if
            end do
        end if
        m = max(0, last - first + 1)
        if (m > 0) then
            w(1:m) = ev(first:last)
            if (wantz) z(1:n, 1:m) = a(:, first:last)
        end if
        if (wantz) then
            q(1:n, 1:n) = 0
            do j = 1, n
                q(j, j) = 1
            end do
            ifail(1:max(1, m)) = 0
        end if
    end subroutine fortnum_dsbevx


    ! Plane rotation with c f + s g = r and c g - s f = 0 (LAPACK DLARTG)
    pure subroutine givens(f, g, c, s, r)
        real(dp), intent(in) :: f, g
        real(dp), intent(out) :: c, s, r
        if (g == 0) then
            c = 1
            s = 0
            r = f
        else if (f == 0) then
            c = 0
            s = 1
            r = g
        else
            r = hypot(f, g)
            c = f/r
            s = g/r
        end if
    end subroutine givens

    ! Eigenvalues rt1 (larger in magnitude) and rt2 of [a b; b c] and the
    ! unit eigenvector (cs1, sn1) of rt1 (LAPACK DLAEV2)
    pure subroutine sym2x2(a, b, c, rt1, rt2, cs1, sn1)
        real(dp), intent(in) :: a, b, c
        real(dp), intent(out) :: rt1, rt2, cs1, sn1
        real(dp) :: sm, df, adf, tb, ab, acmx, acmn, rt, cs, acs, ct, tn
        integer :: sgn1, sgn2
        sm = a + c
        df = a - c
        adf = abs(df)
        tb = b + b
        ab = abs(tb)
        if (abs(a) > abs(c)) then
            acmx = a
            acmn = c
        else
            acmx = c
            acmn = a
        end if
        if (adf > ab) then
            rt = adf*sqrt(1 + (ab/adf)**2)
        else if (adf < ab) then
            rt = ab*sqrt(1 + (adf/ab)**2)
        else
            rt = ab*sqrt(2.0_dp)
        end if
        if (sm < 0) then
            rt1 = 0.5_dp*(sm - rt)
            sgn1 = -1
            rt2 = (acmx/rt1)*acmn - (b/rt1)*b
        else if (sm > 0) then
            rt1 = 0.5_dp*(sm + rt)
            sgn1 = 1
            rt2 = (acmx/rt1)*acmn - (b/rt1)*b
        else
            rt1 = 0.5_dp*rt
            rt2 = -0.5_dp*rt
            sgn1 = 1
        end if
        if (df >= 0) then
            cs = df + rt
            sgn2 = 1
        else
            cs = df - rt
            sgn2 = -1
        end if
        acs = abs(cs)
        if (acs > ab) then
            ct = -tb/cs
            sn1 = 1/sqrt(1 + ct*ct)
            cs1 = ct*sn1
        else if (ab == 0) then
            cs1 = 1
            sn1 = 0
        else
            tn = -cs/tb
            cs1 = 1/sqrt(1 + tn*tn)
            sn1 = tn*cs1
        end if
        if (sgn1 == sgn2) then
            tn = cs1
            cs1 = -sn1
            sn1 = tn
        end if
    end subroutine sym2x2

end module fortnum_lapack_syev
