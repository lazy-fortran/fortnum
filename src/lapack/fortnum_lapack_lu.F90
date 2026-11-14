module fortnum_lapack_lu
    ! LU factorization and solves on fortnum's BLAS: DGETRF/ZGETRF (recursive),
    ! DGETRF2/ZGETRF2,
    ! DGETRS/ZGETRS, DGESV/ZGESV, DGETRI/ZGETRI, DLASWP/ZLASWP and the band
    ! versions DGBTRF/ZGBTRF, DGBTRS/ZGBTRS, DGBSV/ZGBSV, with the arguments
    ! and info codes of reference LAPACK (prefix fortnum_).
    use fortnum_kinds, only: dp
    use fortnum_blas, only: fortnum_dgemm, fortnum_zgemm
    use fortnum_blas3, only: fortnum_dtrsm, fortnum_ztrsm
    implicit none
    private

    public :: fortnum_dgetrf, fortnum_dgetrf2, fortnum_dgetf2, fortnum_dlaswp, fortnum_dgetrs, fortnum_dgesv
    public :: fortnum_dgetri, fortnum_dgbtrf, fortnum_dgbtrs, fortnum_dgbsv
    public :: fortnum_zgetrf, fortnum_zgetrf2, fortnum_zgetf2, fortnum_zlaswp, fortnum_zgetrs, fortnum_zgesv
    public :: fortnum_zgetri, fortnum_zgbtrf, fortnum_zgbtrs, fortnum_zgbsv


contains

#define SCALAR real(dp)
#define ABS1(x) abs(x)
#define CONJ(x) (x)
#define ONE 1.0_dp
#define ZERO 0.0_dp
#define GEMM fortnum_dgemm
#define TRSM fortnum_dtrsm
#define GETRF fortnum_dgetrf
#define GETF2 fortnum_dgetf2
#define GETRF2 fortnum_dgetrf2
#define LASWP fortnum_dlaswp
#define GETRS fortnum_dgetrs
#define GESV fortnum_dgesv
#define GETRI fortnum_dgetri
#define GBTRF fortnum_dgbtrf
#define GBTRS fortnum_dgbtrs
#define GBSV fortnum_dgbsv
#include "fortnum_lapack_lu.inc"
#undef SCALAR
#undef ABS1
#undef CONJ
#undef ONE
#undef ZERO
#undef GEMM
#undef TRSM
#undef GETRF
#undef GETF2
#undef GETRF2
#undef LASWP
#undef GETRS
#undef GESV
#undef GETRI
#undef GBTRF
#undef GBTRS
#undef GBSV
#define SCALAR complex(dp)
#define ABS1(x) (abs(real(x, dp)) + abs(aimag(x)))
#define CONJ(x) conjg(x)
#define ONE (1.0_dp, 0.0_dp)
#define ZERO (0.0_dp, 0.0_dp)
#define GEMM fortnum_zgemm
#define TRSM fortnum_ztrsm
#define GETRF fortnum_zgetrf
#define GETF2 fortnum_zgetf2
#define GETRF2 fortnum_zgetrf2
#define LASWP fortnum_zlaswp
#define GETRS fortnum_zgetrs
#define GESV fortnum_zgesv
#define GETRI fortnum_zgetri
#define GBTRF fortnum_zgbtrf
#define GBTRS fortnum_zgbtrs
#define GBSV fortnum_zgbsv
#include "fortnum_lapack_lu.inc"

end module fortnum_lapack_lu
