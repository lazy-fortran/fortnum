module fortnum_eigen_dp
    use fortnum_eigen_types
    use fortnum_kinds, only: rk => dp
    use fortnum_blas, only: fortnum_zgemm
    use fortnum_blas12, only: fortnum_zgemv
    use fortnum_lapack_gees, only: fortnum_zgees, fortnum_ztrexc, fortnum_ztrevc
    use fortnum_lapack_syev, only: fortnum_zheev
#define GEMM fortnum_zgemm
#define GEMV fortnum_zgemv
#define GEES fortnum_zgees
#define TREXC fortnum_ztrexc
#define TREVC fortnum_ztrevc
#define HEEV fortnum_zheev
    use fortnum_eigen_radius_dp_kernel, only: bound_radius => fortnum_eigen_radius_dp
    use fortnum_eigen_temple_dp_kernel, only: bound_temple => fortnum_eigen_temple_dp
#include "fortnum_eigen_core.inc"
end module

