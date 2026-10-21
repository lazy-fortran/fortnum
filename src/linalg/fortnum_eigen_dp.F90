module fortnum_eigen_dp
    use fortnum_kinds, only: rk => dp
#define GEMM zgemm
#define GEMV zgemv
#define GEES zgees
#define TREXC ztrexc
#define TREVC ztrevc
#define HEEV zheev
    use fortnum_eigen_radius_dp_kernel, only: bound_radius => fortnum_eigen_radius_dp
    use fortnum_eigen_temple_dp_kernel, only: bound_temple => fortnum_eigen_temple_dp
#include "fortnum_eigen_core.inc"
end module

