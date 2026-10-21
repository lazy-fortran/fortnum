module fortnum_eigen_sp
    use fortnum_kinds, only: rk => sp
#define GEMM cgemm
#define GEMV cgemv
#define GEES cgees
#define TREXC ctrexc
#define TREVC ctrevc
#define HEEV cheev
    use fortnum_eigen_radius_sp_kernel, only: bound_radius => fortnum_eigen_radius_sp
    use fortnum_eigen_temple_sp_kernel, only: bound_temple => fortnum_eigen_temple_sp
#include "fortnum_eigen_core.inc"
end module
