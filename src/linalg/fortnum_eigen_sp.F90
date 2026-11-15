module fortnum_eigen_sp
    use fortnum_eigen_types
    use fortnum_kinds, only: rk => sp
    use fortnum_eigen_dense_sp
#define GEMM eigen_sp_gemm
#define GEMV eigen_sp_gemv
#define GEES eigen_sp_gees
#define TREXC eigen_sp_trexc
#define TREVC eigen_sp_trevc
#define HEEV eigen_sp_heev
    use fortnum_eigen_radius_sp_kernel, only: bound_radius => fortnum_eigen_radius_sp
    use fortnum_eigen_temple_sp_kernel, only: bound_temple => fortnum_eigen_temple_sp
#include "fortnum_eigen_core.inc"
end module
