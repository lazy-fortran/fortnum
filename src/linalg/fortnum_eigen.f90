module fortnum_eigen
    use fortnum_eigen_types
    use fortnum_eigen_sp, only: arnoldi_sp => krylov_schur, &
                          lanczos_sp => lanczos, real_arnoldi_sp => real_krylov_schur, &
                           real_lanczos_sp => real_lanczos, bounds_sp => eigen_bounds, &
                                inertia_sp => inertia_count, &
                                eigen_real_operator_sp => real_operator, &
                                eigen_complex_operator_sp => complex_operator, &
              eigen_real_solve_sp => real_solve, eigen_complex_solve_sp => complex_solve
    use fortnum_eigen_dp, only: arnoldi_dp => krylov_schur, &
                          lanczos_dp => lanczos, real_arnoldi_dp => real_krylov_schur, &
                           real_lanczos_dp => real_lanczos, bounds_dp => eigen_bounds, &
                                inertia_dp => inertia_count, &
                                eigen_real_operator_dp => real_operator, &
                                eigen_complex_operator_dp => complex_operator, &
              eigen_real_solve_dp => real_solve, eigen_complex_solve_dp => complex_solve
    implicit none
    private
    public :: krylov_schur, lanczos, eigen_bounds, inertia_count
    public :: eigen_options_t, eigen_stats_t
    public :: EIGEN_OK, EIGEN_NOT_CONVERGED, EIGEN_INVALID, EIGEN_OPERATOR_FAILED
    public :: EIGEN_METRIC_FAILED, EIGEN_DENSE_FAILED
    public :: EIGEN_SMALLEST, EIGEN_LARGEST, EIGEN_NEAREST, EIGEN_LARGEST_MAGNITUDE
    public :: eigen_real_operator_sp, eigen_complex_operator_sp
    public :: eigen_real_operator_dp, eigen_complex_operator_dp
    public :: eigen_real_solve_sp, eigen_complex_solve_sp
    public :: eigen_real_solve_dp, eigen_complex_solve_dp
    interface krylov_schur
        module procedure arnoldi_sp, arnoldi_dp, real_arnoldi_sp, real_arnoldi_dp
    end interface
    interface lanczos
        module procedure lanczos_sp, lanczos_dp, real_lanczos_sp, real_lanczos_dp
    end interface
    interface eigen_bounds
        module procedure bounds_sp, bounds_dp
    end interface
    interface inertia_count
        module procedure inertia_sp, inertia_dp
    end interface
end module
