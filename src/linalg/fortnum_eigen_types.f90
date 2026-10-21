module fortnum_eigen_types
    implicit none
    private
    integer, parameter, public :: EIGEN_OK = 0, EIGEN_NOT_CONVERGED = 1, &
             EIGEN_INVALID = -1, EIGEN_OPERATOR_FAILED = -2, EIGEN_METRIC_FAILED = -3, &
                                  EIGEN_DENSE_FAILED = -4
    integer, parameter, public :: EIGEN_SMALLEST = 1, EIGEN_LARGEST = 2, &
                                  EIGEN_NEAREST = 3, EIGEN_LARGEST_MAGNITUDE = 4
    type, public :: eigen_options_t
        integer :: max_basis = 32, max_restarts = 100, retain = 0
        integer :: which = EIGEN_SMALLEST
    end type
    type, public :: eigen_stats_t
        integer :: restarts = 0, applications = 0, converged = 0, locked = 0
    end type
end module
