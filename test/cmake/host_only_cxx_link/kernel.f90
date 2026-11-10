module host_only_kernel
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
contains
    pure function sin_cos(x) result(y)
        !$omp declare target
        real(dp), intent(in) :: x
        real(dp) :: y
        y = sin(x)*cos(x)
    end function sin_cos
end module host_only_kernel

program host_only_cxx_link
    use host_only_kernel, only: dp, sin_cos
    implicit none
    real(dp), parameter :: x = 0.3_dp
    if (abs(sin_cos(x) - 0.5_dp*sin(2*x)) > 1.0e-15_dp) error stop 1
end program host_only_cxx_link
