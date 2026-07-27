module test_gpu_implicit_root_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_gpu_implicit_root_wrapper, only: implicit_root_jvp_batch
    implicit none
    private

    public :: run_implicit_root_test

contains

    subroutine run_implicit_root_test()
        integer, parameter :: n = 4096
        real(dp), parameter :: tolerance = 8.0e-14_dp
        real(dp) :: parameters(n), directions(n), roots(n), tangents(n)
        real(dp) :: expected_root, expected_tangent, scale
        logical :: reliable(n)
        integer :: i

        do i = 1, n
            parameters(i) = 0.5_dp + 1.5_dp* &
                real(mod(17*i, 4093), dp)/4092.0_dp
            directions(i) = -0.8_dp + 1.6_dp* &
                real(mod(23*i, 4091), dp)/4090.0_dp
        end do

        !$acc data copyin(parameters, directions) &
        !$acc& create(roots, tangents, reliable)
        !$omp target data map(to: parameters, directions) &
        !$omp& map(alloc: roots, tangents, reliable)
        call implicit_root_jvp_batch( &
            n, parameters, directions, roots, tangents, reliable)
        !$omp target update from(roots, tangents, reliable)
        !$omp end target data
        !$acc update self(roots, tangents, reliable)
        !$acc end data

        do i = 1, n
            expected_root = sqrt(parameters(i))
            expected_tangent = directions(i)/(2.0_dp*expected_root)
            scale = max(1.0_dp, abs(expected_root), abs(expected_tangent))
            if (.not. reliable(i)) then
                error stop "implicit GPU root unexpectedly unreliable"
            end if
            if (abs(roots(i) - expected_root) > tolerance*scale) then
                error stop "implicit GPU root disagrees with closed form"
            end if
            if (abs(tangents(i) - expected_tangent) > tolerance*scale) then
                error stop "implicit GPU JVP disagrees with closed form"
            end if
        end do
    end subroutine run_implicit_root_test

end module test_gpu_implicit_root_support
