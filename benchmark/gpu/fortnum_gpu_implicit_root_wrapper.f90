module fortnum_gpu_implicit_root_wrapper
    use fortnum_kinds, only: dp
    use fortnum_generated_implicit_root_residual, only: &
        fortnum_implicit_root_residual_kernel
    use fortnum_roots, only: root_scalar_implicit_jvp_core
    implicit none
    private

    public :: implicit_root_jvp_batch

contains

    subroutine implicit_root_jvp_batch( &
            n, parameters, directions, roots, tangents, reliable)
        integer, intent(in) :: n
        real(dp), intent(in) :: parameters(n), directions(n)
        real(dp), intent(out) :: roots(n), tangents(n)
        logical, intent(out) :: reliable(n)
        real(dp) :: x, residual, f_x, f_p_tp
        integer :: i, iteration

        !$acc parallel loop present(parameters, directions, roots, tangents, reliable)
        !$omp target teams distribute parallel do &
        !$omp& map(to: parameters, directions) &
        !$omp& map(from: roots, tangents, reliable)
        do i = 1, n
            x = 1.0_dp
            do iteration = 1, 8
                x = 0.5_dp*(x + parameters(i)/x)
            end do
            call fortnum_implicit_root_residual_kernel( &
                x, parameters(i), directions(i), residual, f_x, f_p_tp)
            call root_scalar_implicit_jvp_core( &
                f_x, f_p_tp, 1.0e-14_dp, tangents(i), reliable(i))
            roots(i) = x
        end do
    end subroutine implicit_root_jvp_batch

end module fortnum_gpu_implicit_root_wrapper
