module fortnum_fixed_point
    use fortnum_kinds, only: dp
    use fortnum_multiroot, only: multiroot_grad, multiroot_vjp
    use fortnum_status, only: fortnum_status_t
    implicit none
    private

    public :: fixed_point_jvp, fixed_point_vjp

contains

    ! Analytical implicit JVP for x = G(x,p):
    ! (I - G_x) dx = (G_p tp).
    subroutine fixed_point_jvp(map_x, map_p, tp, dx, status)
        real(dp), intent(in) :: map_x(:, :), map_p(:, :), tp(:)
        real(dp), intent(out) :: dx(:)
        type(fortnum_status_t), intent(out) :: status

        real(dp) :: residual_x(size(map_x, 1), size(map_x, 2))
        real(dp) :: residual_p_tangent(size(map_x, 1))
        integer :: i, j

        residual_x = -map_x
        do i = 1, size(map_x, 1)
            residual_x(i, i) = residual_x(i, i) + 1.0_dp
        end do
        residual_p_tangent = 0.0_dp
        do j = 1, size(map_p, 2)
            do i = 1, size(map_p, 1)
                residual_p_tangent(i) = residual_p_tangent(i) - map_p(i, j)*tp(j)
            end do
        end do
        call multiroot_grad(residual_x, residual_p_tangent, dx, status)
    end subroutine fixed_point_jvp

    ! Analytical implicit VJP for x = G(x,p):
    ! (I - G_x)^T lambda = u, p_bar = G_p^T lambda.
    subroutine fixed_point_vjp(map_x, map_p, u, jtu, status)
        real(dp), intent(in) :: map_x(:, :), map_p(:, :), u(:)
        real(dp), intent(out) :: jtu(:)
        type(fortnum_status_t), intent(out) :: status

        real(dp) :: residual_x(size(map_x, 1), size(map_x, 2))
        real(dp) :: residual_p(size(map_p, 1), size(map_p, 2))
        integer :: i

        residual_x = -map_x
        do i = 1, size(map_x, 1)
            residual_x(i, i) = residual_x(i, i) + 1.0_dp
        end do
        residual_p = -map_p
        call multiroot_vjp(residual_x, residual_p, u, jtu, status)
    end subroutine fixed_point_vjp

end module fortnum_fixed_point
