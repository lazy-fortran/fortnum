pure subroutine fortnum_dot_sin_jvp_v(n_dir, n, a, a_d, b, b_d, s, s_d)
    implicit none
    integer, intent(in) :: n_dir
    integer, intent(in) :: n
    real(8), intent(in), dimension(n) :: a
    real(8), intent(in), dimension(n_dir, n) :: a_d
    real(8), intent(in), dimension(n) :: b
    real(8), intent(in), dimension(n_dir, n) :: b_d
    real(8), intent(out) :: s
    real(8), intent(out), dimension(n_dir) :: s_d
    integer :: i

    s_d(:) = 0.0d0
    s = 0.0d0
    do i = 1, n
        s_d(:) = s_d(:) + (a_d(:, i) * sin(b(i)) + a(i) * (cos(b(i)) * b_d(:, i)))
        s = s + a(i) * sin(b(i))
    end do
end subroutine fortnum_dot_sin_jvp_v

