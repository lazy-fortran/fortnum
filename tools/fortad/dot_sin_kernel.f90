! Testbed kernel for fortad, in the shape fortnum's own kernels take: a
! nonlinear reduction over two input arrays with an explicit intent(out)
! result. It is deliberately plain Fortran with no annotations, because that
! is the whole point of a source-transformation tool.
subroutine fortnum_dot_sin(n, a, b, s)
    integer, intent(in) :: n
    real(8), intent(in) :: a(n)
    real(8), intent(in) :: b(n)
    real(8), intent(out) :: s
    integer :: i
    s = 0.0d0
    do i = 1, n
        s = s + a(i)*sin(b(i))
    end do
end subroutine fortnum_dot_sin
