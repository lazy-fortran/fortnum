!> Fourier collocation differentiation on odd periodic grids.
!>
!> For n odd equispaced nodes theta_j = 2 pi (j - 1)/n, the derivative of
!> the trigonometric interpolant through the nodal values is the matrix
!>   D(i, j) = 0.5 (-1)^(i - j)/sin(pi (i - j)/n),  i /= j,   D(i, i) = 0,
!> Trefethen, Spectral Methods in MATLAB, ch. 3; the modes run symmetrically
!> from -(n - 1)/2 to (n - 1)/2, so the sum of every row is zero and the
!> matrix is exactly antisymmetric. Even n has no Nyquist-free symmetric
!> mode set; its matrix has a different closed form and stays unsupported,
!> and callers must reject even n before calling. The routine is linear in
!> the nodal values; derivatives along data follow directly from applying
!> it to the differentiated values.
module fortnum_fourier
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private
    public :: fourier_derivative

    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp

contains

    !> First-derivative collocation matrix d of the n odd periodic nodes.
    pure subroutine fourier_derivative(n, d)
        integer, intent(in) :: n
        real(dp), intent(out) :: d(n, n)
        integer :: i, j
        d = 0
        do j = 1, n
            do i = 1, n
                if (i /= j) d(i, j) = 0.5_dp*(-1)**(i - j)/sin(pi*(i - j)/n)
            end do
        end do
    end subroutine fourier_derivative

end module fortnum_fourier
