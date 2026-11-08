program test_fortnum_fourier
    ! Behavioral tests for the odd-grid Fourier collocation derivative:
    ! agreement with the DFT-mode definition, exact structure, and
    ! differentiation of resolved trigonometric data.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_fourier, only: fourier_derivative
    implicit none

    call test_matches_dft_definition()
    call test_structure()
    call test_trig_data()
    write (*, "(a)") "PASS: fortnum_fourier behavioral tests"
    stop 0

contains

    ! The matrix equals F^-1 diag(i m) F over the modes m = -(n-1)/2..(n-1)/2,
    ! an independent definition of the same derivative operator.
    subroutine test_matches_dft_definition()
        integer, parameter :: lengths(5) = [3, 5, 9, 17, 33]
        real(dp), parameter :: pi = acos(-1.0_dp)
        integer :: n, nt, j, k, q, m
        real(dp), allocatable :: d(:, :), ref(:, :)
        real(dp) :: err
        complex(dp), parameter :: iu = (0.0_dp, 1.0_dp)

        do n = 1, size(lengths)
            nt = lengths(n)
            allocate (d(nt, nt), ref(nt, nt))
            do j = 1, nt
                do k = 1, nt
                    ref(j, k) = 0
                    do q = 1, nt
                        m = merge(q - 1 - nt, q - 1, q > (nt + 1)/2)
                        ref(j, k) = ref(j, k) + iu*m &
                            *exp(iu*2*pi*(q - 1)*(j - k)/real(nt, dp))/nt
                    end do
                end do
            end do
            call fourier_derivative(nt, d)
            err = maxval(abs(d - ref))
            if (err > 1.0e-13_dp) error stop 'matrix disagrees with DFT definition'
            deallocate (d, ref)
        end do
    end subroutine test_matches_dft_definition

    ! Antisymmetry, zero diagonal and zero row sums hold exactly or at
    ! rounding; they distinguish the matrix from nearby wrong formulas.
    subroutine test_structure()
        integer, parameter :: n = 15
        real(dp) :: d(n, n), rowsum, diag
        integer :: i
        call fourier_derivative(n, d)
        if (maxval(abs(d + transpose(d))) > 0) error stop 'matrix is not antisymmetric'
        diag = 0
        rowsum = 0
        do i = 1, n
            diag = max(diag, abs(d(i, i)))
            rowsum = max(rowsum, abs(sum(d(i, :))))
        end do
        if (diag > 0) error stop 'diagonal is not zero'
        if (rowsum > 1.0e-14_dp) error stop 'rows do not annihilate constants'
    end subroutine test_structure

    ! Well-resolved cos and sin data differentiate to -k sin and k cos at
    ! the nodes, to rounding governed by the smallest singular value.
    subroutine test_trig_data()
        integer, parameter :: n = 33
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: d(n, n), theta(n), err
        integer :: k, i
        theta = 2*pi*[(real(i - 1, dp), i=1, n)]/n
        call fourier_derivative(n, d)
        do k = 1, 12
            err = maxval(abs(matmul(d, cos(k*theta)) + k*sin(k*theta)))
            if (err > 1.0e-13_dp*k) error stop 'cos data differentiates wrong'
            err = maxval(abs(matmul(d, sin(k*theta)) - k*cos(k*theta)))
            if (err > 1.0e-13_dp*k) error stop 'sin data differentiates wrong'
        end do
    end subroutine test_trig_data

end program test_fortnum_fourier
