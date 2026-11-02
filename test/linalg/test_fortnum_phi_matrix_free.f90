module phi_stencil_operator
    use fortnum_kinds, only: dp
    use fortnum_krylov, only: real_linear_operator_t
    implicit none
    type, extends(real_linear_operator_t) :: stencil_t
    contains
        procedure :: apply => laplacian
    end type
contains
    subroutine laplacian(self, input, output)
        class(stencil_t), intent(inout) :: self
        real(dp), intent(in) :: input(:)
        real(dp), intent(out) :: output(:)
        integer :: i, n
        associate (unused => self)
        end associate
        n = size(input)
        output = -2.0_dp*input
        if (n < 2) return
        output(1) = output(1) + input(2)
        output(n) = output(n) + input(n - 1)
        !$omp parallel do simd schedule(static)
        do i = 2, n - 1
            output(i) = input(i - 1) - 2.0_dp*input(i) + input(i + 1)
        end do
        !$omp end parallel do simd
    end subroutine
end module

program test_fortnum_phi_matrix_free
    use fortnum_kinds, only: dp
    use fortnum_phi, only: phi_action, phi_options_t, phi_report_t, PHI_OK
    use phi_stencil_operator, only: stencil_t
    use omp_lib, only: omp_get_wtime, omp_get_max_threads
    implicit none
    integer, parameter :: n = 262144
    type(stencil_t) :: operator
    type(phi_options_t) :: options
    type(phi_report_t) :: report
    real(dp), allocatable :: v(:, :), expected(:), actual(:), low(:), high(:)
    real(dp) :: eigenvalue, weight, scalar, term, factorial, error, started
    integer :: i, k, mode, j, info
    allocate (v(n, 0:3), expected(n), actual(n), low(n), high(n))
    do i = 1, n
        low(i) = sin(acos(-1.0_dp)*i/(n + 1))
        high(i) = sin(acos(-1.0_dp)*((n - 1)/3)*i/(n + 1))
    end do
    do k = 0, 3
        v(:, k) = real(k + 1, dp)*low + real(4 - k, dp)*high
    end do
    expected = 0.0_dp
    do mode = 1, 2
        i = 1
        if (mode == 2) i = (n - 1)/3
        eigenvalue = -10.0_dp*sin(acos(-1.0_dp)*i/(2*(n + 1)))**2
        weight = 0.0_dp
        do k = 0, 3
            factorial = 1.0_dp
            do j = 2, k
                factorial = factorial*j
            end do
            scalar = exp(eigenvalue)
            if (abs(eigenvalue) > 1.0_dp) then
                factorial = 1.0_dp
                do j = 1, k
                    scalar = (scalar - 1.0_dp/factorial)/eigenvalue
                    factorial = factorial*j
                end do
            else
                term = 1.0_dp/factorial
                scalar = term
                do j = 1, 100
                    term = term*eigenvalue/(j + k)
                    scalar = scalar + term
                    if (abs(term) < 1.0e-17_dp*abs(scalar)) exit
                end do
            end if
            if (mode == 1) weight = weight + (k + 1)*scalar
            if (mode == 2) weight = weight + (4 - k)*scalar
        end do
        if (mode == 1) expected = expected + weight*low
        if (mode == 2) expected = expected + weight*high
    end do
    options%rtol = 1.0e-13_dp
    options%max_dimension = 32
    started = omp_get_wtime()
    call phi_action(operator, 2.5_dp, v, actual, info, options, report)
    started = omp_get_wtime() - started
    if (info /= PHI_OK) error stop 'matrix-free action failed'
    error = norm2(actual - expected)/norm2(expected)
    write (*, '(a,i0,a,i0,a,es10.2,a,i0,a,i0,a,f9.5,a)') &
        'n=', n, ' threads=', omp_get_max_threads(), ' relative error=', error, &
        ' m=', report%max_dimension, ' substeps=', report%substeps, &
        ' wall time=', started, ' s'
    if (error > 1.0e-10_dp) error stop 'matrix-free spectral oracle disagrees'
end program
