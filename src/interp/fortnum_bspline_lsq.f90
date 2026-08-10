module fortnum_bspline_lsq
    !! Matrix-free least-squares fitting for the Fortnum B-spline workspace.
    !!
    !! The basis matrix is assembled once in caller-independent storage, while
    !! the CGLS iteration applies it and its transpose without forming normal
    !! equations. This keeps the fitting layer on top of the public
    !! fortnum_bspline representation instead of introducing another spline
    !! type.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_bspline, only: bspline_workspace_t, bspline_eval_basis
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none
    private

    public :: bspline_1d_lsq_cgls

contains

    ! Fit f_data at x_data in the least-squares sense using the B-spline basis
    ! in ws. The stopping tolerance is relative to max(||f_data||_2, 1), and
    ! the residual norm is returned even when the iteration limit is reached.
    subroutine bspline_1d_lsq_cgls(ws, x_data, f_data, coeff, max_iterations, &
            tolerance, status, iterations, residual_norm)
        type(bspline_workspace_t), intent(in) :: ws
        real(dp), intent(in) :: x_data(:), f_data(:)
        real(dp), intent(out) :: coeff(:)
        integer, intent(in), optional :: max_iterations
        real(dp), intent(in), optional :: tolerance
        type(fortnum_status_t), intent(out) :: status
        integer, intent(out), optional :: iterations
        real(dp), intent(out), optional :: residual_norm

        integer :: i, max_iter, n_data
        real(dp) :: tol, rhs_norm, target, gamma, gamma_new, alpha, beta
        real(dp) :: denom, current_norm
        real(dp), allocatable :: basis(:, :), residual(:), q(:)
        real(dp), allocatable :: gradient(:), direction(:)
        type(fortnum_status_t) :: eval_status

        call status_set(status, FORTNUM_OK, "")
        if (present(iterations)) iterations = 0
        if (present(residual_norm)) residual_norm = 0.0_dp
        coeff = 0.0_dp

        if (.not. ws%knots_set) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "bspline_lsq: knots must be set before fitting")
            return
        end if
        if (size(x_data) == 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "bspline_lsq: x_data must not be empty")
            return
        end if
        if (size(f_data) /= size(x_data)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "bspline_lsq: x_data and f_data sizes differ")
            return
        end if
        if (size(coeff) /= ws%ncoef) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "bspline_lsq: coeff size /= workspace ncoef")
            return
        end if

        max_iter = 200
        if (present(max_iterations)) max_iter = max_iterations
        if (max_iter < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "bspline_lsq: max_iterations must be positive")
            return
        end if

        tol = 1.0e-10_dp
        if (present(tolerance)) tol = tolerance
        if (tol <= 0.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "bspline_lsq: tolerance must be positive")
            return
        end if

        n_data = size(x_data)
        allocate (basis(ws%ncoef, n_data), residual(n_data), q(n_data))
        allocate (gradient(ws%ncoef), direction(ws%ncoef))

        do i = 1, n_data
            call bspline_eval_basis(ws, x_data(i), basis(:, i), eval_status)
            if (eval_status%code /= FORTNUM_OK) then
                status = eval_status
                return
            end if
        end do

        residual = f_data
        rhs_norm = vector_norm(f_data)
        current_norm = rhs_norm
        target = tol*max(rhs_norm, 1.0_dp)
        if (present(residual_norm)) residual_norm = current_norm
        if (current_norm <= target) return

        call apply_transpose(basis, residual, gradient)
        direction = gradient
        gamma = dot_product(gradient, gradient)
        if (gamma <= 0.0_dp) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                "bspline_lsq: zero initial gradient")
            return
        end if

        do i = 1, max_iter
            call apply_basis(basis, direction, q)
            denom = dot_product(q, q)
            if (denom <= 0.0_dp) then
                if (present(iterations)) iterations = i - 1
                if (present(residual_norm)) residual_norm = current_norm
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "bspline_lsq: singular search direction")
                return
            end if

            alpha = gamma/denom
            coeff = coeff + alpha*direction
            residual = residual - alpha*q
            current_norm = vector_norm(residual)
            if (present(iterations)) iterations = i
            if (present(residual_norm)) residual_norm = current_norm
            if (current_norm <= target) return

            call apply_transpose(basis, residual, gradient)
            gamma_new = dot_product(gradient, gradient)
            if (gamma_new <= 0.0_dp) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "bspline_lsq: normal residual stagnated")
                return
            end if
            beta = gamma_new/gamma
            direction = gradient + beta*direction
            gamma = gamma_new
        end do

        call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
            "bspline_lsq: iteration limit reached")
    end subroutine bspline_1d_lsq_cgls

    subroutine apply_basis(basis, coeff, values)
        real(dp), intent(in) :: basis(:, :), coeff(:)
        real(dp), intent(out) :: values(:)
        integer :: i, j

        values = 0.0_dp
        do i = 1, size(values)
            do j = 1, size(coeff)
                values(i) = values(i) + basis(j, i)*coeff(j)
            end do
        end do
    end subroutine apply_basis

    subroutine apply_transpose(basis, values, coeff)
        real(dp), intent(in) :: basis(:, :), values(:)
        real(dp), intent(out) :: coeff(:)
        integer :: i, j

        coeff = 0.0_dp
        do j = 1, size(coeff)
            do i = 1, size(values)
                coeff(j) = coeff(j) + basis(j, i)*values(i)
            end do
        end do
    end subroutine apply_transpose

    pure function vector_norm(values) result(result)
        real(dp), intent(in) :: values(:)
        real(dp) :: result

        result = sqrt(dot_product(values, values))
    end function vector_norm

end module fortnum_bspline_lsq
