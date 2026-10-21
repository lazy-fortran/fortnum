module fortnum_refinement
    !! Mixed-precision defect correction; the exact operator and residual are fp64.
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_precision, only: precision_policy_t, precision_pair_t, &
        STAGE_PRECONDITION, FP64
    use fortnum_krylov, only: real_matvec_t
    implicit none
    private
    public :: mixed_refine, correction_solver_t
    abstract interface
        subroutine correction_solver_t(rhs, delta, precision, info)
            import :: real64, precision_pair_t
            real(real64), intent(in) :: rhs(:)
            real(real64), intent(out) :: delta(:)
            type(precision_pair_t), intent(in) :: precision
            integer, intent(out) :: info
        end subroutine
    end interface
contains
    subroutine mixed_refine(matvec64, correction, b, x, policy, tolerance, &
            max_steps, info, steps, residual)
        procedure(real_matvec_t) :: matvec64
        procedure(correction_solver_t) :: correction
        real(real64), intent(in) :: b(:), tolerance
        real(real64), intent(inout) :: x(:)
        type(precision_policy_t), intent(inout) :: policy
        integer, intent(in) :: max_steps
        integer, intent(out) :: info, steps
        real(real64), intent(out) :: residual
        real(real64), allocatable :: r(:), delta(:), ax(:)
        real(real64) :: target, previous, scale_rhs
        type(precision_pair_t) :: pair
        integer :: inner_info
        logical :: accepted
        info = -1
        steps = 0
        residual = huge(1.0_real64)
        if (size(x) /= size(b) .or. size(b) < 1) return
        if (.not. ieee_is_finite(tolerance)) return
        if (tolerance <= 0 .or. max_steps < 0) return
        if (.not. all(ieee_is_finite(b))) return
        if (.not. all(ieee_is_finite(x))) return
        allocate(r(size(b)), delta(size(b)), ax(size(b)))
        target = tolerance * maxval(abs(b))
        previous = huge(1.0_real64)
        do
            call matvec64(x, ax)
            r = b - ax
            residual = maxval(abs(r))
            info = 2
            if (.not. all(ieee_is_finite(r))) return
            if (residual <= target) then
                info = 0
                return
            end if
            info = 1
            if (steps == max_steps) return
            pair = policy%get(STAGE_PRECONDITION)
            if (residual >= 0.9_real64 * previous) then
                ! Measured stagnation forces rebuild/retry at a higher format.
                call policy%assess(STAGE_PRECONDITION, 1.0_real64, &
                    0.0_real64, accepted)
                pair = policy%get(STAGE_PRECONDITION)
            end if
            previous = residual
            ! Normalize defect before a half solve to prevent overflow/underflow.
            scale_rhs = residual
            call correction(r / scale_rhs, delta, pair, inner_info)
            steps = steps + 1
            if (inner_info /= 0 .or. .not. all(ieee_is_finite(delta))) then
                if (pair%storage == FP64) then
                    info = 3
                    return
                end if
                call policy%assess(STAGE_PRECONDITION, 1.0_real64, &
                    0.0_real64, accepted)
                ! No update was attempted: do not also count the unchanged
                ! defect as stagnation on the retry.
                previous = huge(1.0_real64)
                cycle
            end if
            x = x + scale_rhs * delta
        end do
    end subroutine
end module fortnum_refinement
