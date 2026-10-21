module fortnum_precision
    !! Per-stage policy. Accepted residuals and certificates cannot be lowered.
    use, intrinsic :: iso_fortran_env, only: real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortgen_precision, only: FP64 => STORAGE_FP64, FP32 => STORAGE_FP32, &
        FP16 => STORAGE_FP16, BF16 => STORAGE_BF16
    implicit none
    private
    public :: FP64, FP32, FP16, BF16
    integer, parameter, public :: STAGE_EVALUATE=1, STAGE_ASSEMBLE=2, &
        STAGE_PRECONDITION=3, STAGE_KRYLOV=4, STAGE_RESIDUAL=5, &
        STAGE_ESTIMATE=6, STAGE_CERTIFY=7
    type, public :: precision_pair_t
        integer :: storage=FP64, compute=FP64
    end type
    type, public :: precision_policy_t
        private
        type(precision_pair_t) :: stages(7)
        real(real64) :: fraction=0.1_real64
    contains
        procedure :: get => policy_get
        procedure :: set => policy_set
        procedure :: assess => policy_assess
        procedure :: set_budget => policy_budget
    end type
    public :: precision_preset, unit_roundoff, valid_storage
contains
    pure logical function valid_storage(format) result(valid)
        integer, intent(in) :: format
        valid = any(format == [FP64, FP32, FP16, BF16])
    end function

    pure real(real64) function unit_roundoff(format) result(u)
        integer, intent(in) :: format
        u = huge(1.0_real64)
        select case(format)
        case(FP64); u = 2.0_real64**(-53)
        case(FP32); u = 2.0_real64**(-24)
        case(FP16); u = 2.0_real64**(-11)
        case(BF16); u = 2.0_real64**(-8)
        end select
    end function

    function precision_preset(name, ok) result(policy)
        character(*), intent(in) :: name
        logical, intent(out) :: ok
        type(precision_policy_t) :: policy
        ok = .true.
        select case(name)
        case('cpu-fp64', 'gpu-hpc')
        case('cpu-mixed', 'gpu-consumer')
            policy%stages(1:4) = precision_pair_t(FP32, FP32)
            policy%stages(STAGE_PRECONDITION) = precision_pair_t(FP16, FP32)
        case default
            ok = .false.
        end select
    end function

    pure function policy_get(self, stage) result(pair)
        class(precision_policy_t), intent(in) :: self
        integer, intent(in) :: stage
        type(precision_pair_t) :: pair
        ! Invalid stages fail closed to the fp64 default.
        if (stage >= 1 .and. stage <= 7) pair = self%stages(stage)
    end function

    subroutine policy_set(self, stage, storage, compute, ok)
        class(precision_policy_t), intent(inout) :: self
        integer, intent(in) :: stage, storage, compute
        logical, intent(out) :: ok
        ok = .false.
        if (stage < 1 .or. stage > 7) return
        if (.not. valid_storage(storage)) return
        if (compute /= FP32 .and. compute /= FP64) return
        if (storage == FP64 .and. compute /= FP64) return
        if (stage >= STAGE_RESIDUAL) then
            if (storage /= FP64 .or. compute /= FP64) return
        end if
        self%stages(stage) = precision_pair_t(storage, compute)
        ok = .true.
    end subroutine

    subroutine policy_budget(self, fraction, ok)
        class(precision_policy_t), intent(inout) :: self
        real(real64), intent(in) :: fraction
        logical, intent(out) :: ok
        ok = .false.
        if (.not. ieee_is_finite(fraction)) return
        if (fraction <= 0 .or. fraction > 1) return
        self%fraction = fraction
        ok = .true.
    end subroutine

    subroutine policy_assess(self, stage, rounding_error, model_error, accepted)
        !! Caller supplies an a posteriori error after sensitivity amplification.
        !! A rejected stage must be rerun, including rebuilding its stored data.
        class(precision_policy_t), intent(inout) :: self
        integer, intent(in) :: stage
        real(real64), intent(in) :: rounding_error, model_error
        logical, intent(out) :: accepted
        accepted = .false.
        if (stage < 1 .or. stage > 7) return
        if (.not. ieee_is_finite(rounding_error)) return
        if (.not. ieee_is_finite(model_error)) return
        if (rounding_error < 0 .or. model_error < 0) return
        accepted = rounding_error <= self%fraction * model_error
        if (accepted .or. stage >= STAGE_RESIDUAL) return
        select case(self%stages(stage)%storage)
        case(FP16, BF16)
            self%stages(stage) = precision_pair_t(FP32, FP32)
        case default
            self%stages(stage) = precision_pair_t(FP64, FP64)
        end select
    end subroutine
end module fortnum_precision
