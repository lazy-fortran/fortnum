program test_fortnum_precision
    use iso_fortran_env, only: real32, real64, int64
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use fortnum_precision
    use fortnum_precision_vector, only: precision_vector_t
    use fortnum_refinement, only: mixed_refine
    implicit none
    type(precision_vector_t) :: v
    type(precision_policy_t) :: policy
    type(precision_pair_t) :: pair
    real(real64) :: a(513), actual(513), x(3), b(3), residual, bound
    real(real32) :: narrow(513)
    integer, parameter :: formats(4)=[FP64,FP32,FP16,BF16], widths(4)=[8,4,2,2]
    integer :: f, i, info, steps, stage
    logical :: ok, reject_half
    do i=1,513
        a(i) = sin(real(i,real64)) * 1.0e100_real64
    end do
    do f=1,4
        call v%store(a,formats(f),ok)
        if (.not. ok) error stop 'store'
        call v%load(actual,ok)
        if (.not. ok) error stop 'load'
        bound = unit_roundoff(formats(f))*maxval(abs(a))
        if (maxval(abs(a-actual)) > bound) error stop 'scaled roundoff'
        if (v%bytes() /= int(513,int64)*widths(f)) error stop 'payload bytes'
        call v%load(narrow,ok)
        if (ok) error stop 'unrepresentable fp32 compute accepted'
    end do
    a = [(real(i,real64)/513, i=1,513)]
    call v%store(a,FP16,ok)
    call v%load(narrow,ok)
    if (.not. ok) error stop 'fp32 compute load'
    if (maxval(abs(real(narrow,real64)-a)) > unit_roundoff(FP16)) &
        error stop 'separate compute'
    actual = a
    actual(1) = ieee_value(0.0_real64,ieee_quiet_nan)
    call v%store(actual,FP16,ok)
    if (ok) error stop 'nan storage accepted'
    call v%load(actual,ok)
    if (maxval(abs(actual-a)) > unit_roundoff(FP16)) &
        error stop 'failed store modified vector'
    actual = huge(1.0_real64)
    call v%store(actual,FP16,ok)
    if (ok) error stop 'half reconstruction overflow accepted'
    policy = precision_preset('cpu-mixed',ok)
    if (.not. ok) error stop 'preset'
    do stage=STAGE_RESIDUAL,STAGE_CERTIFY
        call policy%set(stage,FP32,FP32,ok)
        if (ok) error stop 'fp64 lock'
        pair = policy%get(stage)
        if (pair%storage /= FP64 .or. pair%compute /= FP64) &
            error stop 'fp64 stage changed'
    end do
    call policy%assess(STAGE_PRECONDITION,1.0_real64,1.0_real64,ok)
    if (ok) error stop 'budget not enforced'
    pair = policy%get(STAGE_PRECONDITION)
    if (pair%storage /= FP32) error stop 'half escalation'
    call policy%assess(STAGE_PRECONDITION,1.0_real64,1.0_real64,ok)
    pair = policy%get(STAGE_PRECONDITION)
    if (pair%storage /= FP64) error stop 'single escalation'
    policy = precision_preset('cpu-mixed',ok)
    reject_half = .true.
    x = 0
    b = [2.0_real64, 6.0_real64, 12.0_real64]
    call mixed_refine(product,correction,b,x,policy,32*epsilon(1.0_real64), &
        20,info,steps,residual)
    if (info /= 0) error stop 'refinement convergence'
    pair = policy%get(STAGE_PRECONDITION)
    if (pair%storage /= FP32) error stop 'failed solve did not escalate'
    if (maxval(abs(x-[1.0_real64,2.0_real64,3.0_real64])) > &
        32*epsilon(1.0_real64)) error stop 'refinement accuracy'
    call mixed_refine(product,correction,b,x,policy,-1.0_real64, &
        20,info,steps,residual)
    if (info /= -1) error stop 'invalid tolerance'
    print *, 'precision policy, scaled vector and fp64 defect correction PASS'
contains
    subroutine product(input,output)
        real(real64), intent(in) :: input(:)
        real(real64), intent(out) :: output(:)
        output = [2.0_real64,3.0_real64,4.0_real64]*input
    end subroutine
    subroutine correction(rhs,delta,precision,info)
        real(real64), intent(in) :: rhs(:)
        real(real64), intent(out) :: delta(:)
        type(precision_pair_t), intent(in) :: precision
        integer, intent(out) :: info
        if (reject_half .and. precision%storage == FP16) then
            reject_half = .false.
            delta = 0
            info = 1
            return
        end if
        ! Independent low precision diagonal solve; no iteration differentiation.
        delta = real(real(rhs,real32)/[2.0_real32,3.0_real32,4.0_real32],real64)
        if (precision%compute == FP64) &
            delta = rhs/[2.0_real64,3.0_real64,4.0_real64]
        info = 0
    end subroutine
end program
