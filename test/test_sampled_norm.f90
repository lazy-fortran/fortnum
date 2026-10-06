program test_sampled_norm
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, i64 => int64
    use fortnum_interval, only: interval_t, ipoint
    use fortnum_rng, only: rng_t, rng_seed
    use fortnum_status, only: fortnum_status_t
    use fortnum_sampled_norm
    implicit none
    integer, parameter :: n=257, repetitions=2000
    real(dp) :: vector(n), estimate, upper, exact, rate
    type(rng_t) :: draw
    type(fortnum_status_t) :: status
    integer :: family,k,j,failures
    logical :: ok
    do family=1,3
        vector=1
        if (family==2) vector=[(real(mod(j,7),dp)/6,j=1,n)]
        if (family==3) then
            vector=0; vector(1)=1
        end if
        exact=real(sqrt(sum(real(vector,qp)**2)),dp)
        failures=0
        do k=1,repetitions
            call rng_seed(draw,int(k,i64),status)
            call sampled_norm_upper(n,128,.05_dp,1._dp,draw,component,estimate,upper,ok)
            if (.not. ok) error stop 'sampler failed'
            if (upper<exact) failures=failures+1
        end do
        rate=real(failures,dp)/repetitions
        print '(a,i0,a,f8.5)', 'family=',family,' observed failure rate=',rate
        ! Conservative independent binomial margin: delta + 5 sigma.
        if (rate>.05_dp+5*sqrt(.05_dp*.95_dp/repetitions)) error stop 'coverage'
    end do
    call sampled_norm_upper(n,128,0._dp,1._dp,draw,component,estimate,upper,ok)
    if (ok) error stop 'invalid confidence'
    print *, 'PASS: constant, heterogeneous, sparse norms and coverage'
contains
    subroutine component(index,value,ok)
        integer, intent(in) :: index
        type(interval_t), intent(out) :: value
        logical, intent(out) :: ok
        value=ipoint(vector(index)); ok=.true.
    end subroutine
end program
