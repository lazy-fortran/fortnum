module sampled_vector_cases
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_interval, only: interval_t, ipoint
    use fortnum_sampled_norm, only: coordinate_enclosure_t
    implicit none
    type, extends(coordinate_enclosure_t) :: sample_vector_t
        real(dp), allocatable :: values(:)
    contains
        procedure :: evaluate => component
    end type
contains
    subroutine component(self, index, value, ok)
        class(sample_vector_t), intent(in) :: self
        integer, intent(in) :: index
        type(interval_t), intent(out) :: value
        logical, intent(out) :: ok
        value = ipoint(self%values(index)); ok = .true.
    end subroutine
end module

program test_sampled_norm
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, i64 => int64
    use fortnum_interval, only: interval_t, ipoint
    use fortnum_rng, only: rng_t, rng_seed
    use fortnum_status, only: fortnum_status_t
    use fortnum_sampled_norm
    use sampled_vector_cases, only: sample_vector_t
    implicit none
    integer, parameter :: n=257, repetitions=2000
    real(dp) :: estimate, upper, exact, rate
    type(sample_vector_t) :: sample
    type(rng_t) :: draw
    type(fortnum_status_t) :: status
    integer :: family,k,j,failures
    logical :: ok
    allocate(sample%values(n))
    do family=1,3
        sample%values=1
        if (family==2) sample%values=[(real(mod(j,7),dp)/6,j=1,n)]
        if (family==3) then
            sample%values=0; sample%values(1)=1
        end if
        exact=real(sqrt(sum(real(sample%values,qp)**2)),dp)
        failures=0
        do k=1,repetitions
            call rng_seed(draw,int(k,i64),status)
            call sampled_norm_upper(n,128,.05_dp,1._dp,draw,sample,estimate,upper,ok)
            if (.not. ok) error stop 'sampler failed'
            if (upper<exact) failures=failures+1
        end do
        rate=real(failures,dp)/repetitions
        print '(a,i0,a,f8.5)', 'family=',family,' observed failure rate=',rate
        ! Conservative independent binomial margin: delta + 5 sigma.
        if (rate>.05_dp+5*sqrt(.05_dp*.95_dp/repetitions)) error stop 'coverage'
    end do
    call sampled_norm_upper(n,128,0._dp,1._dp,draw,sample,estimate,upper,ok)
    if (ok) error stop 'invalid confidence'
    print *, 'PASS: constant, heterogeneous, sparse norms and coverage'
end program
