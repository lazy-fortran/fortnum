program test_verified_matvec
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128
    use fortnum_interval, only: interval_t, interval
    use fortnum_verified_matvec
    implicit none
    real(dp) :: a(4,4)
    type(interval_t) :: x(4), y(4)
    type(matvec_bound_t) :: plan
    real(qp) :: exact(4), corner(4)
    integer :: i,j,k
    logical :: ok
    a=reshape([1._dp,1.e100_dp,1._dp,-1._dp, -1._dp,1._dp,1.e-100_dp,2._dp, &
        1._dp,-1.e100_dp,3._dp,1._dp, 2._dp,1._dp,-1._dp,2._dp],[4,4])
    x=[interval(1._dp,nearest(1._dp,1._dp)),interval(-.5_dp,.25_dp), &
        interval(1._dp),interval(-2._dp,-1._dp)]
    call prepare_matvec_bound(a,plan,ok)
    if (.not. ok) error stop 'prepare'
    call enclose_matvec(a,plan,x,y,ok)
    if (.not. ok) error stop 'enclose'
    ! Every box vertex, accumulated independently in binary128.
    do k=0,15
        do j=1,4
            corner(j)=real(x(j)%lo,qp)
            if (btest(k,j-1)) corner(j)=real(x(j)%hi,qp)
        end do
        exact=0
        do j=1,4
            do i=1,4
                exact(i)=exact(i)+real(a(i,j),qp)*corner(j)
            end do
        end do
        if (any(exact<real(y%lo,qp)).or.any(exact>real(y%hi,qp))) &
            error stop 'independent vertex containment'
    end do
    a=tiny(0._dp); x=interval(tiny(0._dp))
    call prepare_matvec_bound(a,plan,ok)
    call enclose_matvec(a,plan,x,y,ok)
    exact=4*real(tiny(0._dp),qp)**2
    if (.not. ok) error stop 'gradual underflow'
    if (any(exact<real(y%lo,qp)).or.any(exact>real(y%hi,qp))) error stop 'underflow'
    a=huge(0._dp)/2; x=interval(4._dp)
    call prepare_matvec_bound(a,plan,ok)
    if (ok) then
        call enclose_matvec(a,plan,x,y,ok)
        if (ok) error stop 'overflow must fail'
    end if
    print *, 'PASS: binary128 vertices, cancellation, underflow, overflow'
end program
