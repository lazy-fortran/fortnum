program test_fortnum_coulomb
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_special_coulomb, only: coulomb_fg0, coulomb_c0_squared
    implicit none
    real(dp) :: ref(4), val(4), eta, rho, err, worst
    integer :: unit, ios, rows
    logical :: ok
    character(1024) :: path
    call get_command_argument(1, path)
    open (newunit=unit, file=trim(path), status='old')
    rows = 0; worst = 0
    do
        read (unit, *, iostat=ios) eta, rho, ref
        if (ios < 0) exit
        if (ios /= 0) error stop 'malformed Coulomb oracle'
        call coulomb_fg0(eta, rho, val(1), val(2), val(3), val(4), ok)
        if (.not. ok) error stop 'Coulomb function rejected oracle point'
        err = maxval(abs(val - ref)/max(1.0_dp, abs(ref)))
        worst = max(worst, err)
        if (err > 2.0e-10_dp) then
            print *, 'eta,rho,error', eta, rho, err
            print *, 'value', val
            print *, 'reference', ref
            error stop 'Coulomb mpmath oracle'
        end if
        if (abs(val(3)*val(2) - val(1)*val(4) - 1) > 2.0e-10_dp) &
            error stop 'Coulomb Wronskian'
        rows = rows + 1
    end do
    close (unit)
    if (rows < 30) error stop 'insufficient Coulomb coverage'
    call coulomb_fg0(31.0_dp, 1.0_dp, val(1), val(2), val(3), val(4), ok)
    if (ok) error stop 'unsupported Coulomb domain accepted'
    if (coulomb_c0_squared(0.0_dp) /= 1) error stop 'neutral limit'
    print *, 'Coulomb references:', rows, 'worst scaled error:', worst
end program
