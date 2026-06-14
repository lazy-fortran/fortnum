program test_fortnum_special_complex_bessel_oracle
    ! Oracle test for fortnum_special_complex_bessel against scipy.special
    ! jv / iv / ive / kv / kve with complex argument.
    ! CSV columns: index,func,n,re_z,im_z,re_expected,im_expected
    !   func 0 -> J_n(z)              bessel_j_complex
    !   func 1 -> I_n(z)              bessel_i_complex (scaled=.false.)
    !   func 2 -> e^{-Re z} I_n(z)    bessel_i_complex (scaled=.true.)
    !   func 3 -> K_n(z)              bessel_k_complex (scaled=.false.)
    !   func 4 -> e^{z} K_n(z)        bessel_k_complex (scaled=.true.)
    !
    ! Tolerance 3e-12 relative.  The binding constraint is J near its
    ! series/asymptotic crossover at |z| ~ 14, where the DLMF 10.2.2 power
    ! series loses ~2e-12 to cancellation; I and K reach machine precision.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_special_complex_bessel, only: bessel_j_complex, &
        bessel_i_complex, bessel_k_complex

    implicit none

    real(dp), parameter :: rel_tol = 3.0e-12_dp

    character(len=4096) :: csv_path
    integer             :: arglen, argstat, unit, ios, row_idx
    character(len=512)  :: line

    integer     :: func, n, nfail
    real(dp)    :: re_z, im_z, re_e, im_e, err, scale
    complex(dp) :: z, expected, got
    real(dp)    :: worst_j, worst_i, worst_k
    type(fortnum_status_t) :: status

    call get_command_argument(1, csv_path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") &
            "usage: test_fortnum_special_complex_bessel_oracle <csv>"
        stop 1
    end if

    open (newunit=unit, file=csv_path(1:arglen), status="old", &
          action="read", iostat=ios)
    if (ios /= 0) then
        write (error_unit, "(a)") "cannot open "//csv_path(1:arglen)
        stop 1
    end if

    nfail   = 0
    row_idx = 0
    worst_j = 0.0_dp
    worst_i = 0.0_dp
    worst_k = 0.0_dp

    do
        read (unit, "(a)", iostat=ios) line
        if (ios /= 0) exit
        if (is_comment(line) .or. len_trim(line) == 0) cycle

        call parse_row(line, row_idx, func, n, re_z, im_z, re_e, im_e, ios)
        if (ios /= 0) then
            write (error_unit, "(a,i0)") "malformed row at index ", row_idx
            nfail = nfail + 1
            cycle
        end if

        z        = cmplx(re_z, im_z, kind=dp)
        expected = cmplx(re_e, im_e, kind=dp)

        select case (func)
        case (0)
            call bessel_j_complex(n, z, got, status)
        case (1)
            call bessel_i_complex(n, z, .false., got, status)
        case (2)
            call bessel_i_complex(n, z, .true., got, status)
        case (3)
            call bessel_k_complex(n, z, .false., got, status)
        case (4)
            call bessel_k_complex(n, z, .true., got, status)
        case default
            write (error_unit, "(a,i0)") "unknown func at row ", row_idx
            nfail = nfail + 1
            cycle
        end select

        if (.not. status_ok(status)) then
            write (error_unit, "(a,i0)") "status error at row ", row_idx
            nfail = nfail + 1
            cycle
        end if

        scale = max(abs(expected), 1.0e-300_dp)
        err   = abs(got - expected)/scale
        select case (func)
        case (0)
            worst_j = max(worst_j, err)
        case (1, 2)
            worst_i = max(worst_i, err)
        case (3, 4)
            worst_k = max(worst_k, err)
        end select

        if (err > rel_tol) then
            nfail = nfail + 1
            write (error_unit, "(a,i0,a,i0,a,i0,a,es12.5,a,es12.5,a,es12.5)") &
                "FAIL row ", row_idx, " func=", func, " n=", n, &
                " re_z=", re_z, " im_z=", im_z, " relerr=", err
        end if
    end do
    close (unit)

    write (*, "(a,es13.6)") "worst J rel err: ", worst_j
    write (*, "(a,es13.6)") "worst I rel err: ", worst_i
    write (*, "(a,es13.6)") "worst K rel err: ", worst_k

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " failures"
        stop 1
    end if

    write (*, "(a)") "PASS"
    stop 0

contains

    pure logical function is_comment(line)
        character(*), intent(in) :: line
        integer :: p
        is_comment = .false.
        p = verify(line, " ")
        if (p > 0) is_comment = (line(p:p) == "#")
    end function is_comment

    subroutine parse_row(line, idx, func, n, re_z, im_z, re_e, im_e, ios)
        character(*), intent(in)  :: line
        integer,      intent(out) :: idx, func, n, ios
        real(dp),     intent(out) :: re_z, im_z, re_e, im_e

        character(len=len(line)) :: buf
        integer :: i

        buf = line
        do i = 1, len(buf)
            if (buf(i:i) == ",") buf(i:i) = " "
        end do
        read (buf, *, iostat=ios) idx, func, n, re_z, im_z, re_e, im_e
    end subroutine parse_row

end program test_fortnum_special_complex_bessel_oracle
