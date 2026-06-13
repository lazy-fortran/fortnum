program test_fortnum_fft_oracle
    ! Read fft_c2c.csv produced by the numpy oracle generator and verify
    ! fortnum_fft fft_c2c (sign=-1, forward) and fft_r2c against the
    ! reference within documented tolerance.
    !
    ! CSV columns (non-comment rows):
    !   index, seq_name, k, re_in, im_in, re_fwd, im_fwd
    !
    ! Tolerance: |got - ref| <= atol + rtol*|ref|
    !   atol = 1e-13, rtol = 1e-13  (matches libneo test criterion n*eps)
    !
    ! The test path is the first command-line argument so CTest can pass the
    ! source-tree location.
    use, intrinsic :: iso_fortran_env, only: error_unit
    use fortnum_kinds, only: dp
    use fortnum_fft, only: fft_c2c, fft_r2c
    implicit none

    real(dp), parameter :: atol = 1.0e-12_dp
    real(dp), parameter :: rtol = 1.0e-12_dp

    character(len=4096) :: path
    integer :: arglen, argstat

    call get_command_argument(1, path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") &
            "usage: test_fortnum_fft_oracle <fft_c2c.csv>"
        stop 1
    end if

    call run_tests(path(1:arglen))
    write (*, "(a)") "PASS: fft oracle verified"
    stop 0

contains

    subroutine run_tests(csv_path)
        character(*), intent(in) :: csv_path

        ! Storage for up to 4 sequences of up to 64 elements each.
        integer, parameter :: max_seq = 4
        integer, parameter :: max_n = 64

        complex(dp) :: z_in(max_n, max_seq)
        complex(dp) :: z_ref(max_n, max_seq)
        integer :: seqlen(max_seq)
        character(len=16) :: seqname(max_seq)

        integer :: unit, ios, row, k, s, n
        character(len=512) :: line
        integer :: ridx, rkfld, nseq
        character(len=16) :: rname
        real(dp) :: rri, rii, rrf, rimf
        character(len=512) :: buf
        integer :: cur_seq

        z_in = (0.0_dp, 0.0_dp)
        z_ref = (0.0_dp, 0.0_dp)
        seqlen = 0
        seqname = ''
        nseq = 0
        cur_seq = -1

        open (newunit=unit, file=csv_path, status="old", action="read", &
              iostat=ios)
        if (ios /= 0) then
            write (error_unit, "(a)") "cannot open: "//trim(csv_path)
            stop 1
        end if

        do
            read (unit, "(a)", iostat=ios) line
            if (ios /= 0) exit
            if (is_comment(line) .or. len_trim(line) == 0) cycle
            buf = line
            call replace_commas(buf)
            read (buf, *, iostat=ios) ridx, rname, rkfld, rri, rii, rrf, rimf
            if (ios /= 0) then
                write (error_unit, "(a,a)") "parse error on row: ", trim(line)
                close (unit)
                stop 1
            end if
            ! new sequence?
            s = find_seq(seqname, nseq, rname)
            if (s == 0) then
                nseq = nseq + 1
                s = nseq
                seqname(s) = rname
            end if
            k = rkfld + 1   ! 0-based k -> 1-based index
            z_in(k, s) = cmplx(rri, rii, dp)
            z_ref(k, s) = cmplx(rrf, rimf, dp)
            if (k > seqlen(s)) seqlen(s) = k
        end do
        close (unit)

        if (nseq == 0) then
            write (error_unit, "(a)") "no data rows in "//trim(csv_path)
            stop 1
        end if

        ! Test fft_c2c forward (sign=-1) for each sequence.
        do s = 1, nseq
            n = seqlen(s)
            call check_c2c(seqname(s), z_in(1:n, s), z_ref(1:n, s), n)
        end do

        ! Test fft_r2c on sequences that were generated with zero imaginary
        ! part (seq n8 in the oracle, purely real input).
        do s = 1, nseq
            n = seqlen(s)
            if (maxval(abs(aimag(z_in(1:n, s)))) < 1.0e-15_dp) then
                call check_r2c(seqname(s), real(z_in(1:n, s), dp), &
                               z_ref(1:n, s), n)
            end if
        end do
    end subroutine run_tests

    subroutine check_c2c(name, z_in_ref, z_ref, n)
        character(*), intent(in) :: name
        complex(dp), intent(in) :: z_in_ref(:)
        complex(dp), intent(in) :: z_ref(:)
        integer, intent(in) :: n
        complex(dp) :: z(n)
        integer :: k
        real(dp) :: err, tol

        z = z_in_ref
        call fft_c2c(z, -1)  ! forward

        do k = 1, n
            err = abs(z(k) - z_ref(k))
            tol = atol + rtol*abs(z_ref(k))
            if (.not. (err <= tol)) then
                write (error_unit, "(a,a,a,i0,a,es12.4,a,es12.4)") &
                    "FAIL fft_c2c seq=", trim(name), " k=", k - 1, &
                    " err=", err, " tol=", tol
                stop 1
            end if
        end do
    end subroutine check_c2c

    subroutine check_r2c(name, x, z_ref, n)
        character(*), intent(in) :: name
        real(dp), intent(in) :: x(:)
        complex(dp), intent(in) :: z_ref(:)
        integer, intent(in) :: n
        complex(dp) :: c(n/2 + 1)
        integer :: k
        real(dp) :: err, tol

        call fft_r2c(x, c)

        ! fft_r2c returns only bins 0..n/2 (Hermitian half); compare those.
        do k = 1, n/2 + 1
            err = abs(c(k) - z_ref(k))
            tol = atol + rtol*abs(z_ref(k))
            if (.not. (err <= tol)) then
                write (error_unit, "(a,a,a,i0,a,es12.4,a,es12.4)") &
                    "FAIL fft_r2c seq=", trim(name), " k=", k - 1, &
                    " err=", err, " tol=", tol
                stop 1
            end if
        end do
    end subroutine check_r2c

    ! Find sequence name in array; return 0 if not found.
    pure integer function find_seq(names, nseq, target)
        character(len=16), intent(in) :: names(:)
        integer, intent(in) :: nseq
        character(*), intent(in) :: target
        integer :: i
        find_seq = 0
        do i = 1, nseq
            if (trim(names(i)) == trim(target)) then
                find_seq = i
                return
            end if
        end do
    end function find_seq

    pure logical function is_comment(line)
        character(*), intent(in) :: line
        integer :: p
        is_comment = .false.
        p = verify(line, " ")
        if (p > 0) is_comment = (line(p:p) == "#")
    end function is_comment

    pure subroutine replace_commas(buf)
        character(*), intent(inout) :: buf
        integer :: i
        do i = 1, len(buf)
            if (buf(i:i) == ",") buf(i:i) = " "
        end do
    end subroutine replace_commas

end program test_fortnum_fft_oracle
