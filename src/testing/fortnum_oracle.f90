module fortnum_oracle
    ! Read reference ("oracle") tables produced by test/oracle/gen_oracle.py
    ! and check a Fortran primal (and, later, derivative) implementation
    ! against them within absolute+relative tolerance. The file format is
    ! frozen so derivative reference tables (issue #40) drop in without an
    ! API or reader change.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none
    private

    public :: oracle_table_t
    public :: oracle_read
    public :: oracle_check
    public :: oracle_primal_fn

    ! Loaded reference table. has_derivative mirrors the file header flag so
    ! callers know whether the derivative column carries meaning.
    type :: oracle_table_t
        character(:),  allocatable :: name
        real(dp),      allocatable :: x(:)
        real(dp),      allocatable :: primal(:)
        real(dp),      allocatable :: derivative(:)
        logical                    :: has_derivative = .false.
    end type oracle_table_t

    ! Shape every checkable primal must present: y = f(x).
    abstract interface
        pure function oracle_primal_fn(x) result(y)
            import :: dp
            real(dp), intent(in) :: x
            real(dp)             :: y
        end function oracle_primal_fn
    end interface

contains

    ! Read a CSV oracle table. Lines beginning with '#' are header/comment;
    ! the "has_derivative:" header sets table%has_derivative. Data rows are
    ! "index,x,primal,derivative". On any I/O or parse failure the status
    ! carries a domain error and the table is returned unallocated.
    subroutine oracle_read(path, table, status)
        character(*),         intent(in)  :: path
        type(oracle_table_t), intent(out) :: table
        type(fortnum_status_t), intent(out) :: status

        integer :: unit, ios, n, i
        character(len=512) :: line
        real(dp), allocatable :: xs(:), ps(:), ds(:)

        call status_set(status, FORTNUM_OK, "")
        table%has_derivative = .false.

        open (newunit=unit, file=path, status="old", action="read", &
            iostat=ios)
        if (ios /= 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "oracle_read: cannot open "//trim(path))
            return
        end if

        ! First pass: count data rows and read header flags.
        n = 0
        do
            read (unit, "(a)", iostat=ios) line
            if (ios /= 0) exit
            if (is_comment(line)) then
                call scan_header(line, table%has_derivative)
            else if (len_trim(line) > 0) then
                n = n + 1
            end if
        end do

        if (n == 0) then
            close (unit)
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "oracle_read: no data rows in "//trim(path))
            return
        end if

        allocate (xs(n), ps(n), ds(n))

        ! Second pass: parse the data rows in order.
        rewind (unit)
        i = 0
        do
            read (unit, "(a)", iostat=ios) line
            if (ios /= 0) exit
            if (is_comment(line) .or. len_trim(line) == 0) cycle
            i = i + 1
            call parse_row(line, xs(i), ps(i), ds(i), ios)
            if (ios /= 0) then
                close (unit)
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "oracle_read: malformed row in "//trim(path))
                return
            end if
        end do
        close (unit)

        table%name = trim(path)
        call move_alloc(xs, table%x)
        call move_alloc(ps, table%primal)
        call move_alloc(ds, table%derivative)
    end subroutine oracle_read

    ! Check a primal function against a loaded table. Every entry must satisfy
    ! |got - expected| <= atol + rtol*|expected|. On the first miss the status
    ! carries a convergence error and a report (index, input, expected, got,
    ! error) is written to error_unit; checking still scans all rows so the
    ! report lists every failure.
    subroutine oracle_check(table, f, atol, rtol, status)
        type(oracle_table_t),    intent(in)  :: table
        procedure(oracle_primal_fn)           :: f
        real(dp),                intent(in)  :: atol
        real(dp),                intent(in)  :: rtol
        type(fortnum_status_t),  intent(out) :: status

        integer  :: i, nfail
        real(dp) :: got, expected, err, tol

        call status_set(status, FORTNUM_OK, "")
        nfail = 0

        do i = 1, size(table%x)
            expected = table%primal(i)
            got = f(table%x(i))
            err = abs(got - expected)
            tol = atol + rtol*abs(expected)
            if (.not. (err <= tol)) then
                nfail = nfail + 1
                if (nfail == 1) then
                    write (error_unit, "(a)") &
                        "oracle_check FAIL: "//trim(table%name)
                    write (error_unit, "(a)") &
                        "  idx            input         expected"// &
                        "              got            abserr"
                end if
                write (error_unit, "(2x,i4,4(1x,es24.16e3))") &
                    i - 1, table%x(i), expected, got, err
            end if
        end do

        if (nfail > 0) then
            call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                "oracle_check: entries outside tolerance")
        end if
    end subroutine oracle_check

    ! .true. for header/comment lines (first non-blank char is '#').
    pure logical function is_comment(line)
        character(*), intent(in) :: line
        integer :: p
        is_comment = .false.
        p = verify(line, " ")
        if (p > 0) is_comment = (line(p:p) == "#")
    end function is_comment

    ! Set has_derivative from a "# has_derivative: <0|1>" header line. Other
    ! header lines are ignored.
    subroutine scan_header(line, has_derivative)
        character(*), intent(in)    :: line
        logical,      intent(inout) :: has_derivative
        integer :: k
        k = index(line, "has_derivative:")
        if (k > 0) then
            has_derivative = (index(line(k:), "1") > 0)
        end if
    end subroutine scan_header

    ! Parse "index,x,primal,derivative" into the three floats (index is
    ! positional and discarded). ios /= 0 on any conversion failure.
    subroutine parse_row(line, x, primal, derivative, ios)
        character(*), intent(in)  :: line
        real(dp),     intent(out) :: x, primal, derivative
        integer,      intent(out) :: ios
        character(len=len(line)) :: buf

        ! List-directed read needs separators; turn commas into spaces and let
        ! the runtime convert the four fields. The leading index is read into
        ! a throwaway integer.
        integer :: idx
        buf = line
        call replace_commas(buf)
        read (buf, *, iostat=ios) idx, x, primal, derivative
    end subroutine parse_row

    ! In-place comma -> space so list-directed input can tokenise the row.
    pure subroutine replace_commas(buf)
        character(*), intent(inout) :: buf
        integer :: i
        do i = 1, len(buf)
            if (buf(i:i) == ",") buf(i:i) = " "
        end do
    end subroutine replace_commas

end module fortnum_oracle
