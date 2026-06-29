program test_fortnum_gen_laguerre_oracle
    ! Oracle test for fortnum_quadrature::gauss_gen_laguerre.  Compares nodes
    ! and weights against scipy.special.roots_genlaguerre (Golub-Welsch) for
    ! the weight x^alpha exp(-x) on [0,inf) at n=32, alpha in {5/2, 7/2} -- the
    ! orders/alphas libneo transport.f90 calc_D_one_over_nu uses.
    !
    ! CSV layout: alpha_id,alpha,n,i,node,weight  (nodes ascending).
    !
    ! Beyond table agreement the test asserts the rule's defining property
    ! directly: quadrature exactness on the monomials
    !   integral_0^inf x^(alpha+k) exp(-x) dx = Gamma(alpha+k+1), k=0..2n-1,
    ! and the zeroth moment mu0 = Gamma(alpha+1) (DLMF Table 18.3.1).
    !
    ! Tolerance: 1e-12 relative on nodes/weights (the QL eigensolve and scipy's
    ! tridiagonal solve agree to ~1e-13 at n=32); 1e-11 relative on the
    ! monomial integrals (float64-cancellation limited at high k).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_quadrature, only: gauss_gen_laguerre
    implicit none

    character(len=4096) :: path
    integer             :: arglen, argstat, unit, ios
    character(len=512)  :: line

    integer, parameter  :: max_rows = 128
    integer  :: csv_aid(max_rows), csv_n(max_rows), csv_i(max_rows)
    real(dp) :: csv_alpha(max_rows), csv_node(max_rows), csv_weight(max_rows)
    integer  :: nrows, row

    integer, parameter  :: n_rule = 32, n_alpha = 2
    real(dp), parameter :: alphas(n_alpha) = [2.5_dp, 3.5_dp]
    real(dp), parameter :: tol_nw = 1.0e-12_dp
    real(dp), parameter :: tol_mono = 1.0e-11_dp

    integer  :: j, n, i, k, nfail
    real(dp) :: x(n_rule), w(n_rule), alpha
    real(dp) :: err_x, err_w, approx, exact, relerr

    call get_command_argument(1, path, arglen, argstat)
    if (argstat /= 0 .or. arglen == 0) then
        write (error_unit, "(a)") &
            "usage: test_fortnum_gen_laguerre_oracle <gauss_gen_laguerre.csv>"
        stop 1
    end if

    open (newunit=unit, file=path(1:arglen), status="old", action="read", &
        iostat=ios)
    if (ios /= 0) then
        write (error_unit, "(a)") "cannot open: "//path(1:arglen)
        stop 1
    end if
    nrows = 0
    do
        read (unit, "(a)", iostat=ios) line
        if (ios /= 0) exit
        if (len_trim(line) == 0) cycle
        if (line(verify(line, " "):verify(line, " ")) == "#") cycle
        nrows = nrows + 1
        if (nrows > max_rows) then
            write (error_unit, "(a)") "oracle: too many rows"
            stop 1
        end if
        call replace_commas(line)
        read (line, *, iostat=ios) csv_aid(nrows), csv_alpha(nrows), &
            csv_n(nrows), csv_i(nrows), csv_node(nrows), csv_weight(nrows)
        if (ios /= 0) then
            write (error_unit, "(a,i0)") "oracle: parse error at row ", nrows
            stop 1
        end if
    end do
    close (unit)

    nfail = 0

    do j = 1, n_alpha
        alpha = alphas(j)
        n = n_rule
        call gauss_gen_laguerre(n, alpha, x, w)

        ! Table agreement.
        do row = 1, nrows
            if (abs(csv_alpha(row) - alpha) > 1.0e-30_dp) cycle
            if (csv_n(row) /= n) cycle
            i = csv_i(row)
            if (i < 1 .or. i > n) then
                write (error_unit, "(a,i0)") "oracle: bad index i=", i
                nfail = nfail + 1
                cycle
            end if
            err_x = abs(x(i) - csv_node(row))
            err_w = abs(w(i) - csv_weight(row))
            if (err_x > tol_nw*abs(csv_node(row)) .or. &
                err_w > tol_nw*abs(csv_weight(row))) then
                nfail = nfail + 1
                if (nfail == 1) write (error_unit, "(a)") &
                    "oracle FAIL: alpha,i  ours_node ref_node  ours_w ref_w"
                write (error_unit, "(f5.2,1x,i4,1x,4(es22.14,1x))") &
                    alpha, i, x(i), csv_node(row), w(i), csv_weight(row)
            end if
        end do

        ! Monomial exactness int_0^inf x^(alpha+k) e^-x = Gamma(alpha+k+1).
        do k = 0, 2*n - 1
            approx = 0.0_dp
            do i = 1, n
                approx = approx + w(i)*x(i)**k
            end do
            exact = gamma(alpha + real(k, dp) + 1.0_dp)
            relerr = abs(approx - exact)/exact
            if (relerr > tol_mono) then
                nfail = nfail + 1
                write (error_unit, "(a,f5.2,a,i0,a,es12.4)") &
                    "oracle monomial FAIL: alpha=", alpha, " k=", k, &
                    " relerr=", relerr
            end if
        end do
    end do

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " oracle check(s) failed"
        stop 1
    end if

    write (*, "(a,i0,a)") "oracle passed: ", nrows, &
        " table rows + monomial exactness verified"
    stop 0

contains

    subroutine replace_commas(buf)
        character(*), intent(inout) :: buf
        integer :: k
        do k = 1, len(buf)
            if (buf(k:k) == ",") buf(k:k) = " "
        end do
    end subroutine replace_commas

end program test_fortnum_gen_laguerre_oracle
