program test_fortnum_quadrature
    ! Behavioral tests for fortnum_quadrature.
    ! (1) Weights sum to 2 on [-1,1] for several n.
    ! (2) Nodes lie strictly in (-1, 1) and are symmetric.
    ! (3) gauss_legendre_ab integrates polynomials exactly up to degree 2n-1.
    ! (4) gauss_legendre_ab maps interval correctly (integral of 1 = b - a).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_quadrature, only: gauss_legendre, gauss_legendre_ab
    implicit none

    integer,  parameter :: n_cases(6) = [1, 2, 3, 5, 8, 16]
    real(dp), parameter :: tol = 4.0_dp*epsilon(1.0_dp)   ! near machine eps
    real(dp), parameter :: poly_tol = 1.0e-13_dp           ! polynomial exactness
    integer             :: j, k, n
    real(dp)            :: sum_w, err
    logical             :: ok

    ok = .true.

    ! --- (1) & (2): weight sum and node symmetry ---
    do j = 1, size(n_cases)
        n = n_cases(j)
        block
            real(dp) :: x(n), w(n)
            call gauss_legendre(n, x, w)

            ! Weights sum to 2 (integral of 1 on [-1,1]).
            sum_w = sum(w)
            err   = abs(sum_w - 2.0_dp)
            if (err > tol) then
                write (error_unit, "(a,i0,a,es12.4)") &
                    "FAIL weight sum n=", n, " err=", err
                ok = .false.
            end if

            ! Nodes strictly interior and symmetric.
            do k = 1, n
                if (abs(x(k)) >= 1.0_dp) then
                    write (error_unit, "(a,i0,a,i0)") &
                        "FAIL node outside (-1,1) n=", n, " k=", k
                    ok = .false.
                end if
                if (abs(x(k) + x(n + 1 - k)) > tol) then
                    write (error_unit, "(a,i0,a,i0)") &
                        "FAIL node symmetry n=", n, " k=", k
                    ok = .false.
                end if
            end do

            ! Nodes are ascending.
            do k = 1, n - 1
                if (x(k) >= x(k + 1)) then
                    write (error_unit, "(a,i0,a,i0)") &
                        "FAIL nodes not ascending n=", n, " k=", k
                    ok = .false.
                end if
            end do
        end block
    end do

    ! --- (3): exact integration of monomials x^p, p = 0..2n-1 ---
    ! Exact value: integral_{-1}^{1} x^p dx = 2/(p+1) for even p, 0 for odd p.
    do j = 1, size(n_cases)
        n = n_cases(j)
        block
            real(dp)    :: x(n), w(n), quad_val, exact
            integer     :: p
            call gauss_legendre(n, x, w)
            do p = 0, 2*n - 1
                quad_val = sum(w * x**p)
                if (mod(p, 2) == 0) then
                    exact = 2.0_dp/real(p + 1, dp)
                else
                    exact = 0.0_dp
                end if
                err = abs(quad_val - exact)
                if (err > poly_tol) then
                    write (error_unit, "(a,i0,a,i0,a,es12.4)") &
                        "FAIL exactness n=", n, " degree=", p, " err=", err
                    ok = .false.
                end if
            end do
        end block
    end do

    ! --- (4): gauss_legendre_ab: interval mapping ---
    block
        integer,  parameter :: na = 5
        real(dp), parameter :: a = 0.3_dp, b = 2.7_dp
        real(dp) :: x(na), w(na)
        call gauss_legendre_ab(na, a, b, x, w)

        ! Integral of 1 = b - a.
        err = abs(sum(w) - (b - a))
        if (err > poly_tol) then
            write (error_unit, "(a,es12.4)") "FAIL ab weight sum err=", err
            ok = .false.
        end if

        ! Nodes lie in (a, b).
        do k = 1, na
            if (x(k) <= a .or. x(k) >= b) then
                write (error_unit, "(a,i0,a,f12.8)") &
                    "FAIL ab node out of range k=", k, " x=", x(k)
                ok = .false.
            end if
        end do

        ! Integral of x over [a,b] = (b^2 - a^2)/2.
        err = abs(sum(w*x) - 0.5_dp*(b*b - a*a))
        if (err > poly_tol) then
            write (error_unit, "(a,es12.4)") "FAIL ab integral x err=", err
            ok = .false.
        end if
    end block

    if (.not. ok) then
        write (error_unit, "(a)") "test_fortnum_quadrature: FAILED"
        stop 1
    end if
    write (*, "(a)") "test_fortnum_quadrature: all tests passed"
    stop 0

end program test_fortnum_quadrature
