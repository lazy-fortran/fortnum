program test_fortnum_region_roots_noisy
    ! complex_region_roots on an analytic function sampled with relative noise.
    !
    ! Oracle: f(z) = exp(z) (z - z1)(z - z2)(z - z3)^2 (z - z4) (z - z5) with
    ! the zeros placed by construction; z1..z4 lie in the box [-1,1]^2 (z3 is
    ! double), z5 outside.  Every evaluation gets eta |exp(z)| (u + i v) added,
    ! eta = 1e-7 and u, v a deterministic hash of z in [-1, 1]: noise at a
    ! fixed fraction of the function scale, as in an ODE-based dispersion
    ! determinant, which does not vanish at the zeros.
    !
    ! Case A (defaults): the tight default contour quadrature cannot settle
    ! on the noise, so the default search either fails or exhausts the
    ! evaluation budget.  It is only run once and is capped by the budget.
    ! Case B (contour_tol, fd_step, newtonf for the noise level): exactly
    ! the four zeros, multiplicities 1,1,2,1, within the noise-limited
    ! accuracy, inside the budget.
    ! Case C (exact f, derivative callback df, contour_tol): the same zeros to
    ! 1e-10 with fewer evaluations than the finite-difference variant.
    ! Case D (zero-free box with the options): no zeros, and the zero count
    ! pass alone settles it (far fewer evaluations than a box with zeros).
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_roots_complex, only: complex_region_roots
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    real(dp), parameter :: eta = 1.0e-7_dp
    integer,  parameter :: budget = 40000
    complex(dp), parameter :: zexact(5) = [ &
        (0.3_dp, 0.4_dp), (-0.5_dp, -0.2_dp), (0.1_dp, -0.6_dp), &
        (0.98_dp, 0.25_dp), (1.5_dp, 0.0_dp)]
    integer, parameter :: mexact(5) = [1, 1, 2, 1, 1]
    integer, parameter :: ninside = 4

    integer, save :: ncalls = 0
    logical, save :: noisy = .true.
    integer :: nfail

    nfail = 0
    call case_defaults_fail(nfail)
    call case_noisy_options(nfail)
    call case_derivative_callback(nfail)
    call case_zero_free(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " noisy region-roots check(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"

contains

    pure function hash_unit(x) result(u)
        real(dp), intent(in) :: x
        real(dp) :: u, s
        s = sin(x)*43758.5453123_dp
        u = 2.0_dp*(s - floor(s)) - 1.0_dp
    end function hash_unit

    function exact_f(z) result(fz)
        complex(dp), intent(in) :: z
        complex(dp) :: fz
        integer :: k
        fz = exp(z)
        do k = 1, size(zexact)
            fz = fz*(z - zexact(k))**mexact(k)
        end do
    end function exact_f

    subroutine f_sampled(z, fz, ctx)
        complex(dp), intent(in)  :: z
        complex(dp), intent(out) :: fz
        class(*),    intent(in), optional :: ctx
        real(dp) :: u, v
        ncalls = ncalls + 1
        if (ncalls > budget) then
            fz = (1.0_dp, 0.0_dp)
            return
        end if
        fz = exact_f(z)
        if (noisy) then
            u = hash_unit(12.9898_dp*1.0e4_dp*real(z) + 78.233_dp*1.0e4_dp*aimag(z))
            v = hash_unit(39.3468_dp*1.0e4_dp*real(z) + 11.135_dp*1.0e4_dp*aimag(z))
            fz = fz + eta*abs(exp(z))*cmplx(u, v, dp)
        end if
        if (present(ctx)) continue
    end subroutine f_sampled

    ! f' = f (1 + sum_k m_k/(z - z_k)) for the exact f.
    subroutine df_exact(z, dfz, ctx)
        complex(dp), intent(in)  :: z
        complex(dp), intent(out) :: dfz
        class(*),    intent(in), optional :: ctx
        complex(dp) :: lsum
        integer :: k
        ncalls = ncalls + 1
        lsum = (1.0_dp, 0.0_dp)
        do k = 1, size(zexact)
            lsum = lsum + mexact(k)/(z - zexact(k))
        end do
        dfz = exact_f(z)*lsum
        if (present(ctx)) continue
    end subroutine df_exact

    subroutine case_defaults_fail(nfail)
        integer, intent(inout) :: nfail
        complex(dp), allocatable :: roots(:), fvals(:)
        integer, allocatable :: mult(:)
        integer :: nfound
        type(fortnum_status_t) :: status

        noisy = .true.
        ncalls = 0
        call complex_region_roots(f_sampled, (-1.0_dp, -1.0_dp), (1.0_dp, 1.0_dp), &
            roots, fvals, mult, nfound, status)
        write (*, "(a,i0,a,i0,2a)") "A defaults: nfound=", nfound, " calls=", ncalls, &
            " status: ", trim(status%msg)
        if (status%code == FORTNUM_OK .and. ncalls <= budget) then
            nfail = nfail + 1
            write (error_unit, "(a)") &
                "FAIL A: defaults settled the noisy contour within the budget"
        end if
    end subroutine case_defaults_fail

    subroutine case_noisy_options(nfail)
        integer, intent(inout) :: nfail
        complex(dp), allocatable :: roots(:), fvals(:)
        integer, allocatable :: mult(:)
        integer :: nfound
        type(fortnum_status_t) :: status

        noisy = .true.
        ncalls = 0
        call complex_region_roots(f_sampled, (-1.0_dp, -1.0_dp), (1.0_dp, 1.0_dp), &
            roots, fvals, mult, nfound, status, newtonz=1.0e-9_dp, newtonf=1.0e-6_dp, &
            contour_tol=1.0e-4_dp, fd_step=1.0e-3_dp)
        call check_zero_set("B noisy+options", roots, mult, nfound, status, &
            [1.0e-5_dp, 1.0e-5_dp, 3.0e-3_dp, 1.0e-5_dp], nfail)
        if (ncalls > budget) then
            nfail = nfail + 1
            write (error_unit, "(a,i0)") "FAIL B: evaluation budget exceeded: ", ncalls
        end if
    end subroutine case_noisy_options

    subroutine case_derivative_callback(nfail)
        integer, intent(inout) :: nfail
        complex(dp), allocatable :: roots(:), fvals(:)
        integer, allocatable :: mult(:)
        integer :: nfound, calls_df, calls_fd
        type(fortnum_status_t) :: status

        noisy = .false.
        ncalls = 0
        call complex_region_roots(f_sampled, (-1.0_dp, -1.0_dp), (1.0_dp, 1.0_dp), &
            roots, fvals, mult, nfound, status, newtonf=1.0e-13_dp, &
            contour_tol=1.0e-8_dp, df=df_exact)
        calls_df = ncalls
        call check_zero_set("C exact+df", roots, mult, nfound, status, &
            [1.0e-10_dp, 1.0e-10_dp, 1.0e-7_dp, 1.0e-10_dp], nfail)

        ncalls = 0
        call complex_region_roots(f_sampled, (-1.0_dp, -1.0_dp), (1.0_dp, 1.0_dp), &
            roots, fvals, mult, nfound, status, newtonf=1.0e-13_dp, &
            contour_tol=1.0e-8_dp)
        calls_fd = ncalls
        call check_zero_set("C exact+fd", roots, mult, nfound, status, &
            [1.0e-10_dp, 1.0e-10_dp, 1.0e-7_dp, 1.0e-10_dp], nfail)
        write (*, "(a,i0,a,i0)") "C calls with df=", calls_df, " with fd=", calls_fd
        if (calls_df >= calls_fd) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL C: derivative callback did not save evaluations"
        end if
    end subroutine case_derivative_callback

    subroutine case_zero_free(nfail)
        integer, intent(inout) :: nfail
        complex(dp), allocatable :: roots(:), fvals(:)
        integer, allocatable :: mult(:)
        integer :: nfound, calls_free
        type(fortnum_status_t) :: status

        noisy = .true.
        ncalls = 0
        call complex_region_roots(f_sampled, (-2.0_dp, 1.0_dp), (-1.0_dp, 2.0_dp), &
            roots, fvals, mult, nfound, status, newtonf=1.0e-5_dp, &
            contour_tol=1.0e-4_dp, fd_step=1.0e-3_dp)
        calls_free = ncalls
        write (*, "(a,i0,a,i0,2a)") "D zero-free: nfound=", nfound, " calls=", calls_free, &
            " status: ", trim(status%msg)
        if (status%code /= FORTNUM_OK .or. nfound /= 0) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL D: zero-free box reported zeros or failed"
        end if
        if (calls_free > budget/10) then
            nfail = nfail + 1
            write (error_unit, "(a,i0)") "FAIL D: zero-free box too expensive: ", calls_free
        end if
    end subroutine case_zero_free

    ! Exactly the ninside constructed zeros, each once, with its
    ! multiplicity, within tol(k) of zexact(k); nothing else.
    subroutine check_zero_set(label, roots, mult, nfound, status, tol, nfail)
        character(len=*), intent(in) :: label
        complex(dp), intent(in) :: roots(:)
        integer, intent(in) :: mult(:), nfound
        type(fortnum_status_t), intent(in) :: status
        real(dp), intent(in) :: tol(ninside)
        integer, intent(inout) :: nfail
        integer :: j, k, hits, khit

        write (*, "(2a,i0,a,i0,2a)") label, ": nfound=", nfound, " calls=", ncalls, &
            " status: ", trim(status%msg)
        do k = 1, nfound
            write (*, "(2x,2es24.15,a,i0)") real(roots(k)), aimag(roots(k)), &
                "  mult ", mult(k)
        end do
        if (status%code /= FORTNUM_OK) then
            nfail = nfail + 1
            write (error_unit, "(3a)") "FAIL ", label, ": status not OK"
        end if
        if (nfound /= ninside) then
            nfail = nfail + 1
            write (error_unit, "(3a)") "FAIL ", label, ": wrong distinct-zero count"
        end if
        do j = 1, ninside
            hits = 0
            khit = 0
            do k = 1, nfound
                if (abs(roots(k) - zexact(j)) <= tol(j)) then
                    hits = hits + 1
                    khit = k
                end if
            end do
            if (hits /= 1) then
                nfail = nfail + 1
                write (error_unit, "(3a,i0)") "FAIL ", label, ": zero not found once: ", j
            else if (mult(khit) /= mexact(j)) then
                nfail = nfail + 1
                write (error_unit, "(3a,i0)") "FAIL ", label, ": wrong multiplicity: ", j
            end if
        end do
    end subroutine check_zero_set

end program test_fortnum_region_roots_noisy
