program test_fft_ad
    ! Adjoint-consistency gate for the fft derivative products. Wraps the
    ! complex transforms as real vector maps (interleaved re/im) so the shared
    ! harness can FD-check the JVP and enforce the dot-product identity
    ! u.(Jv)=v.(J^T u) to machine precision.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_fft, only: fft_c2c, fft_c2c_jvp, fft_c2c_vjp, &
        fft_r2c, fft_r2c_jvp, fft_r2c_vjp
    use fortnum_ad_test_utils, only: check_jvp_vs_fd, dot_product_identity
    implicit none

    integer, parameter :: NC = 6   ! c2c radix length
    integer, parameter :: NB = 7   ! c2c Bluestein length
    integer, parameter :: NR = 8   ! r2c length

    integer :: nfail
    nfail = 0

    call test_c2c_jvp_fd(nfail)
    call test_c2c_adjoint(nfail, -1)
    call test_c2c_adjoint(nfail, +1)
    call test_c2c_adjoint_n(nfail, 7)    ! Bluestein length
    call test_r2c_jvp_fd(nfail)
    call test_r2c_adjoint(nfail, 8)
    call test_r2c_adjoint(nfail, 6)
    call test_r2c_adjoint(nfail, 5)      ! odd length

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    ! ---- c2c forward (sign=-1) as a real map R^{2n} -> R^{2n} ----------------

    subroutine c2c_prim(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        complex(dp) :: z(NC)
        call pack_c(x, z)
        call fft_c2c(z, -1)
        call unpack_c(z, y)
    end subroutine c2c_prim

    subroutine c2c_jvp(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        complex(dp) :: dz(NC)
        call pack_c(v, dz)
        call fft_c2c_jvp(dz, -1)
        call unpack_c(dz, jv)
    end subroutine c2c_jvp

    subroutine c2c_vjp(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        complex(dp) :: w(NC)
        call pack_c(u, w)
        call fft_c2c_vjp(w, -1)
        call unpack_c(w, jtu)
    end subroutine c2c_vjp

    subroutine test_c2c_jvp_fd(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x(2*NC), v(2*NC)
        integer :: i
        do i = 1, 2*NC
            x(i) = 0.3_dp*i - 1.1_dp
            v(i) = sin(0.7_dp*i)
        end do
        if (.not. check_jvp_vs_fd("c2c", c2c_prim, c2c_jvp, x, v, &
                tol=1.0e-7_dp)) nfail = nfail + 1
    end subroutine test_c2c_jvp_fd

    subroutine test_c2c_adjoint(nfail, sgn)
        integer, intent(inout) :: nfail
        integer, intent(in) :: sgn
        real(dp) :: x(2*NC), u(2*NC), v(2*NC)
        integer :: i
        do i = 1, 2*NC
            x(i) = 0.0_dp
            u(i) = cos(0.4_dp*i + 0.2_dp)
            v(i) = 0.9_dp*sin(1.3_dp*i)
        end do
        if (sgn == -1) then
            if (.not. dot_product_identity("c2c-adj-fwd", c2c_jvp, c2c_vjp, &
                    x, u, v, tol=1.0e-12_dp)) nfail = nfail + 1
        else
            if (.not. dot_product_identity("c2c-adj-inv", c2c_jvp_p, &
                    c2c_vjp_p, x, u, v, tol=1.0e-12_dp)) nfail = nfail + 1
        end if
    end subroutine test_c2c_adjoint

    ! sign=+1 variants for the adjoint check
    subroutine c2c_jvp_p(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        complex(dp) :: dz(NC)
        call pack_c(v, dz)
        call fft_c2c_jvp(dz, +1)
        call unpack_c(dz, jv)
    end subroutine c2c_jvp_p

    subroutine c2c_vjp_p(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        complex(dp) :: w(NC)
        call pack_c(u, w)
        call fft_c2c_vjp(w, +1)
        call unpack_c(w, jtu)
    end subroutine c2c_vjp_p

    ! ---- Bluestein-length c2c adjoint (n=7) ----------------------------------

    subroutine c2c_jvp_b(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        complex(dp) :: dz(NB)
        call pack_n(v, dz, NB)
        call fft_c2c_jvp(dz, -1)
        call unpack_n(dz, jv, NB)
    end subroutine c2c_jvp_b

    subroutine c2c_vjp_b(x, u, jtu)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: u(:)
        real(dp), intent(out) :: jtu(:)
        complex(dp) :: w(NB)
        call pack_n(u, w, NB)
        call fft_c2c_vjp(w, -1)
        call unpack_n(w, jtu, NB)
    end subroutine c2c_vjp_b

    subroutine test_c2c_adjoint_n(nfail, n)
        integer, intent(inout) :: nfail
        integer, intent(in) :: n
        real(dp) :: x(2*n), u(2*n), v(2*n)
        integer :: i
        do i = 1, 2*n
            x(i) = 0.0_dp
            u(i) = cos(0.4_dp*i + 0.2_dp)
            v(i) = 0.9_dp*sin(1.3_dp*i)
        end do
        if (.not. dot_product_identity("c2c-adj-bluestein", c2c_jvp_b, &
                c2c_vjp_b, x, u, v, tol=1.0e-12_dp)) nfail = nfail + 1
    end subroutine test_c2c_adjoint_n

    ! ---- r2c as a real map R^n -> R^{2*(n/2+1)} ------------------------------

    subroutine r2c_prim(x, y)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(out) :: y(:)
        complex(dp) :: c(NR/2 + 1)
        call fft_r2c(x, c)
        call unpack_c(c, y)
    end subroutine r2c_prim

    subroutine r2c_jvp(x, v, jv)
        real(dp), intent(in)  :: x(:)
        real(dp), intent(in)  :: v(:)
        real(dp), intent(out) :: jv(:)
        complex(dp) :: dc(NR/2 + 1)
        call fft_r2c_jvp(v, dc)
        call unpack_c(dc, jv)
    end subroutine r2c_jvp

    subroutine test_r2c_jvp_fd(nfail)
        integer, intent(inout) :: nfail
        real(dp) :: x(NR), v(NR), jv((NR/2 + 1)*2)
        integer :: i
        do i = 1, NR
            x(i) = 0.5_dp*i - 1.7_dp
            v(i) = cos(0.6_dp*i)
        end do
        ! check_jvp_vs_fd assumes a square map (m=n); the r2c map is
        ! rectangular, so call fd directly via a square embedding is awkward.
        ! Instead verify JVP against central FD component-wise here.
        call fd_check_rect(nfail, "r2c", r2c_prim, r2c_jvp, x, v, &
            (NR/2 + 1)*2)
    end subroutine test_r2c_jvp_fd

    ! component-wise central-FD check for a rectangular map
    subroutine fd_check_rect(nfail, label, f, jvp, x, v, m)
        integer, intent(inout) :: nfail
        character(*), intent(in) :: label
        interface
            subroutine f(x, y)
                import :: dp
                real(dp), intent(in)  :: x(:)
                real(dp), intent(out) :: y(:)
            end subroutine f
            subroutine jvp(x, v, jv)
                import :: dp
                real(dp), intent(in)  :: x(:)
                real(dp), intent(in)  :: v(:)
                real(dp), intent(out) :: jv(:)
            end subroutine jvp
        end interface
        real(dp), intent(in) :: x(:), v(:)
        integer, intent(in) :: m
        real(dp), allocatable :: yp(:), ym(:), jvfd(:), jvad(:)
        real(dp) :: h, worst, e
        integer :: i
        h = 1.0e-6_dp
        allocate (yp(m), ym(m), jvfd(m), jvad(m))
        call f(x + h*v, yp)
        call f(x - h*v, ym)
        jvfd = (yp - ym)/(2.0_dp*h)
        call jvp(x, v, jvad)
        worst = 0.0_dp
        do i = 1, m
            e = abs(jvad(i) - jvfd(i))/max(abs(jvfd(i)), 1.0_dp)
            worst = max(worst, e)
        end do
        if (worst > 1.0e-7_dp) then
            write (error_unit, '(a,a,a,es12.4)') "FAIL [", label, &
                "] r2c jvp-vs-fd worst=", worst
            nfail = nfail + 1
        end if
    end subroutine fd_check_rect

    ! r2c adjoint identity for various lengths
    subroutine test_r2c_adjoint(nfail, n)
        integer, intent(inout) :: nfail
        integer, intent(in) :: n
        real(dp) :: u(2*(n/2 + 1)), v(n), jv(2*(n/2 + 1)), jtu(n)
        real(dp) :: lhs, rhs, e
        complex(dp) :: dc(n/2 + 1), uc(n/2 + 1)
        integer :: i
        do i = 1, n
            v(i) = sin(0.8_dp*i + 0.1_dp)
        end do
        do i = 1, 2*(n/2 + 1)
            u(i) = cos(0.5_dp*i - 0.3_dp)
        end do
        ! Jv = fft_r2c_jvp(v)
        call fft_r2c_jvp(v, dc)
        call unpack_c(dc, jv)
        ! J^T u = fft_r2c_vjp(u)
        call pack_c(u, uc)
        call fft_r2c_vjp(uc, jtu)
        lhs = dot_product(u, jv)
        rhs = dot_product(v, jtu)
        e = abs(lhs - rhs)/max(abs(lhs), abs(rhs), 1.0_dp)
        if (e > 1.0e-12_dp) then
            write (error_unit, '(a,i0,a,es24.16,a,es24.16,a,es12.4)') &
                "FAIL [r2c-adj n=", n, "] u.(Jv)=", lhs, " v.(J^T u)=", rhs, &
                " rel_err=", e
            nfail = nfail + 1
        end if
    end subroutine test_r2c_adjoint

    ! ---- interleaving helpers ------------------------------------------------

    subroutine pack_c(r, z)
        real(dp), intent(in)  :: r(:)
        complex(dp), intent(out) :: z(:)
        integer :: k
        do k = 1, size(z)
            z(k) = cmplx(r(2*k - 1), r(2*k), dp)
        end do
    end subroutine pack_c

    subroutine unpack_c(z, r)
        complex(dp), intent(in) :: z(:)
        real(dp), intent(out) :: r(:)
        integer :: k
        do k = 1, size(z)
            r(2*k - 1) = real(z(k), dp)
            r(2*k) = aimag(z(k))
        end do
    end subroutine unpack_c

    subroutine pack_n(r, z, n)
        real(dp), intent(in)  :: r(:)
        complex(dp), intent(out) :: z(:)
        integer, intent(in) :: n
        integer :: k
        do k = 1, n
            z(k) = cmplx(r(2*k - 1), r(2*k), dp)
        end do
    end subroutine pack_n

    subroutine unpack_n(z, r, n)
        complex(dp), intent(in) :: z(:)
        real(dp), intent(out) :: r(:)
        integer, intent(in) :: n
        integer :: k
        do k = 1, n
            r(2*k - 1) = real(z(k), dp)
            r(2*k) = aimag(z(k))
        end do
    end subroutine unpack_n

end program test_fft_ad
