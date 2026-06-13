module fortnum_fft
    ! 1D discrete Fourier transforms.
    !
    ! Convention (FFTW / numpy.fft): the forward transform is
    !   c(k+1) = sum_{j=0}^{n-1} x(j+1) exp(-2 pi i j k / n),  unnormalized.
    ! The inverse c2c applies exp(+2 pi i j k / n) without the 1/n factor;
    ! the caller normalizes when needed.
    ! Real input uses the half-length complex packing of Cooley, Lewis and
    ! Welch, J. Sound Vib. 12 (1970) 315.
    !
    ! Algorithm: mixed-radix (2,3,4,5) self-sorting Stockham passes after
    ! Temperton, J. Comput. Phys. 52 (1983) 1; Bluestein chirp-z fallback for
    ! all other lengths (Bluestein, IEEE Trans. Audio Electroacoust. 18 (1970)
    ! 451).
    !
    ! Derivative policy (ad.md §1, §4):
    !   Default class: transparent.
    !   The DFT is linear; Enzyme differentiates the implementation directly.
    !   Active arguments: z (fft_c2c), x and c (fft_r2c).
    !   Inactive arguments: sign, plan, all integers.
    !   Derivative entry points (ad.md §2): fft_c2c_jvp, fft_c2c_vjp,
    !   fft_r2c_jvp, fft_r2c_vjp — none implemented yet (M1 primal only).
    !
    ! No module-level state. Pass a caller-owned fortnum_fft_plan_t to reuse
    ! twiddle tables across transforms of the same length (thread safe: plans
    ! are read-only during the transform).

    use, intrinsic :: iso_fortran_env, only: int64
    use fortnum_kinds, only: dp
    implicit none
    private

    public :: fft_r2c, fft_c2c, fortnum_fft_plan_t, fft_plan_init
    public :: fft_c2c_jvp, fft_c2c_vjp, fft_r2c_jvp, fft_r2c_vjp

    real(dp), parameter :: pi = 3.141592653589793238462643383279502884_dp
    real(dp), parameter :: sin60 = 0.866025403784438646763723170752936183_dp
    real(dp), parameter :: cos72 = 0.309016994374947424102293417182819059_dp
    real(dp), parameter :: sin72 = 0.951056516295153572116439333379382143_dp
    real(dp), parameter :: cos144 = -0.809016994374947424102293417182819059_dp
    real(dp), parameter :: sin144 = 0.587785252292473129168705954639072769_dp
    integer, parameter :: max_fac = 40

    type :: fortnum_fft_plan_t
        integer :: n = 0
        integer :: ncpx = 0
        integer :: nfac = 0
        integer :: m = 0
        integer :: nfacm = 0
        logical :: even = .false.
        logical :: bluestein = .false.
        integer :: fac(max_fac) = 0
        integer :: facm(max_fac) = 0
        complex(dp), allocatable :: w(:)
        complex(dp), allocatable :: wr(:)
        complex(dp), allocatable :: wm(:)
        complex(dp), allocatable :: chirp(:)
        complex(dp), allocatable :: bfft(:)
    end type fortnum_fft_plan_t

contains

    subroutine fft_plan_init(plan, n)
        ! Initialize a reusable plan for transforms of length n.
        type(fortnum_fft_plan_t), intent(out) :: plan
        integer, intent(in) :: n
        integer :: k

        if (n < 1) error stop "fft_plan_init: n must be positive"
        plan%n = n
        plan%even = mod(n, 2) == 0
        if (plan%even) then
            plan%ncpx = n/2
            allocate (plan%wr(0:n/2))
            do k = 0, n/2
                plan%wr(k) = twiddle(k, n)
            end do
        else
            plan%ncpx = n
        end if
        call init_core(plan, plan%ncpx)
    end subroutine fft_plan_init

    subroutine fft_r2c(x, c, plan)
        ! Real-to-complex DFT. c has size n/2+1; c(k+1) is the k-th bin.
        ! If plan is supplied it must have been built for size(x).
        real(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: c(:)
        type(fortnum_fft_plan_t), intent(in), optional :: plan
        type(fortnum_fft_plan_t) :: local_plan
        integer :: n

        n = size(x)
        if (n < 1) error stop "fft_r2c: empty input"
        if (size(c) /= n/2 + 1) error stop "fft_r2c: size(c) /= size(x)/2 + 1"
        if (present(plan)) then
            if (plan%n /= n) error stop "fft_r2c: plan built for different length"
            call r2c_planned(plan, x, c)
        else
            call fft_plan_init(local_plan, n)
            call r2c_planned(local_plan, x, c)
        end if
    end subroutine fft_r2c

    subroutine fft_c2c(z, sign)
        ! In-place complex-to-complex DFT.
        ! sign = -1: forward (exp(-2 pi i j k / n), FFTW / numpy convention).
        ! sign = +1: inverse (exp(+2 pi i j k / n)), unnormalized.
        complex(dp), intent(inout) :: z(:)
        integer, intent(in) :: sign
        type(fortnum_fft_plan_t) :: plan

        if (size(z) < 1) error stop "fft_c2c: empty input"
        if (abs(sign) /= 1) error stop "fft_c2c: sign must be -1 or +1"
        plan%n = size(z)
        plan%ncpx = size(z)
        call init_core(plan, plan%ncpx)
        ! Inverse via conj(forward(conj(.))): conjugate, forward, conjugate.
        if (sign == 1) z = conjg(z)
        call c2c_forward(plan, z)
        if (sign == 1) z = conjg(z)
    end subroutine fft_c2c

    ! --- derivative products (ad.md sec 2, sec 4: transparent/linear) ---------
    !
    ! The DFT is C-linear: c = F z with F(k,j) = exp(sign * 2 pi i j k / n).
    ! Hence J = F and the forward product (JVP) is the transform itself applied
    ! to the tangent. The reverse product (VJP) is the real-adjoint of z -> F z
    ! under the real inner product <a,b> = Re(sum conjg(a) b) on C^n ~ R^2n: it
    ! is multiplication by F^H, the conjugate-transpose. For F(k,j) above,
    ! F^H(j,k) = conjg(F(k,j)) = exp(-sign * 2 pi i j k / n), i.e. fft_c2c with
    ! the opposite sign. The dot-product identity u.(Jv)=v.(J^T u) then holds to
    ! machine precision because both sides equal Re(u^H F v).
    !
    ! Active argument: z (the tangent / cotangent). Inactive: sign, all sizes
    ! (ad.md sec 3). HVP is deferred: the map is linear, so the Hessian is zero
    ! and fft_*_hvp would carry no information.

    subroutine fft_c2c_jvp(dz, sign)
        ! Forward product J(z) dz for the in-place complex DFT. J = F is the
        ! transform, so the JVP is fft_c2c applied to the tangent dz.
        complex(dp), intent(inout) :: dz(:)
        integer, intent(in) :: sign
        call fft_c2c(dz, sign)
    end subroutine fft_c2c_jvp

    subroutine fft_c2c_vjp(u, sign)
        ! Reverse product J(z)^H u = F^H u. F^H is the conjugate-transpose, the
        ! transform with the opposite sign (conjugate kernel, unnormalized).
        complex(dp), intent(inout) :: u(:)
        integer, intent(in) :: sign
        call fft_c2c(u, -sign)
    end subroutine fft_c2c_vjp

    subroutine fft_r2c_jvp(dx, dc, plan)
        ! Forward product for the real-to-complex DFT. The map x -> c is
        ! R-linear with the same coefficients as the primal, so the JVP is
        ! fft_r2c applied to the tangent dx.
        real(dp), intent(in) :: dx(:)
        complex(dp), intent(out) :: dc(:)
        type(fortnum_fft_plan_t), intent(in), optional :: plan
        if (present(plan)) then
            call fft_r2c(dx, dc, plan)
        else
            call fft_r2c(dx, dc)
        end if
    end subroutine fft_r2c_jvp

    subroutine fft_r2c_vjp(u, xbar, plan)
        ! Reverse product for the real-to-complex DFT. The primal is the real
        ! map (R x)(k) = sum_j x(j) exp(-2 pi i j k / n), k = 0..n/2, with x
        ! real. Under <p,q>_R on the complex output and the ordinary real inner
        ! product on x, the adjoint is xbar(j) = Re( sum_k conjg(R(k,j)) u(k) ),
        ! i.e. xbar(j) = Re( sum_{k=0}^{n/2} exp(+2 pi i j k / n) u(k) ).
        ! Computed by Hermitian-extending u to a full length-n spectrum and
        ! taking one length-n inverse-direction transform; the real part is the
        ! adjoint. xbar is sized n = size(x); u is sized n/2+1.
        complex(dp), intent(in) :: u(:)
        real(dp), intent(out) :: xbar(:)
        type(fortnum_fft_plan_t), intent(in), optional :: plan
        complex(dp), allocatable :: full(:)
        integer :: n, nh, k

        n = size(xbar)
        if (n < 1) error stop "fft_r2c_vjp: empty output"
        if (size(u) /= n/2 + 1) error stop "fft_r2c_vjp: size(u) /= n/2 + 1"
        if (present(plan)) then
            if (plan%n /= n) error stop &
                "fft_r2c_vjp: plan built for different length"
        end if

        ! xbar(j) = Re( sum_{k=0}^{n/2} exp(+i 2 pi j k / n) u(k) ). Build the
        ! length-n vector g(k)=u(k) for k=0..n/2, g(k)=0 otherwise, then apply
        ! the sign=+1 (conjugate-kernel) transform; Re of the result is xbar.
        allocate (full(n))
        full = (0.0_dp, 0.0_dp)
        do k = 0, n/2
            full(k + 1) = u(k + 1)
        end do
        call fft_c2c(full, 1)
        do k = 1, n
            xbar(k) = real(full(k), dp)
        end do
    end subroutine fft_r2c_vjp

    ! --- private implementation below ----------------------------------------

    subroutine r2c_planned(p, x, c)
        type(fortnum_fft_plan_t), intent(in) :: p
        real(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: c(:)

        if (p%even) then
            call r2c_even(p, x, c)
        else
            call r2c_odd(p, x, c)
        end if
    end subroutine r2c_planned

    subroutine r2c_even(p, x, c)
        type(fortnum_fft_plan_t), intent(in) :: p
        real(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: c(:)
        complex(dp) :: z(0:p%ncpx - 1)
        complex(dp) :: zk, zn, even_part, odd_part
        integer :: nh, k

        nh = p%ncpx
        do k = 0, nh - 1
            z(k) = cmplx(x(2*k + 1), x(2*k + 2), dp)
        end do
        call c2c_forward(p, z)
        c(1) = cmplx(real(z(0)) + aimag(z(0)), 0.0_dp, dp)
        c(nh + 1) = cmplx(real(z(0)) - aimag(z(0)), 0.0_dp, dp)
        do k = 1, nh - 1
            zk = z(k)
            zn = conjg(z(nh - k))
            even_part = 0.5_dp*(zk + zn)
            odd_part = cmplx(0.0_dp, -0.5_dp, dp)*(zk - zn)
            c(k + 1) = even_part + p%wr(k)*odd_part
        end do
    end subroutine r2c_even

    subroutine r2c_odd(p, x, c)
        type(fortnum_fft_plan_t), intent(in) :: p
        real(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: c(:)
        complex(dp) :: z(0:p%ncpx - 1)
        integer :: k

        do k = 0, p%n - 1
            z(k) = cmplx(x(k + 1), 0.0_dp, dp)
        end do
        call c2c_forward(p, z)
        do k = 0, p%n/2
            c(k + 1) = z(k)
        end do
    end subroutine r2c_odd

    subroutine c2c_forward(p, z)
        type(fortnum_fft_plan_t), intent(in) :: p
        complex(dp), intent(inout) :: z(0:)
        complex(dp) :: work(0:p%ncpx - 1)

        if (p%bluestein) then
            call bluestein_forward(p, z)
        else
            call stockham_forward(p%ncpx, p%nfac, p%fac, p%w, z, work)
        end if
    end subroutine c2c_forward

    subroutine init_core(plan, nc)
        type(fortnum_fft_plan_t), intent(inout) :: plan
        integer, intent(in) :: nc
        integer :: rem, k

        call factorize(nc, plan%fac, plan%nfac, rem)
        plan%bluestein = rem /= 1
        if (plan%bluestein) then
            call init_bluestein(plan, nc)
        else
            allocate (plan%w(0:nc - 1))
            do k = 0, nc - 1
                plan%w(k) = twiddle(k, nc)
            end do
        end if
    end subroutine init_core

    subroutine init_bluestein(plan, nc)
        type(fortnum_fft_plan_t), intent(inout) :: plan
        integer, intent(in) :: nc
        complex(dp), allocatable :: work(:)
        integer :: m, rem, k

        m = 1
        do while (m < 2*nc - 1)
            m = 2*m
        end do
        plan%m = m
        call factorize(m, plan%facm, plan%nfacm, rem)
        allocate (plan%wm(0:m - 1), plan%chirp(0:nc - 1), plan%bfft(0:m - 1))
        allocate (work(0:m - 1))
        do k = 0, m - 1
            plan%wm(k) = twiddle(k, m)
        end do
        do k = 0, nc - 1
            plan%chirp(k) = chirp_factor(k, nc)
        end do
        plan%bfft = (0.0_dp, 0.0_dp)
        plan%bfft(0:nc - 1) = conjg(plan%chirp)
        do k = 1, nc - 1
            plan%bfft(m - k) = conjg(plan%chirp(k))
        end do
        call stockham_forward(m, plan%nfacm, plan%facm, plan%wm, plan%bfft, work)
        ! absorb 1/m of the length-m inverse into the Bluestein filter
        plan%bfft = plan%bfft/real(m, dp)
    end subroutine init_bluestein

    subroutine bluestein_forward(p, z)
        type(fortnum_fft_plan_t), intent(in) :: p
        complex(dp), intent(inout) :: z(0:)
        complex(dp), allocatable :: a(:), work(:)
        integer :: nc, m, k

        nc = p%ncpx
        m = p%m
        allocate (a(0:m - 1), work(0:m - 1))
        do k = 0, nc - 1
            a(k) = z(k)*p%chirp(k)
        end do
        a(nc:m - 1) = (0.0_dp, 0.0_dp)
        call stockham_forward(m, p%nfacm, p%facm, p%wm, a, work)
        ! convolve: inverse via conj(forward(conj(.))), with 1/m inside bfft
        a = conjg(a*p%bfft)
        call stockham_forward(m, p%nfacm, p%facm, p%wm, a, work)
        do k = 0, nc - 1
            z(k) = p%chirp(k)*conjg(a(k))
        end do
    end subroutine bluestein_forward

    subroutine stockham_forward(nc, nfac, fac, w, a, work)
        integer, intent(in) :: nc, nfac
        integer, intent(in) :: fac(:)
        complex(dp), intent(in) :: w(0:nc - 1)
        complex(dp), intent(inout) :: a(0:nc - 1)
        complex(dp), intent(inout) :: work(0:nc - 1)
        integer :: la, lrem, kf
        logical :: in_a

        la = 1
        in_a = .true.
        do kf = 1, nfac
            lrem = nc/(la*fac(kf))
            if (in_a) then
                call radix_pass(fac(kf), la, lrem, nc, w, a, work)
            else
                call radix_pass(fac(kf), la, lrem, nc, w, work, a)
            end if
            in_a = .not. in_a
            la = la*fac(kf)
        end do
        if (.not. in_a) a = work
    end subroutine stockham_forward

    subroutine radix_pass(r, la, lrem, nc, w, src, dst)
        integer, intent(in) :: r, la, lrem, nc
        complex(dp), intent(in) :: w(0:nc - 1), src(0:nc - 1)
        complex(dp), intent(out) :: dst(0:nc - 1)

        select case (r)
        case (4)
            call pass4(la, lrem, nc, w, src, dst)
        case (2)
            call pass2(la, lrem, nc, w, src, dst)
        case (3)
            call pass3(la, lrem, nc, w, src, dst)
        case (5)
            call pass5(la, lrem, nc, w, src, dst)
        case default
            error stop "fortnum_fft: unsupported radix"
        end select
    end subroutine radix_pass

    subroutine pass2(la, lrem, nc, w, src, dst)
        integer, intent(in) :: la, lrem, nc
        complex(dp), intent(in) :: w(0:nc - 1), src(0:nc - 1)
        complex(dp), intent(out) :: dst(0:nc - 1)
        complex(dp) :: u0, u1
        integer :: q, k, i0, i1, o0

        if (la == 1) then
            do q = 0, lrem - 1
                u0 = src(q)
                u1 = src(q + lrem)
                dst(2*q) = u0 + u1
                dst(2*q + 1) = u0 - u1
            end do
            return
        end if
        do q = 0, lrem - 1
            i0 = q*la
            i1 = i0 + lrem*la
            o0 = 2*q*la
            do k = 0, la - 1
                u0 = src(i0 + k)
                u1 = w(k*lrem)*src(i1 + k)
                dst(o0 + k) = u0 + u1
                dst(o0 + la + k) = u0 - u1
            end do
        end do
    end subroutine pass2

    subroutine pass3(la, lrem, nc, w, src, dst)
        integer, intent(in) :: la, lrem, nc
        complex(dp), intent(in) :: w(0:nc - 1), src(0:nc - 1)
        complex(dp), intent(out) :: dst(0:nc - 1)
        complex(dp) :: u0, u1, u2, t1, t2, t3
        integer :: q, k, i0, i1, i2, o0

        do q = 0, lrem - 1
            i0 = q*la
            i1 = i0 + lrem*la
            i2 = i1 + lrem*la
            o0 = 3*q*la
            do k = 0, la - 1
                u0 = src(i0 + k)
                u1 = w(k*lrem)*src(i1 + k)
                u2 = w(2*k*lrem)*src(i2 + k)
                t1 = u1 + u2
                t2 = u0 - 0.5_dp*t1
                t3 = sin60*(u1 - u2)
                t3 = cmplx(aimag(t3), -real(t3), dp)
                dst(o0 + k) = u0 + t1
                dst(o0 + la + k) = t2 + t3
                dst(o0 + 2*la + k) = t2 - t3
            end do
        end do
    end subroutine pass3

    subroutine pass4(la, lrem, nc, w, src, dst)
        integer, intent(in) :: la, lrem, nc
        complex(dp), intent(in) :: w(0:nc - 1), src(0:nc - 1)
        complex(dp), intent(out) :: dst(0:nc - 1)
        complex(dp) :: u0, u1, u2, u3, s0, s1, d0, d1
        integer :: q, k, i0, i1, i2, i3, o0

        if (la == 1) then
            do q = 0, lrem - 1
                u0 = src(q)
                u1 = src(q + lrem)
                u2 = src(q + 2*lrem)
                u3 = src(q + 3*lrem)
                s0 = u0 + u2
                d0 = u0 - u2
                s1 = u1 + u3
                d1 = u1 - u3
                d1 = cmplx(aimag(d1), -real(d1), dp)
                dst(4*q) = s0 + s1
                dst(4*q + 1) = d0 + d1
                dst(4*q + 2) = s0 - s1
                dst(4*q + 3) = d0 - d1
            end do
            return
        end if
        do q = 0, lrem - 1
            i0 = q*la
            i1 = i0 + lrem*la
            i2 = i1 + lrem*la
            i3 = i2 + lrem*la
            o0 = 4*q*la
            do k = 0, la - 1
                u0 = src(i0 + k)
                u1 = w(k*lrem)*src(i1 + k)
                u2 = w(2*k*lrem)*src(i2 + k)
                u3 = w(3*k*lrem)*src(i3 + k)
                s0 = u0 + u2
                d0 = u0 - u2
                s1 = u1 + u3
                d1 = u1 - u3
                d1 = cmplx(aimag(d1), -real(d1), dp)
                dst(o0 + k) = s0 + s1
                dst(o0 + la + k) = d0 + d1
                dst(o0 + 2*la + k) = s0 - s1
                dst(o0 + 3*la + k) = d0 - d1
            end do
        end do
    end subroutine pass4

    subroutine pass5(la, lrem, nc, w, src, dst)
        integer, intent(in) :: la, lrem, nc
        complex(dp), intent(in) :: w(0:nc - 1), src(0:nc - 1)
        complex(dp), intent(out) :: dst(0:nc - 1)
        complex(dp) :: u0, u1, u2, u3, u4
        integer :: q, k, i0, i1, i2, i3, i4, o0

        if (la == 1) then
            do q = 0, lrem - 1
                call butterfly5(src(q), src(q + lrem), src(q + 2*lrem), &
                                src(q + 3*lrem), src(q + 4*lrem), &
                                dst(5*q), dst(5*q + 1), dst(5*q + 2), &
                                dst(5*q + 3), dst(5*q + 4))
            end do
            return
        end if
        do q = 0, lrem - 1
            i0 = q*la
            i1 = i0 + lrem*la
            i2 = i1 + lrem*la
            i3 = i2 + lrem*la
            i4 = i3 + lrem*la
            o0 = 5*q*la
            do k = 0, la - 1
                u0 = src(i0 + k)
                u1 = w(k*lrem)*src(i1 + k)
                u2 = w(2*k*lrem)*src(i2 + k)
                u3 = w(3*k*lrem)*src(i3 + k)
                u4 = w(4*k*lrem)*src(i4 + k)
                call butterfly5(u0, u1, u2, u3, u4, dst(o0 + k), &
                                dst(o0 + la + k), dst(o0 + 2*la + k), &
                                dst(o0 + 3*la + k), dst(o0 + 4*la + k))
            end do
        end do
    end subroutine pass5

    pure subroutine butterfly5(u0, u1, u2, u3, u4, d0, d1, d2, d3, d4)
        complex(dp), intent(in) :: u0, u1, u2, u3, u4
        complex(dp), intent(out) :: d0, d1, d2, d3, d4
        complex(dp) :: t1, t2, t3, t4, m1, m2, m3, m4

        t1 = u1 + u4
        t2 = u2 + u3
        t3 = u1 - u4
        t4 = u2 - u3
        m1 = u0 + cos72*t1 + cos144*t2
        m2 = u0 + cos144*t1 + cos72*t2
        m3 = sin72*t3 + sin144*t4
        m3 = cmplx(aimag(m3), -real(m3), dp)
        m4 = sin144*t3 - sin72*t4
        m4 = cmplx(aimag(m4), -real(m4), dp)
        d0 = u0 + t1 + t2
        d1 = m1 + m3
        d2 = m2 + m4
        d3 = m2 - m4
        d4 = m1 - m3
    end subroutine butterfly5

    subroutine factorize(nc, fac, nfac, rem)
        integer, intent(in) :: nc
        integer, intent(out) :: fac(:)
        integer, intent(out) :: nfac, rem

        fac = 0
        nfac = 0
        rem = nc
        do while (mod(rem, 4) == 0)
            nfac = nfac + 1
            fac(nfac) = 4
            rem = rem/4
        end do
        if (mod(rem, 2) == 0) then
            nfac = nfac + 1
            fac(nfac) = 2
            rem = rem/2
        end if
        do while (mod(rem, 3) == 0)
            nfac = nfac + 1
            fac(nfac) = 3
            rem = rem/3
        end do
        do while (mod(rem, 5) == 0)
            nfac = nfac + 1
            fac(nfac) = 5
            rem = rem/5
        end do
    end subroutine factorize

    pure function twiddle(j, n) result(wj)
        integer, intent(in) :: j, n
        complex(dp) :: wj
        real(dp) :: angle

        angle = -2.0_dp*pi*real(j, dp)/real(n, dp)
        wj = cmplx(cos(angle), sin(angle), dp)
    end function twiddle

    pure function chirp_factor(j, n) result(ch)
        integer, intent(in) :: j, n
        complex(dp) :: ch
        real(dp) :: angle
        integer(int64) :: j2

        ! exp(-i pi j^2 / n); reduce j^2 mod 2n to keep the argument small
        j2 = modulo(int(j, int64)*int(j, int64), 2_int64*int(n, int64))
        angle = -pi*real(j2, dp)/real(n, dp)
        ch = cmplx(cos(angle), sin(angle), dp)
    end function chirp_factor

end module fortnum_fft
