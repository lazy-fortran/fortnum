!> Ball Fourier series against real128 oracles. A series (balls a_n, tail
!> tau) represents every function sum c_n e^{i n theta} with
!> sum_n excess_n <= tau, excess_n = max(0, abs(c_n - mid a_n) - rad a_n)
!> inside the kept window and abs(c_n) outside (weighted by
!> exp(sigma (abs(m) + abs(n))) in 2D). The tests draw concrete members of
!> the input series, form the exact real128 result, and check membership of
!> the result series.
program test_fortnum_fseries
    use, intrinsic :: iso_fortran_env, only: dp => real64, qp => real128, &
        error_unit
    use fortnum_ball, only: ball_t, ball
    use fortnum_cfft_rigorous, only: rigorous_fft_plan_t, rigorous_fft_plan_init
    use fortnum_fseries, only: fseries_t, fs_zero, fs_trig, fs_add, fs_sub, &
        fs_scale, fs_l1_hi, fs_mul, fs_mul_fast, fs_dtheta, fs_inner_w, &
        fs_inner_w_fast, fs_inverse, fs_fft_size, fseries2d_t, fs2_zero, &
        fs2_trig, fs2_trig_tail, fs2_add, fs2_sup_hi, fs2_wnorm_hi, fs2_mul, &
        fs2_mul_fast, fs2_dtheta, fs2_dphi, fs2_inner_w, fs2_inner_w_fast, &
        fs2_inverse
    implicit none

    integer, parameter :: kt = 40
    integer :: nfail, it, n, seed_size, j
    integer, allocatable :: seed(:)
    type(fseries_t) :: f, g, h, hf, w, y
    type(fseries2d_t) :: f2, g2, h2, hf2, y2
    type(rigorous_fft_plan_t) :: plan
    type(ball_t) :: s, sf
    complex(qp) :: cf(-kt:kt), cg(-kt:kt), cw(-kt:kt), cp(-2*kt:2*kt), z, v
    complex(qp) :: c2f(-8:8, -6:6), c2g(-8:8, -6:6), c2p(-16:16, -12:12)
    real(qp) :: th, ph
    logical :: ok, ok2
    complex(dp) :: cc(-3:3), c2(-2:2, -1:1)

    nfail = 0
    call random_seed(size=seed_size)
    allocate (seed(seed_size))
    seed = 31415
    call random_seed(put=seed)

    ! 1D products, direct and fast, with radii and tails.
    ok = .true.
    ok2 = .true.
    do it = 1, 40
        f = random_series(12, 1.0e-9_dp, 1.0e-7_dp)
        g = random_series(9, 1.0e-12_dp, 0.0_dp)
        call member(f, 12, cf)
        call member(g, 9, cg)
        cp = (0.0_qp, 0.0_qp)
        call conv1(cf, cg, cp)
        h = fs_mul(f, g, 15)
        hf = fs_mul_fast(f, g, 15)
        if (.not. contains1(h, cp, 2*kt)) ok = .false.
        if (.not. contains1(hf, cp, 2*kt)) ok2 = .false.
    end do
    call require(ok, "fs_mul result contains the exact real128 product", nfail)
    call require(ok2, "fs_mul_fast result contains the exact real128 product", &
        nfail)
    ! Coherent worst case: positive centres, members on the ball edges and a
    ! positive tail, so every error term adds up; a missing propagation term
    ! makes the check fail.
    ok = .true.
    do it = 1, 4
        w = worst_series(10, 1.0e-6_dp, merge(1.0e-5_dp, 0.0_dp, it > 2))
        y = worst_series(7, 1.0e-6_dp, merge(1.0e-5_dp, 0.0_dp, mod(it, 2) == 0))
        call member_edge(w, 10, cf)
        call member_edge(y, 7, cg)
        call conv1(cf, cg, cp)
        if (.not. contains1(fs_mul(w, y, 12), cp, 2*kt)) ok = .false.
        if (.not. contains1(fs_mul_fast(w, y, 12), cp, 2*kt)) ok = .false.
    end do
    call require(ok, "products contain the coherent worst case", nfail)
    call member(f, 12, cf)
    call member(g, 9, cg)
    call conv1(cf, cg, cp)
    call rigorous_fft_plan_init(plan, fs_fft_size(20, 20))
    hf = fs_mul_fast(f, g, 30, plan)
    call require(contains1(hf, cp, 2*kt), "fs_mul_fast with a reused plan", nfail)

    ! Sums, scaling, l1 bound.
    call member(f, 12, cf)
    call member(g, 9, cg)
    h = fs_sub(fs_add(f, fs_scale(g, ball((2.0_dp, -1.0_dp), 0.0_dp))), g)
    cp = (0.0_qp, 0.0_qp)
    cp(-kt:kt) = cf + cg*cmplx(2.0_qp, -1.0_qp, qp) - cg
    call require(contains1(h, cp, 2*kt), "add, sub and scale", nfail)
    call require(real(fs_l1_hi(f), qp) >= sum(abs(cf)), "l1 upper bound", nfail)
    cc = [(1, 1), (0, 2), (3, 0), (1, 0), (0, 0), (0, -1), (2, 2)]
    h = fs_dtheta(fs_trig(cc, 3))
    cp = (0.0_qp, 0.0_qp)
    do n = -3, 3
        cp(n) = cmplx(cc(n), kind=qp)*cmplx(0, n, qp)
    end do
    call require(contains1(h, cp, 2*kt) .and. h%tau <= 0.0_dp, "derivative", nfail)
    h = fs_dtheta(f)
    call require(h%tau >= huge(1.0_dp), "derivative of a tail is unbounded", nfail)

    ! Weighted inner product (1/2pi) int conj(F) G w.
    ok = .true.
    do it = 1, 20
        f = random_series(10, 1.0e-10_dp, 1.0e-8_dp)
        g = random_series(8, 1.0e-10_dp, 1.0e-8_dp)
        w = random_series(5, 0.0_dp, 1.0e-9_dp)
        call member(f, 10, cf)
        call member(g, 8, cg)
        call member(w, 5, cw)
        v = (0.0_qp, 0.0_qp)
        do n = -kt, kt
            do j = -kt, kt
                if (abs(n - j) <= kt) v = v + conjg(cf(n))*cg(j)*cw(n - j)
            end do
        end do
        s = fs_inner_w(f, g, w, fs_l1_hi(w))
        sf = fs_inner_w_fast(f, g, w, fs_l1_hi(w), 20)
        if (.not. (abs(v - cmplx(s%c, kind=qp)) <= real(s%r, qp))) ok = .false.
        if (.not. (abs(v - cmplx(sf%c, kind=qp)) <= real(sf%r, qp))) ok = .false.
    end do
    call require(ok, "weighted inner products contain the real128 value", nfail)

    ! Wiener reciprocal: 1/F at sample points lies within the result.
    cc = [(0.0_dp, 0.0_dp), (0.0_dp, 0.15_dp), (0.5_dp, 0.0_dp), (2.0_dp, 0.0_dp), &
        (0.5_dp, 0.0_dp), (0.0_dp, -0.15_dp), (0.0_dp, 0.0_dp)]
    f = fs_trig(cc, 3)
    f%tau = 1.0e-10_dp
    y = fs_inverse(f, 40)
    ok = y%tau < 1.0e-8_dp
    do it = 0, 63
        th = 2.0_qp*acos(-1.0_qp)*real(it, qp)/64.0_qp
        z = (0.0_qp, 0.0_qp)
        do n = -3, 3
            z = z + cmplx(cc(n), kind=qp)*exp(cmplx(0.0_qp, n*th, qp))
        end do
        z = z + 1.0e-10_qp*cos(5*th)
        v = (0.0_qp, 0.0_qp)
        do n = -40, 40
            v = v + cmplx(y%a(n)%c, kind=qp)*exp(cmplx(0.0_qp, n*th, qp))
        end do
        if (abs(1.0_qp/z - v) > real(fs_l1_hi(y), qp) - sum(abs( &
            cmplx(y%a%c, kind=qp)))) ok = .false.
    end do
    call require(ok, "Wiener reciprocal encloses 1/F with a small tail", nfail)
    cc = [(0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
        (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), (0.0_dp, 0.0_dp)]
    y = fs_inverse(fs_trig(cc, 3), 10)
    call require(y%tau >= huge(1.0_dp), "reciprocal of a vanishing F is refused", &
        nfail)

    ! 2D products with sigma-weighted tails.
    ok = .true.
    ok2 = .true.
    do it = 1, 10
        f2 = random_series2(4, 3, 0.2_dp, 1.0e-9_dp, 1.0e-8_dp)
        g2 = random_series2(3, 2, 0.2_dp, 1.0e-10_dp, 1.0e-9_dp)
        call member2(f2, c2f)
        call member2(g2, c2g)
        call conv2(c2f, c2g, c2p)
        h2 = fs2_mul(f2, g2, 5, 3)
        hf2 = fs2_mul_fast(f2, g2, 5, 3)
        if (.not. contains2(h2, c2p)) ok = .false.
        if (.not. contains2(hf2, c2p)) ok2 = .false.
    end do
    call require(ok, "fs2_mul contains the exact weighted product", nfail)
    call require(ok2, "fs2_mul_fast contains the exact weighted product", nfail)
    call require(fs2_wnorm_hi(f2) >= fs2_sup_hi(f2), "weighted norm dominates", &
        nfail)
    s = fs2_inner_w(f2, g2, h2, fs2_sup_hi(h2))
    sf = fs2_inner_w_fast(f2, g2, h2, fs2_sup_hi(h2), 8, 6)
    call require(abs(s%c - sf%c) <= s%r + sf%r, "2D inner products agree", nfail)
    h2 = fs2_add(fs2_dtheta(f2), fs2_dphi(g2))
    call require(h2%tau >= huge(1.0_dp), "2D derivative of a tail is unbounded", &
        nfail)

    ! 2D reciprocal of B = 1 + 0.2 cos(theta) + 0.1 cos(theta - 5 phi).
    c2 = (0.0_dp, 0.0_dp)
    c2(0, 0) = (1.0_dp, 0.0_dp)
    c2(1, 0) = (0.1_dp, 0.0_dp)
    c2(-1, 0) = (0.1_dp, 0.0_dp)
    c2(1, 1) = (0.05_dp, 0.0_dp)
    c2(-1, -1) = (0.05_dp, 0.0_dp)
    f2 = fs2_trig_tail(c2, 2, 1, 1, 1, 5, 0.3_dp)
    y2 = fs2_inverse(f2, 12, 12)
    ok = y2%tau < 1.0e-6_dp
    do it = 0, 15
        th = 2.0_qp*acos(-1.0_qp)*real(it, qp)/16.0_qp
        ph = 0.37_qp*real(it, qp)
        z = 1.0_qp + 0.2_qp*cos(th) + 0.1_qp*cos(th - 5*ph)
        v = (0.0_qp, 0.0_qp)
        do j = -12, 12
            do n = -12, 12
                v = v + cmplx(y2%a(n, j)%c, kind=qp)*exp(cmplx(0.0_qp, &
                    n*th - j*5*ph, qp))
            end do
        end do
        if (abs(1.0_qp/z - v) > real(y2%tau, qp) + 1.0e-13_qp) ok = .false.
    end do
    call require(ok, "2D Wiener reciprocal encloses 1/B", nfail)

    deallocate (seed)
    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " fseries test(s) failed"
        error stop 1
    end if
    print '(a)', "fortnum_fseries: all tests passed"

contains

    function random_series(k, rmax, tau) result(fs)
        integer, intent(in) :: k
        real(dp), intent(in) :: rmax, tau
        type(fseries_t) :: fs
        real(dp) :: r(3)
        integer :: i

        fs = fs_zero(k)
        do i = -k, k
            call random_number(r)
            fs%a(i) = ball(cmplx(r(1) - 0.5_dp, r(2) - 0.5_dp, dp)*0.7_dp**abs(i), &
                rmax*r(3))
        end do
        fs%tau = tau
    end function random_series

    function worst_series(k, r, tau) result(fs)
        integer, intent(in) :: k
        real(dp), intent(in) :: r, tau
        type(fseries_t) :: fs
        integer :: i

        fs = fs_zero(k)
        do i = -k, k
            fs%a(i) = ball(cmplx(0.7_dp**abs(i), 0.0_dp, dp), r)
        end do
        fs%tau = tau
    end function worst_series

    subroutine member_edge(fs, k, c)
        type(fseries_t), intent(in) :: fs
        integer, intent(in) :: k
        complex(qp), intent(out) :: c(-kt:kt)
        integer :: i

        c = (0.0_qp, 0.0_qp)
        do i = -k, k
            c(i) = cmplx(fs%a(i)%c, kind=qp) + real(fs%a(i)%r, qp)*(1 - 1.0e-6_qp)
        end do
        c(k + 1) = real(fs%tau, qp)*(1 - 1.0e-6_qp)/2
        c(-k - 1) = c(k + 1)
    end subroutine member_edge

    !> A concrete member: points inside the balls and a tail of l1 norm
    !> just below tau on the indices k+1..k+3 of both signs.
    subroutine member(fs, k, c)
        type(fseries_t), intent(in) :: fs
        integer, intent(in) :: k
        complex(qp), intent(out) :: c(-kt:kt)
        integer :: i

        c = (0.0_qp, 0.0_qp)
        do i = -k, k
            c(i) = cmplx(fs%a(i)%c, kind=qp) + disc(fs%a(i)%r)
        end do
        do i = k + 1, k + 3
            c(i) = disc(fs%tau/6.0_dp, .true.)
            c(-i) = disc(fs%tau/6.0_dp, .true.)
        end do
    end subroutine member

    !> A random point of the disc of radius r (on its boundary if edge),
    !> shrunk by a relative 1e-6.
    function disc(r, edge) result(z)
        real(dp), intent(in) :: r
        logical, intent(in), optional :: edge
        complex(qp) :: z
        real(dp) :: s(2)
        real(qp) :: rho

        call random_number(s)
        rho = real(r, qp)*real(s(1), qp)
        if (present(edge)) rho = real(r, qp)
        rho = rho*(1.0_qp - 1.0e-6_qp)
        z = rho*exp(cmplx(0.0_qp, 2.0_qp*acos(-1.0_qp)*real(s(2), qp), qp))
    end function disc

    subroutine conv1(a, b, p)
        complex(qp), intent(in) :: a(-kt:kt), b(-kt:kt)
        complex(qp), intent(out) :: p(-2*kt:2*kt)
        integer :: i, jj

        p = (0.0_qp, 0.0_qp)
        do jj = -kt, kt
            do i = -kt, kt
                p(i + jj) = p(i + jj) + a(i)*b(jj)
            end do
        end do
    end subroutine conv1

    !> Membership of the exact coefficients c in the series hs.
    logical function contains1(hs, c, kc) result(inside)
        type(fseries_t), intent(in) :: hs
        integer, intent(in) :: kc
        complex(qp), intent(in) :: c(-kc:kc)
        real(qp) :: excess
        integer :: i

        excess = 0.0_qp
        do i = -kc, kc
            if (abs(i) <= hs%k) then
                excess = excess + max(0.0_qp, abs(c(i) - cmplx(hs%a(i)%c, kind=qp)) &
                    - real(hs%a(i)%r, qp))
            else
                excess = excess + abs(c(i))
            end if
        end do
        inside = excess <= real(hs%tau, qp) + 1.0e-28_qp
    end function contains1

    function random_series2(mk, nk, sigma, rmax, tau) result(fs)
        integer, intent(in) :: mk, nk
        real(dp), intent(in) :: sigma, rmax, tau
        type(fseries2d_t) :: fs
        real(dp) :: r(3)
        integer :: i, jj

        fs = fs2_zero(mk, nk, 5, sigma)
        do jj = -nk, nk
            do i = -mk, mk
                call random_number(r)
                fs%a(i, jj) = ball(cmplx(r(1) - 0.5_dp, r(2) - 0.5_dp, dp) &
                    *0.6_dp**(abs(i) + abs(jj)), rmax*r(3))
            end do
        end do
        fs%tau = tau
    end function random_series2

    !> A member: points in the balls plus a tail of weighted l1 norm just
    !> below tau at (mk + 1, nk + 1) and (-mk - 2, 0).
    subroutine member2(fs, c)
        type(fseries2d_t), intent(in) :: fs
        complex(qp), intent(out) :: c(-8:8, -6:6)
        integer :: i, jj

        c = (0.0_qp, 0.0_qp)
        do jj = -fs%nk, fs%nk
            do i = -fs%mk, fs%mk
                c(i, jj) = cmplx(fs%a(i, jj)%c, kind=qp) + disc(fs%a(i, jj)%r)
            end do
        end do
        c(fs%mk + 1, fs%nk + 1) = disc(fs%tau/2.0_dp, .true.) &
            /exp(real(fs%sigma, qp)*(fs%mk + fs%nk + 2))
        c(-fs%mk - 2, 0) = disc(fs%tau/2.0_dp, .true.) &
            /exp(real(fs%sigma, qp)*(fs%mk + 2))
    end subroutine member2

    subroutine conv2(a, b, p)
        complex(qp), intent(in) :: a(-8:8, -6:6), b(-8:8, -6:6)
        complex(qp), intent(out) :: p(-16:16, -12:12)
        integer :: i1, j1, i2, j2

        p = (0.0_qp, 0.0_qp)
        do j2 = -6, 6
            do i2 = -8, 8
                do j1 = -6, 6
                    do i1 = -8, 8
                        p(i1 + i2, j1 + j2) = p(i1 + i2, j1 + j2) + a(i1, j1)*b(i2, j2)
                    end do
                end do
            end do
        end do
    end subroutine conv2

    logical function contains2(hs, c) result(inside)
        type(fseries2d_t), intent(in) :: hs
        complex(qp), intent(in) :: c(-16:16, -12:12)
        real(qp) :: excess, wq
        integer :: i, jj

        excess = 0.0_qp
        do jj = -12, 12
            do i = -16, 16
                wq = exp(real(hs%sigma, qp)*(abs(i) + abs(jj)))
                if (abs(i) <= hs%mk .and. abs(jj) <= hs%nk) then
                    excess = excess + wq*max(0.0_qp, abs(c(i, jj) &
                        - cmplx(hs%a(i, jj)%c, kind=qp)) - real(hs%a(i, jj)%r, qp))
                else
                    excess = excess + wq*abs(c(i, jj))
                end if
            end do
        end do
        inside = excess <= real(hs%tau, qp) + 1.0e-28_qp
    end function contains2

    subroutine require(cond, msg, nfail)
        logical, intent(in) :: cond
        character(*), intent(in) :: msg
        integer, intent(inout) :: nfail

        if (.not. cond) then
            write (error_unit, '(a,a)') "FAIL: ", msg
            nfail = nfail + 1
        end if
    end subroutine require
end program test_fortnum_fseries
