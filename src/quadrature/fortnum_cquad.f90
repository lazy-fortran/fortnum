module fortnum_cquad
    ! Doubly-adaptive Clenshaw-Curtis quadrature.
    !
    ! Independent implementation of the doubly-adaptive scheme of Gonnet,
    ! "Increasing the Reliability of Adaptive Quadrature Using Explicit
    ! Interpolants" (ACM TOMS 37(3), 2010): on each subinterval the integrand
    ! is sampled at Chebyshev-Lobatto nodes, a Clenshaw-Curtis rule integrates
    ! the interpolant, and a nested lower-degree rule on the even-indexed subset
    ! gives the local error. The globally worst subinterval is bisected until the
    ! summed error meets the tolerance. Nodes and weights are built here from the
    ! Chebyshev definition (DLMF 3.5(v)); no external quadrature source.
    !
    ! The rule pair is degree 32 (33 Lobatto nodes) with the degree-16 subset
    ! (17 nodes) as the nested error estimator. Lobatto nodes are nested, so the
    ! coarse estimate reuses the fine samples with no extra integrand calls.
    !
    ! The local error is the L2 norm of the Chebyshev-coefficient difference
    ! between the fine and coarse interpolants (Gonnet section 4), not the
    ! difference of the two rule sums. The rule-sum difference degenerates to
    ! zero for integrands whose alias pattern makes both rules agree while the
    ! interpolant is still wrong (e.g. T40 - c*T44 on the top interval), which
    ! faked convergence with errors ten orders above the request. Segments
    ! whose coefficient error reaches the roundoff floor of the sampled
    ! magnitude are excluded from further bisection so cancellation-heavy
    ! integrands terminate instead of splitting forever.
    !
    ! No module-level mutable state: the node cosines and both weight vectors are
    ! built into call-local arrays, and the subinterval stack is caller-scoped,
    ! so concurrent calls are independent.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_integrate, only: integrate_integrand_t
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    implicit none
    private

    public :: integrate_cquad

    integer,  parameter :: NFINE = 32 ! fine rule degree (33 nodes)
    integer,  parameter :: NCOARSE = 16 ! nested coarse degree (17 nodes)
    integer,  parameter :: DEFAULT_LIMIT = 500 ! max subintervals
    real(dp), parameter :: PI = acos(-1.0_dp)
    real(dp), parameter :: EPMACH = epsilon(1.0_dp)
    ! A segment whose coefficient-difference norm is below this multiple of
    ! machine epsilon times its coefficient norm carries only roundoff; further
    ! bisection cannot lower it.
    real(dp), parameter :: FLOOR_FAC = 16.0_dp

contains

    ! I = integral of f over [a, b] to max(epsabs, epsrel*|I|).
    subroutine integrate_cquad(f, a, b, value, status, epsabs, epsrel, &
            abserr, limit, ctx)
        procedure(integrate_integrand_t)     :: f
        real(dp),               intent(in)    :: a, b
        real(dp),               intent(out)   :: value
        type(fortnum_status_t), intent(out)   :: status
        real(dp), intent(in),  optional :: epsabs, epsrel
        real(dp), intent(out), optional :: abserr
        integer,  intent(in),  optional :: limit
        class(*), intent(in),  optional :: ctx

        real(dp) :: epsabs_loc, epsrel_loc
        integer  :: limit_loc
        real(dp) :: cosf(0:NFINE, 0:NFINE) ! cos(j k pi / NFINE)
        real(dp) :: wfine(0:NFINE) ! degree-32 Clenshaw-Curtis
        real(dp), allocatable :: sa(:), sb(:), si(:), se(:), sn(:)
        real(dp) :: tot, toterr, tol, worst
        integer  :: nseg, iworst, i

        value = 0.0_dp
        if (present(abserr)) abserr = 0.0_dp
        call status_set(status, FORTNUM_OK, "")

        epsabs_loc = 0.0_dp
        epsrel_loc = 1.0e-8_dp
        if (present(epsabs)) epsabs_loc = epsabs
        if (present(epsrel)) epsrel_loc = epsrel
        limit_loc = DEFAULT_LIMIT
        if (present(limit)) limit_loc = limit

        if (.not. (a == a .and. b == b)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "cquad: non-finite interval")
            return
        end if
        if (a == b) return

        call build_rule(cosf, wfine)

        allocate (sa(limit_loc), sb(limit_loc), si(limit_loc), se(limit_loc), &
            sn(limit_loc))
        nseg = 1
        sa(1) = a
        sb(1) = b
        call panel(f, cosf, wfine, sa(1), sb(1), si(1), se(1), sn(1), ctx)
        tot = si(1)
        toterr = se(1)

        do
            tol = max(epsabs_loc, epsrel_loc*abs(tot))
            if (toterr <= tol) exit
            if (nseg >= limit_loc) then
                call status_set(status, FORTNUM_CONVERGENCE_ERROR, &
                    "cquad: subinterval limit reached before tolerance")
                exit
            end if

            ! Pick the worst segment still above its roundoff floor; segments
            ! at the floor cannot improve by bisection.
            iworst = 0
            worst = 0.0_dp
            do i = 1, nseg
                if (se(i) <= FLOOR_FAC*EPMACH*sn(i)) cycle
                if (se(i) > worst) then
                    worst = se(i)
                    iworst = i
                end if
            end do
            if (iworst == 0) exit

            nseg = nseg + 1
            sa(nseg) = 0.5_dp*(sa(iworst) + sb(iworst))
            sb(nseg) = sb(iworst)
            sb(iworst) = sa(nseg)

            tot = tot - si(iworst)
            toterr = toterr - se(iworst)
            call panel(f, cosf, wfine, sa(iworst), sb(iworst), &
                si(iworst), se(iworst), sn(iworst), ctx)
            call panel(f, cosf, wfine, sa(nseg), sb(nseg), &
                si(nseg), se(nseg), sn(nseg), ctx)
            tot = tot + si(iworst) + si(nseg)
            toterr = toterr + se(iworst) + se(nseg)
        end do

        value = tot
        if (present(abserr)) abserr = toterr
        deallocate (sa, sb, si, se, sn)
    end subroutine integrate_cquad

    ! One subinterval over [alpha,beta]: fine integral, coefficient-difference
    ! error estimate, and the coefficient-norm scale for the roundoff floor.
    subroutine panel(f, cosf, wfine, alpha, beta, ival, eval, fnorm, ctx)
        procedure(integrate_integrand_t)  :: f
        real(dp), intent(in)  :: cosf(0:NFINE, 0:NFINE)
        real(dp), intent(in)  :: wfine(0:NFINE)
        real(dp), intent(in)  :: alpha, beta
        real(dp), intent(out) :: ival, eval, fnorm
        class(*), intent(in), optional :: ctx

        real(dp) :: fx(0:NFINE), half, mid, ifine
        real(dp) :: cfine(0:NFINE), ccoarse(0:NCOARSE)
        real(dp) :: ncdiff, nc, t
        integer  :: j, k

        half = 0.5_dp*(beta - alpha)
        mid = 0.5_dp*(beta + alpha)
        do j = 0, NFINE
            fx(j) = f(mid + half*cosf(j, 1), ctx) ! cosf(j,1) = cos(j pi/NFINE)
            ! Non-finite samples (integrable endpoint singularities, 0*inf in
            ! integrand prefactors) are zeroed, as in GSL's CQUAD. The zeroed
            ! node leaves a kink the coefficient estimator sees as a large
            ! error, so bisection isolates the point and its contribution
            ! vanishes with the interval width.
            if (.not. (abs(fx(j)) <= huge(0.0_dp))) fx(j) = 0.0_dp
        end do

        ifine = 0.0_dp
        do j = 0, NFINE
            ifine = ifine + wfine(j)*fx(j)
        end do

        ! Chebyshev coefficients (scaled so p = sum_k c_k T_k) of the 33-node
        ! interpolant and of the nested 17-node interpolant on the even subset.
        ! cosf(NFINE,k) = cos(k pi) is the endpoint factor for both grids.
        do k = 0, NFINE
            t = 0.5_dp*(fx(0) + cosf(NFINE, k)*fx(NFINE))
            do j = 1, NFINE - 1
                t = t + fx(j)*cosf(j, k)
            end do
            cfine(k) = (2.0_dp/real(NFINE, dp))*t
            if (k == 0 .or. k == NFINE) cfine(k) = 0.5_dp*cfine(k)
        end do
        do k = 0, NCOARSE
            t = 0.5_dp*(fx(0) + cosf(NFINE, k)*fx(NFINE))
            do j = 1, NCOARSE - 1
                t = t + fx(2*j)*cosf(2*j, k) ! cos(j k pi/NCOARSE)
            end do
            ccoarse(k) = (2.0_dp/real(NCOARSE, dp))*t
            if (k == 0 .or. k == NCOARSE) ccoarse(k) = 0.5_dp*ccoarse(k)
        end do

        ! Gonnet-style estimate: L2 norm of the interpolant difference over
        ! the interval, with the divergence guard on strongly disagreeing
        ! interpolants (mirrors the published algorithm's err update).
        ncdiff = 0.0_dp
        do k = 0, NCOARSE
            t = ccoarse(k) - cfine(k)
            ncdiff = ncdiff + t*t
        end do
        nc = 0.0_dp
        do k = 0, NFINE
            nc = nc + cfine(k)*cfine(k)
            if (k > NCOARSE) ncdiff = ncdiff + cfine(k)*cfine(k)
        end do
        ncdiff = sqrt(ncdiff)
        nc = sqrt(nc)

        ival = half*ifine
        eval = 2.0_dp*half*ncdiff
        if (ncdiff > 0.1_dp*nc .and. eval < 2.0_dp*half*nc) then
            eval = 2.0_dp*half*nc
        end if
        fnorm = 2.0_dp*half*nc
    end subroutine panel

    ! Clenshaw-Curtis nodes/weights on [-1,1] for the degree-32 rule.
    ! Weight of node j: w_j = (2/N) s_j sum_{k even} mu_k cos(j k pi/N),
    ! s_j the endpoint-halving factor, mu_0 = 1, mu_k = 2/(1-k^2).
    ! cosf(j,1) doubles as the node abscissa cos(j pi/N).
    subroutine build_rule(cosf, wfine)
        real(dp), intent(out) :: cosf(0:NFINE, 0:NFINE)
        real(dp), intent(out) :: wfine(0:NFINE)
        integer :: j, k

        do k = 0, NFINE
            do j = 0, NFINE
                cosf(j, k) = cos(PI*real(j*k, dp)/real(NFINE, dp))
            end do
        end do

        call cc_weights(NFINE, wfine)
    end subroutine build_rule

    subroutine cc_weights(n, w)
        integer,  intent(in)  :: n
        real(dp), intent(out) :: w(0:n)
        real(dp) :: s, mu, ang
        integer  :: j, k

        do j = 0, n
            s = 0.0_dp
            do k = 0, n, 2
                if (k == 0) then
                    mu = 1.0_dp
                else
                    mu = 2.0_dp/(1.0_dp - real(k, dp)**2)
                end if
                ang = PI*real(j*k, dp)/real(n, dp)
                s = s + mu*cos(ang)
            end do
            w(j) = (2.0_dp/real(n, dp))*s
            if (j == 0 .or. j == n) w(j) = 0.5_dp*w(j)
        end do
    end subroutine cc_weights

end module fortnum_cquad
