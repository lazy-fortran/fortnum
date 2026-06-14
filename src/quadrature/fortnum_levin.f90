module fortnum_levin
    ! Levin u sequence-acceleration / nonlinear series transform.
    !
    ! Replaces the reference the Levin-u accelerator family used by MEPHIT (hyper1F1.c
    ! hypergeometric1f1_kummer_modified_0_accel) and KAMEL KiLCA
    ! (math/hyper/hyper1F1.cpp) to sum the Kummer-series ratio of a confluent
    ! hypergeometric function.  The caller passes a finite, recorded sequence
    ! of series terms a_0 .. a_{n-1}; the routine forms the partial sums and
    ! returns the Levin-u accelerated limit together with an error estimate.
    !
    ! Algorithm (clean-room from the published papers, no the external backend source):
    !   D. Levin, Int. J. Comput. Math. B3 (1973) 371-388.
    !   E. J. Weniger, Comput. Phys. Rep. 10 (1989) 189-371, eqs 7.2-8, 7.3-9.
    !   Fessler, Ford, Smith, ACM TOMS 9 (1983) 346-354 (the reference's algorithm).
    !
    !   With remainder estimate omega_n = (beta + n) * a_n (beta = 1, the
    !   "u" variant, Weniger 7.3-9) the L-transformation
    !       L_k^{(n)} = N_k^{(n)} / D_k^{(n)}
    !   is built from the column recurrences (Weniger 7.2-8)
    !       N_0^{(n)} = s_n / omega_n,   D_0^{(n)} = 1 / omega_n
    !       X_k^{(n)} = X_{k-1}^{(n+1)} - r * X_{k-1}^{(n)},  X in {N, D}
    !       r = ((beta + n) / (beta + n + k))^{k-1}.
    !   The fully accelerated value uses terms 0..n-1 and is N_{n-1}^{(0)} /
    !   D_{n-1}^{(0)}.  As in the accelerator routine the table is grown one term
    !   at a time and the order whose diagonal value moves the least, scaled by
    !   a rounding floor, is returned as the best estimate.
    !
    ! DERIVATIVE POLICY (ad.md): primal_only.  The accelerated value is a
    !   nonlinear rational transform of the term sequence, and the order that
    !   minimises the error estimate is data-dependent (a recorded trace).  No
    !   transparent linear or closed-form analytic derivative rule applies, so
    !   no jvp/vjp/grad routine is provided.  Active: terms.  Inactive: n.
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    implicit none
    private

    public :: levin_u_accel

    ! beta = 1 selects the u variant (remainder estimate (beta+n)*a_n).
    real(dp), parameter :: beta = 1.0_dp
    ! Rounding floor for the per-order error estimate, in units of eps.
    real(dp), parameter :: round_floor = 16.0_dp

contains

    ! Levin-u acceleration of the partial sums of the series whose terms are
    ! terms(1:n).  Returns the best accelerated value in sum_accel and an
    ! absolute error estimate in abserr.  status carries FORTNUM_DOMAIN_ERROR
    ! when n < 1 or a term that would divide by zero (a_i == 0) is supplied;
    ! otherwise FORTNUM_OK.  Call separately for the real and imaginary parts
    ! of a complex series, matching the consumer pattern.
    !
    ! Active: terms.  Inactive: n (the term count, fixed by the caller).
    subroutine levin_u_accel(terms, n, sum_accel, abserr, status)
        integer,                intent(in)  :: n
        real(dp),               intent(in)  :: terms(n)
        real(dp),               intent(out) :: sum_accel
        real(dp),               intent(out) :: abserr
        type(fortnum_status_t), intent(out) :: status

        real(dp) :: n0(n), d0(n)       ! k=0 column: s_m/omega_m, 1/omega_m
        real(dp) :: qnum(n), qden(n)   ! working N_k, D_k columns
        real(dp) :: psum, omega, ratio, fact, val, prev, est, best, best_err
        integer  :: i, k, j, terms_used

        sum_accel = 0.0_dp
        abserr    = 0.0_dp

        if (n < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "levin_u_accel: n must be >= 1")
            return
        end if
        do i = 1, n
            if (terms(i) == 0.0_dp) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "levin_u_accel: zero term breaks the u remainder estimate")
                return
            end if
        end do

        psum     = 0.0_dp
        prev     = 0.0_dp
        best     = 0.0_dp
        best_err = huge(1.0_dp)
        terms_used = 0

        do i = 1, n
            ! Index i corresponds to series index (i-1); omega_{i-1}.
            psum  = psum + terms(i)
            omega = (beta + real(i - 1, dp))*terms(i)
            n0(i) = psum/omega
            d0(i) = 1.0_dp/omega

            ! Restore the k=0 column, then build columns k = 1 .. i-1 over the
            ! live window 1..i.  Each column reads the previous one (index j+1
            ! before j is overwritten), so a single ascending j pass is exact.
            qnum(1:i) = n0(1:i)
            qden(1:i) = d0(1:i)
            do k = 1, i - 1
                do j = 1, i - k
                    ratio = (beta + real(j - 1, dp))/(beta + real(j - 1 + k, dp))
                    fact  = ratio**(k - 1)
                    qnum(j) = qnum(j + 1) - fact*qnum(j)
                    qden(j) = qden(j + 1) - fact*qden(j)
                end do
            end do

            val = qnum(1)/qden(1)

            ! Error estimate: change against the previous order's diagonal
            ! value, floored by the relative rounding level (the accelerator scheme).
            if (i == 1) then
                est = abs(val)
            else
                est = abs(val - prev)
            end if
            est = max(est, round_floor*epsilon(1.0_dp)*abs(val))
            if (est < best_err) then
                best_err   = est
                best       = val
                terms_used = i
            end if
            prev = val
        end do

        sum_accel = best
        abserr    = best_err
        call status_set(status, FORTNUM_OK, "")
    end subroutine levin_u_accel

end module fortnum_levin
