module fortnum_special
    ! Public umbrella module re-exporting special functions.
    !
    ! This module aggregates the public APIs of the special function families:
    ! - Modified Bessel functions (fortnum_special_bessel)
    ! - Dawson integral (fortnum_special_dawson)
    ! - Lower incomplete gamma function (fortnum_special_gamma)
    ! - Confluent hypergeometric 1F1 (fortnum_special_hypergeometric_1f1)
    !
    ! Re-export layout preserves each module's derivative policy for future
    ! derivative product entry points (issue #40). Each routine documents its
    ! derivative policy in the originating module; this umbrella merely
    ! re-exports without modification.

    use fortnum_kinds, only: dp
    use fortnum_special_bessel, only: &
        bessel_in, &
        bessel_in_array, &
        bessel_kn, &
        bessel_in_jvp, &
        bessel_kn_jvp
    use fortnum_special_dawson, only: &
        dawson, &
        dawson_jvp, &
        dawson_grad
    use fortnum_special_gamma, only: &
        gamma_lower, &
        gamma_reg_p, &
        gamma_lower_jvp
    use fortnum_special_hypergeometric_1f1, only: &
        hyperg_1f1, &
        hyperg_1f1_a1, &
        hyperg_1f1_a1_jvp, &
        hyperg_1f1_a1_vjp

    implicit none
    private

    ! Bessel functions: modified I_n and K_n (integer order)
    ! Derivative policy: analytic_rule (ad.md §1, §4).
    ! Active argument: x only (n is inactive order selector).
    public :: bessel_in, bessel_in_array, bessel_kn
    public :: bessel_in_jvp, bessel_kn_jvp

    ! Dawson integral F(x) = exp(-x^2) * integral_0^x exp(t^2) dt
    ! Derivative policy: analytic_rule; F'(x) = 1 - 2*x*F(x).
    ! Active argument: x.
    public :: dawson
    public :: dawson_jvp, dawson_grad

    ! Incomplete gamma functions
    ! Derivative policy: analytic_rule.
    ! Active arguments: a (shape), x (integration limit).
    ! gamma_lower(a,x) = P(a,x)*Gamma(a), unnormalized lower incomplete.
    ! gamma_reg_p(a,x) = gamma_lower(a,x)/Gamma(a), regularized form.
    public :: gamma_lower, gamma_reg_p
    ! d/dx only; d/da deferred (requires digamma, see gamma_lower module header).
    public :: gamma_lower_jvp

    ! Confluent hypergeometric 1F1(a;b;z) (Kummer M), complex a, b, z.
    ! Derivative policy: analytic_rule; d/dz M = (a/b) M(a+1,b+1,z).
    ! Active argument: z.  Inactive: a, b.  hyperg_1f1_a1 fixes a = 1.
    public :: hyperg_1f1, hyperg_1f1_a1
    public :: hyperg_1f1_a1_jvp, hyperg_1f1_a1_vjp

end module fortnum_special
