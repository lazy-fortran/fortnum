/* C smoke test for the fortnum C ABI.
 *
 * Includes the installed header, links the fortnum library, and checks a
 * representative slice of the surface against known values: a real Bessel, the
 * lower incomplete gamma, an adaptive quadrature with a C callback, a complex
 * Bessel, and the confluent hypergeometric 1F1. Exits nonzero on the first
 * failed assertion so CTest reports a hard fail.
 */
#include "fortnum.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

static int failures = 0;

static void check(const char *name, double got, double want, double tol) {
    double err = fabs(got - want);
    double rel = err / (fabs(want) + 1e-300);
    int ok = (err <= tol) || (rel <= tol);
    printf("%-28s got=%.15g want=%.15g abs_err=%.3g %s\n",
           name, got, want, err, ok ? "PASS" : "FAIL");
    if (!ok) failures++;
}

/* Integrand for the quadrature check: f(x) = exp(-x^2), int_0^1 = 0.5*sqrt(pi)*erf(1). */
static double gaussian(double x, void *ctx) {
    (void)ctx;
    return exp(-x * x);
}

/* Integrand that reads its parameter through ctx: f(x) = scale * x.
 * Used to verify the C ABI forwards the opaque void* unchanged to the
 * callback. int_0^1 scale*x dx = scale/2. */
static double scaled_linear(double x, void *ctx) {
    double scale = *(const double *)ctx;
    return scale * x;
}

int main(void) {
    /* I_0(1.0) = 1.2660658777520084 (DLMF 10.25). */
    check("bessel_in(0,1)", fortnum_bessel_in(0, 1.0),
          1.2660658777520084, 1e-12);

    /* K_1(2.0) = 0.13986588181652243. */
    check("bessel_kn(1,2)", fortnum_bessel_kn(1, 2.0),
          0.13986588181652243, 1e-11);

    /* gamma_lower(2, 1) = 1 - 2/e = 0.26424111765711533. */
    check("gamma_lower(2,1)", fortnum_gamma_lower(2.0, 1.0),
          1.0 - 2.0 / exp(1.0), 1e-12);

    /* P(2,1) = gamma_lower(2,1)/Gamma(2) = same since Gamma(2)=1. */
    check("gamma_reg_p(2,1)", fortnum_gamma_reg_p(2.0, 1.0),
          1.0 - 2.0 / exp(1.0), 1e-12);

    /* dawson(1.0) = 0.5380795069127684. */
    check("dawson(1)", fortnum_dawson(1.0), 0.5380795069127684, 1e-12);

    /* QAGS of exp(-x^2) on [0,1] = 0.5*sqrt(pi)*erf(1) = 0.7468241328124271. */
    {
        double val = 0.0, abserr = 0.0;
        int code = fortnum_integrate_qags(gaussian, 0.0, 1.0, 0.0, 1e-10,
                                          &val, &abserr, NULL);
        printf("integrate_qags status=%d\n", code);
        if (code != FORTNUM_OK) failures++;
        check("integrate_qags(gauss)", val, 0.7468241328124271, 1e-9);
    }

    /* QAG (rule 21) of the same integrand. */
    {
        double val = 0.0, abserr = 0.0;
        int code = fortnum_integrate_qag(gaussian, 0.0, 1.0, 0.0, 1e-10, 21,
                                         &val, &abserr, NULL);
        printf("integrate_qag status=%d\n", code);
        if (code != FORTNUM_OK) failures++;
        check("integrate_qag(gauss)", val, 0.7468241328124271, 1e-9);
    }

    /* Context forwarding: integrate scale*x with scale passed through ctx.
     * int_0^1 scale*x dx = scale/2; here scale=3 so the result is 1.5. */
    {
        double scale = 3.0, val = 0.0, abserr = 0.0;
        int code = fortnum_integrate_qag(scaled_linear, 0.0, 1.0, 0.0, 1e-10,
                                         21, &val, &abserr, &scale);
        printf("integrate_qag(ctx) status=%d\n", code);
        if (code != FORTNUM_OK) failures++;
        check("integrate_qag(ctx scale*x)", val, 1.5, 1e-9);
    }

    /* Complex I_0(1+0i) must match the real I_0(1). */
    {
        fortnum_complex z, r;
#ifdef __cplusplus
        z = fortnum_complex(1.0, 0.0);
#else
        z = 1.0 + 0.0 * I;
#endif
        int code = fortnum_bessel_i_complex(0, z, 0, &r);
        printf("bessel_i_complex status=%d\n", code);
        if (code != FORTNUM_OK) failures++;
#ifdef __cplusplus
        double re = r.real();
#else
        double re = creal(r);
#endif
        check("bessel_i_complex(0,1)", re, 1.2660658777520084, 1e-10);
    }

    /* 1F1(1; 1; 1) = e (M(a;a;z) = e^z). */
    {
        fortnum_complex a, b, zz, r;
#ifdef __cplusplus
        a = fortnum_complex(1.0, 0.0);
        b = fortnum_complex(1.0, 0.0);
        zz = fortnum_complex(1.0, 0.0);
#else
        a = 1.0 + 0.0 * I;
        b = 1.0 + 0.0 * I;
        zz = 1.0 + 0.0 * I;
#endif
        int code = fortnum_hyperg_1f1(a, b, zz, &r);
        printf("hyperg_1f1 status=%d\n", code);
        if (code != FORTNUM_OK) failures++;
#ifdef __cplusplus
        double re = r.real();
#else
        double re = creal(r);
#endif
        check("hyperg_1f1(1;1;1)", re, exp(1.0), 1e-12);
    }

    /* erf(1) = 0.8427007929497149; erfc(1) = 1 - erf(1). */
    check("erf(1)", fortnum_erf(1.0), 0.8427007929497149, 1e-12);
    check("erfc(1)", fortnum_erfc(1.0), 1.0 - 0.8427007929497149, 1e-12);

    /* B-spline partition of unity: order-4 (cubic) on [0,1,2,3,4], the basis
     * functions at any interior x must sum to 1. */
    {
        double breaks[5] = {0.0, 1.0, 2.0, 3.0, 4.0};
        void *bs = fortnum_bspline_create(4, 5);
        if (!bs) {
            printf("bspline_create returned NULL\n");
            failures++;
        } else {
            int code = fortnum_bspline_set_knots(bs, 5, breaks);
            if (code != FORTNUM_OK) failures++;
            int nc = fortnum_bspline_ncoef(bs);
            double *vals = (double *)malloc((size_t)nc * sizeof(double));
            code = fortnum_bspline_eval_basis(bs, 1.7, nc, vals);
            if (code != FORTNUM_OK) failures++;
            double sum = 0.0;
            for (int i = 0; i < nc; ++i) sum += vals[i];
            printf("bspline ncoef=%d basis status=%d\n", nc, code);
            check("bspline partition_of_unity", sum, 1.0, 1e-12);
            free(vals);
            fortnum_bspline_destroy(bs);
        }
    }

    if (failures == 0) {
        printf("ALL CAPI SMOKE CHECKS PASSED\n");
        return 0;
    }
    printf("%d CAPI SMOKE CHECK(S) FAILED\n", failures);
    return 1;
}
