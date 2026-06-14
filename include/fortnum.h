/* fortnum C ABI.
 *
 * C and C++ link target for the fortnum numerical library: a drop-in source
 * for the reference routines that KAMEL (KiLCA + QL-Balance) and the MEPHIT C
 * sources use. Every prototype here matches a bind(c) wrapper in
 * src/bindings/fortnum_capi.f90.
 *
 * Status convention: routines that can fail return an int status code
 * (FORTNUM_OK == 0; positive values are domain or convergence failures, see
 * the FORTNUM_* enum). Pure-domain scalar routines return the value directly.
 *
 * Complex arguments use the C99 _Complex double type, which is layout- and
 * ABI-compatible with Fortran complex(c_double_complex). C++ translation units
 * include <complex> and alias to double _Complex via the guard below.
 *
 * Arrays are passed by pointer; the caller owns the storage and the documented
 * length. Callback-driven routines take a function pointer plus an opaque
 * void* context forwarded unchanged to every callback invocation.
 */
#ifndef FORTNUM_H
#define FORTNUM_H

#ifdef __cplusplus
#include <complex>
typedef std::complex<double> fortnum_complex;
extern "C" {
#else
#include <complex.h>
typedef double _Complex fortnum_complex;
#endif

/* Status codes (mirror fortnum_status FORTNUM_*). */
enum {
    FORTNUM_OK = 0,
    FORTNUM_DOMAIN_ERROR = 1,
    FORTNUM_CONVERGENCE_ERROR = 2,
    FORTNUM_NOT_IMPLEMENTED = 3
};

/* Callback ABIs. */
typedef double (*fortnum_scalar_fn)(double x, void *ctx);
typedef void   (*fortnum_vector_fn)(int n, const double *x, double *f,
                                    void *ctx);
typedef void   (*fortnum_ode_rhs)(double t, int n, const double *y,
                                  double *dydt, void *ctx);

/* --- Special functions (real) --- */

/* Modified Bessel I_n(x), integer order, real argument. */
double fortnum_bessel_in(int n, double x);

/* Fill values[0..nmax] with I_0(x) .. I_nmax(x) in one pass. */
void fortnum_bessel_in_array(int nmax, double x, double *values);

/* Modified Bessel K_n(x), integer order, x > 0. */
double fortnum_bessel_kn(int n, double x);

/* Dawson integral F(x) = exp(-x^2) int_0^x exp(t^2) dt. */
double fortnum_dawson(double x);

/* Lower incomplete gamma gamma_lower(a,x) = int_0^x t^{a-1} e^{-t} dt. */
double fortnum_gamma_lower(double a, double x);

/* Regularized lower incomplete gamma P(a,x) = gamma_lower(a,x)/Gamma(a). */
double fortnum_gamma_reg_p(double a, double x);

/* Error function and complementary error function (C-callable erf/erfc). */
double fortnum_erf(double x);
double fortnum_erfc(double x);

/* --- Special functions (complex Bessel) --- */

/* J_n(z), integer order, complex z. Returns status; value in *result. */
int fortnum_bessel_j_complex(int order, fortnum_complex z,
                             fortnum_complex *result);

/* I_n(z); scaled != 0 returns e^{-Re z} I_n(z). */
int fortnum_bessel_i_complex(int order, fortnum_complex z, int scaled,
                             fortnum_complex *result);

/* Contiguous order sequence I_{order0..order0+nseq-1}(z) into result[nseq]. */
int fortnum_bessel_i_complex_array(int order0, int nseq, fortnum_complex z,
                                   int scaled, fortnum_complex *result);

/* K_n(z), Re z > 0; scaled != 0 returns e^{z} K_n(z). */
int fortnum_bessel_k_complex(int order, fortnum_complex z, int scaled,
                             fortnum_complex *result);

/* Contiguous order sequence K_{order0..order0+nseq-1}(z) into result[nseq]. */
int fortnum_bessel_k_complex_array(int order0, int nseq, fortnum_complex z,
                                   int scaled, fortnum_complex *result);

/* --- Confluent hypergeometric 1F1(a;b;z) (Kummer M) --- */

int fortnum_hyperg_1f1(fortnum_complex a, fortnum_complex b, fortnum_complex z,
                       fortnum_complex *result);

/* Specialization a = 1 used by MEPHIT/KiLCA. */
int fortnum_hyperg_1f1_a1(fortnum_complex b, fortnum_complex z,
                          fortnum_complex *result);

/* --- Fixed-rule Gauss-Legendre quadrature --- */

/* Nodes x[n] (ascending) and weights w[n] on [-1, 1]. */
void fortnum_gauss_legendre(int n, double *x, double *w);

/* Nodes x[n] and weights w[n] mapped to [a, b]. */
void fortnum_gauss_legendre_ab(int n, double a, double b, double *x, double *w);

/* Levin-u acceleration of the partial sums of terms[0..n-1]. */
int fortnum_levin_u_accel(const double *terms, int n, double *sum_accel,
                          double *abserr);

/* --- Adaptive quadrature (QUADPACK pattern) --- */

/* Globally adaptive Gauss-Kronrod, rule key in {15,21,31,61}, no extrapolation. */
int fortnum_integrate_qag(fortnum_scalar_fn f, double a, double b,
                          double epsabs, double epsrel, int key,
                          double *value, double *abserr);

/* Adaptive bisection plus Wynn epsilon extrapolation (GK21). */
int fortnum_integrate_qags(fortnum_scalar_fn f, double a, double b,
                           double epsabs, double epsrel,
                           double *value, double *abserr);

/* QAGS seeded with npts user break points so each known singularity starts
 * on a panel boundary. */
int fortnum_integrate_qagp(fortnum_scalar_fn f, double a, double b,
                           const double *points, int npts,
                           double epsabs, double epsrel,
                           double *value, double *abserr);

/* Semi-infinite / doubly infinite interval; inf in {-1, +1, +2}. */
int fortnum_integrate_qagiu(fortnum_scalar_fn f, double bound, int inf,
                            double epsabs, double epsrel,
                            double *value, double *abserr);

/* --- Root finding --- */

/* Brent's method on a bracketed interval [a, b]. */
int fortnum_root_brent(fortnum_scalar_fn f, double a, double b,
                       double xtol, double ftol, int max_iter, double *root);

/* Central finite-difference first derivative with Richardson error estimate. */
int fortnum_deriv_central(fortnum_scalar_fn f, double x, double h,
                          double *result, double *abserr);

/* Ascending index sort: perm[k] is the 0-based index into x giving the k-th
 * smallest value (x[perm] nondecreasing). perm has length n. */
void fortnum_argsort(const double *x, int n, int *perm);

/* Multidimensional Newton root finding F(x)=0 in R^n. fdf is a residual-only
 * callback; the wrapper builds the Jacobian by central differences. x0[n] is
 * the start, x[n] receives the converged root. */
int fortnum_multiroot_hybrid(fortnum_vector_fn fdf, int n, const double *x0,
                             double xtol, double ftol, int max_iter,
                             double *x);

/* --- ODE integration (Prince-Dormand RK8(7)13M / dop853) --- */

/* Integrate rhs from t0 to t1. The accepted-step mesh is copied into the
 * caller's buffers, capped at npts_cap; *npts reports the count actually
 * written. t_out has length npts_cap, y_out is column-major neq x npts_cap
 * (y_out[i + neq*k] is component i at mesh point k). */
int fortnum_ode_integrate_dop(fortnum_ode_rhs rhs, int neq, double t0,
                              double t1, const double *y0, double rtol,
                              double atol, int max_steps, int npts_cap,
                              double *t_out, double *y_out, int *npts);

/* Flat dop853 solve with default max_steps; same buffer-capped output. */
int fortnum_ode_solve_dop(fortnum_ode_rhs rhs, int neq, double t0, double t1,
                          const double *y0, double rtol, double atol,
                          int npts_cap, double *t_out, double *y_out,
                          int *npts);

/* --- B-spline basis (clamped knot vector, clamped-knot workflow) ---
 *
 * The workspace is opaque: create returns a handle, eval/deriv take it back,
 * destroy frees it. order = degree + 1; ncoef = nbreak + order - 2. */

/* Allocate a workspace; returns NULL on a domain error (order<2 or nbreak<2). */
void *fortnum_bspline_create(int order, int nbreak);

/* Build the clamped knot vector from nbreak strictly increasing breakpoints. */
int fortnum_bspline_set_knots(void *handle, int nbreak, const double *breakpts);

/* Number of basis functions; 0 for a NULL handle. */
int fortnum_bspline_ncoef(void *handle);

/* Evaluate all ncoef basis functions B_{i,k}(x) into values[ncoef]. */
int fortnum_bspline_eval_basis(void *handle, double x, int ncoef,
                               double *values);

/* Basis functions and derivatives up to order nderiv. dvalues is row-major
 * (nderiv+1) x ncoef: dvalues[d*ncoef + i] is the d-th derivative of basis
 * function i (d = 0 is the value). */
int fortnum_bspline_eval_deriv(void *handle, double x, int nderiv, int ncoef,
                               double *dvalues);

/* 1-based knot span index containing x. */
int fortnum_bspline_span_index(void *handle, double x);

/* Free the workspace behind the handle. */
void fortnum_bspline_destroy(void *handle);

#ifdef __cplusplus
}
#endif

#endif /* FORTNUM_H */
