#include <stddef.h>

typedef struct {
    double primal;
    double tangent;
} enzyme_pair;

extern double fortnum_bessel_i0_kernel(double x);
extern double fortnum_bessel_i1_kernel(double x);
extern void fortnum_enzyme_rule_counter_record(void);

enzyme_pair fortnum_bessel_i0_kernel_derivative(double x, double dx)
{
    const double value = fortnum_bessel_i0_kernel(x);
    const double derivative = fortnum_bessel_i1_kernel(x);
    enzyme_pair result = {value, derivative*dx};
    fortnum_enzyme_rule_counter_record();
    return result;
}
