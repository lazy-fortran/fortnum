#include <stddef.h>

typedef struct {
    double primal;
    double tangent;
} enzyme_pair;

extern double fortnum_dawson_kernel(double x);
extern void fortnum_enzyme_rule_counter_record(void);

enzyme_pair fortnum_dawson_kernel_derivative(double x, double dx)
{
    const double f = fortnum_dawson_kernel(x);
    enzyme_pair result = {f, (1.0 - 2.0*x*f)*dx};
    fortnum_enzyme_rule_counter_record();
    return result;
}
