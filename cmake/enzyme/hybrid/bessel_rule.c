#include <stddef.h>

typedef struct {
    double primal;
    double tangent;
} enzyme_pair;

extern double fortnum_bessel_i0_kernel(double x);
extern double fortnum_bessel_i1_kernel(double x);

static int rule_calls;
static int count_enabled = 1;

enzyme_pair fortnum_bessel_i0_kernel_derivative(double x, double dx)
{
    const double value = fortnum_bessel_i0_kernel(x);
    const double derivative = fortnum_bessel_i1_kernel(x);
    enzyme_pair result = {value, derivative*dx};
    if (count_enabled) {
        ++rule_calls;
    }
    return result;
}

void fortnum_bessel_rule_reset(void)
{
    rule_calls = 0;
    count_enabled = 1;
}

int fortnum_bessel_rule_calls(void)
{
    return rule_calls;
}

void fortnum_bessel_rule_disable_count(void)
{
    count_enabled = 0;
}
