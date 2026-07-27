#include <stddef.h>

typedef struct {
    double primal;
    double tangent;
} enzyme_pair;

extern double fortnum_dawson_kernel(double x);

static int rule_calls;
static int count_enabled = 1;

enzyme_pair fortnum_dawson_kernel_derivative(double x, double dx)
{
    const double f = fortnum_dawson_kernel(x);
    enzyme_pair result = {f, (1.0 - 2.0*x*f)*dx};
    if (count_enabled) {
        ++rule_calls;
    }
    return result;
}

void fortnum_dawson_rule_reset(void)
{
    rule_calls = 0;
    count_enabled = 1;
}

int fortnum_dawson_rule_calls(void)
{
    return rule_calls;
}

void fortnum_dawson_rule_disable_count(void)
{
    count_enabled = 0;
}
