#include <stddef.h>

typedef struct {
    double primal;
    double tangent;
} enzyme_pair;

extern double fortnum_gamma_reg_p_kernel(double x);
extern double fortnum_gamma_reg_p_kernel_jvp(double x, double dx);
extern void fortnum_enzyme_rule_counter_record(void);

enzyme_pair fortnum_gamma_reg_p_kernel_derivative(double x, double dx)
{
    enzyme_pair result = {
        fortnum_gamma_reg_p_kernel(x),
        fortnum_gamma_reg_p_kernel_jvp(x, dx)
    };
    fortnum_enzyme_rule_counter_record();
    return result;
}
