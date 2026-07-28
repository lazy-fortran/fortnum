typedef struct {
    double primal;
    double tangent;
} enzyme_pair;

extern double fortnum_erf_kernel(double x);
extern double fortnum_erf_kernel_jvp(double x, double dx);
extern void fortnum_enzyme_rule_counter_record(void);

enzyme_pair fortnum_erf_kernel_derivative(double x, double dx)
{
    enzyme_pair result = {
        fortnum_erf_kernel(x),
        fortnum_erf_kernel_jvp(x, dx)
    };
    fortnum_enzyme_rule_counter_record();
    return result;
}
