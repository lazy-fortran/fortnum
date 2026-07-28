typedef struct {
    double primal;
    double tangent;
} enzyme_pair;

extern double fortnum_hyperg_kernel_hybrid(double x);
extern double fortnum_hyperg_kernel_jvp(double x, double dx);
extern void fortnum_enzyme_rule_counter_record(void);

enzyme_pair fortnum_hyperg_kernel_hybrid_derivative(double x, double dx)
{
    enzyme_pair result = {
        fortnum_hyperg_kernel_hybrid(x),
        fortnum_hyperg_kernel_jvp(x, dx)
    };
    fortnum_enzyme_rule_counter_record();
    return result;
}
