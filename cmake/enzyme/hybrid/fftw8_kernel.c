#include <fftw3.h>

typedef struct {
    double primal;
    double tangent;
} enzyme_pair;

extern void fortnum_enzyme_rule_counter_record(void);

static fftw_plan plan;
static fftw_complex plan_input[8];
static fftw_complex plan_output[8];

void fortnum_fftw8_init(void)
{
    if (plan == 0) {
        plan = fftw_plan_dft_1d(
            8, plan_input, plan_output, FFTW_FORWARD,
            FFTW_ESTIMATE | FFTW_UNALIGNED);
    }
}

void fortnum_fftw8_finalize(void)
{
    if (plan != 0) {
        fftw_destroy_plan(plan);
        plan = 0;
    }
}

static void transform(const double *input, fftw_complex output[8])
{
    fftw_complex values[8];
    int i;

    fortnum_fftw8_init();
    for (i = 0; i < 8; ++i) {
        values[i][0] = input[2*i];
        values[i][1] = input[2*i + 1];
    }
    fftw_execute_dft(plan, values, output);
}

static double objective_from_transform(const fftw_complex values[8])
{
    double result = 0.0;
    int i;

    for (i = 0; i < 8; ++i) {
        result += 0.5*values[i][0]*values[i][0]
            + 0.25*values[i][1]*values[i][1];
    }
    return result;
}

double fortnum_fftw8_objective_array(const double *input)
{
    fftw_complex values[8];

    transform(input, values);
    return objective_from_transform(values);
}

double fortnum_fftw8_objective_jvp(
    const double *input, const double *direction)
{
    fftw_complex values[8];
    fftw_complex tangent[8];
    double product = 0.0;
    int i;

    transform(input, values);
    transform(direction, tangent);
    for (i = 0; i < 8; ++i) {
        product += values[i][0]*tangent[i][0]
            + 0.5*values[i][1]*tangent[i][1];
    }
    return product;
}

double fortnum_fftw8_objective_scalar(
    double x1, double x2, double x3, double x4,
    double x5, double x6, double x7, double x8,
    double x9, double x10, double x11, double x12,
    double x13, double x14, double x15, double x16)
{
    const double input[16] = {
        x1, x2, x3, x4, x5, x6, x7, x8,
        x9, x10, x11, x12, x13, x14, x15, x16
    };
    return fortnum_fftw8_objective_array(input);
}

enzyme_pair fortnum_fftw8_objective_scalar_derivative(
    double x1, double dx1, double x2, double dx2,
    double x3, double dx3, double x4, double dx4,
    double x5, double dx5, double x6, double dx6,
    double x7, double dx7, double x8, double dx8,
    double x9, double dx9, double x10, double dx10,
    double x11, double dx11, double x12, double dx12,
    double x13, double dx13, double x14, double dx14,
    double x15, double dx15, double x16, double dx16)
{
    const double input[16] = {
        x1, x2, x3, x4, x5, x6, x7, x8,
        x9, x10, x11, x12, x13, x14, x15, x16
    };
    const double direction[16] = {
        dx1, dx2, dx3, dx4, dx5, dx6, dx7, dx8,
        dx9, dx10, dx11, dx12, dx13, dx14, dx15, dx16
    };
    enzyme_pair result = {
        fortnum_fftw8_objective_array(input),
        fortnum_fftw8_objective_jvp(input, direction)
    };
    fortnum_enzyme_rule_counter_record();
    return result;
}
