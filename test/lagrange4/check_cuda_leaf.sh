#!/usr/bin/env bash
# Demonstrate that the committed CUDA Lagrange-4 JVP leaf agrees numerically
# with the independent cubic oracle, and that it is a launchable device leaf.
#
# The leaf (src/generated/cuda/fortnum_lagrange4_jvp_kernel_cuda.cu) is a
# bare extern "C" __device__ __forceinline__ function from
# fortsym_kernel_emit::emit_cuda_device_ir. It knows nothing about launch
# geometry; a backend-owned __global__ wrapper must call it.
#
# Two checks:
#   1. Host compile (always): macro the CUDA qualifiers away and compile the
#      exact leaf body with g++, then compare against the cubic oracle at the
#      kernel's own 1e-13 relative tolerance. This exercises the emitted
#      arithmetic on this runner even though nvcc is absent.
#   2. Device compile (only when nvcc is present): compile the real .cu with
#      nvcc, launch the leaf through a tiny __global__ wrapper, and compare
#      against the same oracle.
set -euo pipefail

repository_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
cuda_leaf="${repository_dir}/src/generated/cuda/fortnum_lagrange4_jvp_kernel_cuda.cu"
tmp_dir=$(mktemp -d)
trap 'rm -rf -- "$tmp_dir"' EXIT

if [[ ! -f "$cuda_leaf" ]]; then
    echo "committed CUDA leaf missing: $cuda_leaf" >&2
    exit 1
fi

cat > "$tmp_dir/cuda_host_driver.cpp" <<'CPP'
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>

// Strip the CUDA qualifiers so the exact emitted leaf body compiles as a
// plain host function with g++. The arithmetic is unchanged.
#define __device__
#define __forceinline__
#include "fortnum_lagrange4_jvp_kernel_cuda.cu"
#undef __device__
#undef __forceinline__

static double primal_polynomial(double z, int index) {
    double shift = double(index % 17) / 31.0;
    return (0.3 + shift) - 0.4*z + 0.2*z*z - 0.05*z*z*z;
}
static double tangent_polynomial(double z, int index) {
    double shift = double(index % 13) / 29.0;
    return (-0.2 + shift) + 0.3*z - 0.1*z*z + 0.04*z*z*z;
}
// Cubic interpolant coefficients through nodes {-1,0,1,2}.
static void cubic_coefficients(const double* s, double* c) {
    c[0] = s[1];
    c[1] = -s[0]/3.0 - s[1]/2.0 + s[2] - s[3]/6.0;
    c[2] =  s[0]/2.0 - s[1] + s[2]/2.0;
    c[3] = -s[0]/6.0 + s[1]/2.0 - s[2]/2.0 + s[3]/6.0;
}

int main() {
    const int n = 4096;
    const double tol = 1.0e-13;
    const double nodes[4] = {-1.0, 0.0, 1.0, 2.0};
    int failures = 0;
    for (int i = 1; i <= n; ++i) {
        double x  = -0.8 + 2.6*double((17*i) % 4093)/4092.0;
        double tx = -0.7 + 1.4*double((19*i) % 4091)/4090.0;
        double samples[4], ty[4];
        for (int j = 0; j < 4; ++j) {
            samples[j] = primal_polynomial(nodes[j], i);
            ty[j]      = tangent_polynomial(nodes[j], i);
        }
        double value, jvp;
        fortnum_lagrange4_jvp_kernel_cuda(
            x, samples[0], samples[1], samples[2], samples[3],
            tx, ty[0], ty[1], ty[2], ty[3], &value, &jvp);

        double c[4];
        cubic_coefficients(samples, c);
        double exp_value = c[0] + x*(c[1] + x*(c[2] + x*c[3]));
        double exp_jvp = tx*(c[1] + x*(2.0*c[2] + 3.0*x*c[3])) +
                         tangent_polynomial(x, i);
        double vscale = std::max(1.0, std::fabs(exp_value));
        double jscale = std::max(1.0, std::fabs(exp_jvp));
        if (std::fabs(value - exp_value) > tol*vscale) ++failures;
        if (std::fabs(jvp - exp_jvp) > tol*jscale) ++failures;
    }
    if (failures != 0) {
        std::fprintf(stderr, "CUDA leaf host check: %d mismatches\n", failures);
        return 1;
    }
    std::printf("CUDA leaf host check: leaf agrees with cubic oracle\n");
    return 0;
}
CPP

g++ -O2 -I "$(dirname "$cuda_leaf")" -o "$tmp_dir/cuda_host_check" \
    "$tmp_dir/cuda_host_driver.cpp"
"$tmp_dir/cuda_host_check"

if command -v nvcc >/dev/null 2>&1; then
    cat > "$tmp_dir/cuda_device_driver.cu" <<'CUDA'
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include "fortnum_lagrange4_jvp_kernel_cuda.cu"

static double primal_polynomial(double z, int index) {
    double shift = double(index % 17) / 31.0;
    return (0.3 + shift) - 0.4*z + 0.2*z*z - 0.05*z*z*z;
}
static double tangent_polynomial(double z, int index) {
    double shift = double(index % 13) / 29.0;
    return (-0.2 + shift) + 0.3*z - 0.1*z*z + 0.04*z*z*z;
}
static void cubic_coefficients(const double* s, double* c) {
    c[0] = s[1];
    c[1] = -s[0]/3.0 - s[1]/2.0 + s[2] - s[3]/6.0;
    c[2] =  s[0]/2.0 - s[1] + s[2]/2.0;
    c[3] = -s[0]/6.0 + s[1]/2.0 - s[2]/2.0 + s[3]/6.0;
}

__global__ void invoke(const double* args, double* out) {
    // args: x,y1..y4,tx,ty1..ty4; out: value,jvp
    fortnum_lagrange4_jvp_kernel_cuda(
        args[0], args[1], args[2], args[3], args[4],
        args[5], args[6], args[7], args[8], args[9],
        &out[0], &out[1]);
}

int main() {
    const int n = 4096;
    const double tol = 1.0e-13;
    const double nodes[4] = {-1.0, 0.0, 1.0, 2.0};
    double* d_args; double* d_out;
    if (cudaMalloc(&d_args, 10*sizeof(double)) != cudaSuccess) return 10;
    if (cudaMalloc(&d_out, 2*sizeof(double)) != cudaSuccess) return 11;
    int failures = 0;
    for (int i = 1; i <= n; ++i) {
        double x  = -0.8 + 2.6*double((17*i) % 4093)/4092.0;
        double tx = -0.7 + 1.4*double((19*i) % 4091)/4090.0;
        double samples[4], ty[4];
        for (int j = 0; j < 4; ++j) {
            samples[j] = primal_polynomial(nodes[j], i);
            ty[j]      = tangent_polynomial(nodes[j], i);
        }
        double h_args[10] = {x, samples[0], samples[1], samples[2], samples[3],
                             tx, ty[0], ty[1], ty[2], ty[3]};
        double h_out[2] = {0.0, 0.0};
        cudaMemcpy(d_args, h_args, 10*sizeof(double), cudaMemcpyHostToDevice);
        invoke<<<1, 1>>>(d_args, d_out);
        if (cudaDeviceSynchronize() != cudaSuccess) return 12;
        cudaMemcpy(h_out, d_out, 2*sizeof(double), cudaMemcpyDeviceToHost);

        double c[4];
        cubic_coefficients(samples, c);
        double exp_value = c[0] + x*(c[1] + x*(c[2] + x*c[3]));
        double exp_jvp = tx*(c[1] + x*(2.0*c[2] + 3.0*x*c[3])) +
                         tangent_polynomial(x, i);
        double vscale = std::max(1.0, std::fabs(exp_value));
        double jscale = std::max(1.0, std::fabs(exp_jvp));
        if (std::fabs(h_out[0] - exp_value) > tol*vscale) ++failures;
        if (std::fabs(h_out[1] - exp_jvp) > tol*jscale) ++failures;
    }
    cudaFree(d_args); cudaFree(d_out);
    if (failures != 0) {
        std::fprintf(stderr, "CUDA leaf device check: %d mismatches\n", failures);
        return 1;
    }
    std::printf("CUDA leaf device check: nvcc leaf agrees with cubic oracle\n");
    return 0;
}
CUDA
    nvcc -O2 -std=c++17 -I "$(dirname "$cuda_leaf")" \
        -o "$tmp_dir/cuda_device_check" "$tmp_dir/cuda_device_driver.cu"
    "$tmp_dir/cuda_device_check"
else
    echo "nvcc unavailable; CUDA device compile/run oracle skipped"
fi
