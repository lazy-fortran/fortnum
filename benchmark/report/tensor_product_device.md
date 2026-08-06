# OpenACC resident tensor-product evidence

This benchmark exercises the OpenACC path added to
`fortnum_tensor_product`. `enter_data(status, n_rhs)` copies the three square
factors and allocates persistent vector and four-right-hand-side workspaces.
The resident timing keeps the input, output, factors, and contraction work on
the GPU across all repeated products. The transfer timing opens a data region
for each product, so it includes input/output transfers but reuses the
operator-owned device data. The host timing calls the ordinary CPU tensor
contraction with the same float64 factors and four right-hand sides.

The small (8^3=512) case is checked against an explicitly assembled dense
Kronecker matrix. The (16^3=4096) case is checked by the same contraction
setup and does not materialize its dense matrix. Both cases use one warm-up and
the repetitions shown below, pinned to CPU 4 for host control. The independent
behavioral test also runs both vector and multi-RHS device products against a
hand-built dense oracle.

| samples | host (ms) | device + transfer (ms) | device resident (ms) | host / resident | transfer / resident |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 512 | 0.0369 | 0.0811 | 0.0652 | 0.57x | 1.24x |
| 4096 | 0.8427 | 0.1779 | 0.1261 | 6.68x | 1.41x |

At 512 samples, GPU launch overhead dominates this small contraction. At 4096
samples, the resident GPU product is 6.68x faster than the nvfortran host
product, and the transfer-inclusive product is 4.74x faster. `NV_ACC_NOTIFY=1`
confirmed CUDA kernel launches for the copy, factor-contraction, and output
stages. Peak RSS is 118,436 KiB, the clean nvfortran build took 1.33 s, and
the linked executable is 335,762 bytes. Raw rows are in
[`tensor_product_device.csv`](data/tensor_product_device.csv).

This is a real accelerator result for the structured contraction, not an
OpenACC compilation-only claim. Toeplitz FFT products and the generic CG
recurrence still need device-resident implementations before the whole
structured GP solve can be called GPU-complete.
