program test_openmp_fft8_offload
    use omp_lib, only: omp_is_initial_device
    use test_gpu_fft8_support, only: run_fft8_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$omp target map(from: on_device)
    on_device = .not. omp_is_initial_device()
    !$omp end target
    if (.not. on_device) error stop "OpenMP FFT8 did not run on target"
    call run_fft8_test()
end program test_openmp_fft8_offload
