program test_openmp_dawson_variants_offload
    use omp_lib, only: omp_is_initial_device
    use test_gpu_dawson_variants_support, only: run_dawson_variants_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$omp target map(from: on_device)
    on_device = .not. omp_is_initial_device()
    !$omp end target
    if (.not. on_device) error stop "OpenMP variants did not run on target"
    call run_dawson_variants_test()
end program test_openmp_dawson_variants_offload
