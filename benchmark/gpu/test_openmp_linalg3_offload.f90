program test_openmp_linalg3_offload
    use omp_lib, only: omp_is_initial_device
    use test_gpu_linalg3_support, only: run_linalg3_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$omp target map(from: on_device)
    on_device = .not. omp_is_initial_device()
    !$omp end target
    if (.not. on_device) error stop "OpenMP linalg3 did not run on target"
    call run_linalg3_test()
end program test_openmp_linalg3_offload
