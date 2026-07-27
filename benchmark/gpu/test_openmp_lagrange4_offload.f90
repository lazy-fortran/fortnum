program test_openmp_lagrange4_offload
    use omp_lib, only: omp_is_initial_device
    use test_gpu_lagrange4_support, only: run_lagrange4_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$omp target map(from: on_device)
    on_device = .not. omp_is_initial_device()
    !$omp end target
    if (.not. on_device) error stop "OpenMP Lagrange did not run on target"
    call run_lagrange4_test()
end program test_openmp_lagrange4_offload
