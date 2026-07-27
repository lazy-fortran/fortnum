program test_openmp_implicit_root_offload
    use omp_lib, only: omp_is_initial_device
    use test_gpu_implicit_root_support, only: run_implicit_root_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$omp target map(from: on_device)
    on_device = .not. omp_is_initial_device()
    !$omp end target
    if (.not. on_device) error stop "OpenMP implicit root did not run on target"
    call run_implicit_root_test()
end program test_openmp_implicit_root_offload
