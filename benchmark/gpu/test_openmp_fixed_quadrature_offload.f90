program test_openmp_fixed_quadrature_offload
    use omp_lib, only: omp_is_initial_device
    use test_gpu_fixed_quadrature_support, only: run_fixed_quadrature_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$omp target map(from: on_device)
    on_device = .not. omp_is_initial_device()
    !$omp end target
    if (.not. on_device) error stop "OpenMP quadrature did not run on target"
    call run_fixed_quadrature_test()
end program test_openmp_fixed_quadrature_offload
