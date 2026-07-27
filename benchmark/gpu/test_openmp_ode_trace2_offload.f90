program test_openmp_ode_trace2_offload
    use omp_lib, only: omp_is_initial_device
    use test_gpu_ode_trace2_support, only: run_ode_trace2_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$omp target map(from: on_device)
    on_device = .not. omp_is_initial_device()
    !$omp end target
    if (.not. on_device) error stop "OpenMP ODE trace did not run on target"
    call run_ode_trace2_test()
end program test_openmp_ode_trace2_offload
