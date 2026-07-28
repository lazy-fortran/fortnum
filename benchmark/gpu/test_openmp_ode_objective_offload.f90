program test_openmp_ode_objective_offload
    use omp_lib, only: omp_is_initial_device
    use test_gpu_ode_objective_support, only: run_ode_objective_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$omp target map(from: on_device)
    on_device = .not. omp_is_initial_device()
    !$omp end target
    if (.not. on_device) error stop "OpenMP ODE objective did not run on target"
    call run_ode_objective_test()
end program test_openmp_ode_objective_offload
