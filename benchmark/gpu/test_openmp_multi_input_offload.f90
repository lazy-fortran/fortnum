program test_openmp_multi_input_offload
    use omp_lib, only: omp_is_initial_device
    use test_gpu_multi_input_support, only: run_multi_input_tests
    implicit none

    logical :: executed_on_device

    executed_on_device = .false.
    !$omp target map(from: executed_on_device)
    executed_on_device = .not. omp_is_initial_device()
    !$omp end target
    if (.not. executed_on_device) then
        error stop "OpenMP target region executed on the initial device"
    end if
    call run_multi_input_tests()
end program test_openmp_multi_input_offload
