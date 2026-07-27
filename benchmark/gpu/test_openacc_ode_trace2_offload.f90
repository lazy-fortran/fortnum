program test_openacc_ode_trace2_offload
    use openacc, only: acc_device_nvidia, acc_on_device
    use test_gpu_ode_trace2_support, only: run_ode_trace2_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$acc serial copyout(on_device)
    on_device = acc_on_device(acc_device_nvidia)
    !$acc end serial
    if (.not. on_device) error stop "OpenACC ODE trace did not run on NVIDIA"
    call run_ode_trace2_test()
end program test_openacc_ode_trace2_offload
