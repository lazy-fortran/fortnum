program test_openacc_rk54_offload
    use openacc, only: acc_device_nvidia, acc_on_device
    use test_gpu_rk54_support, only: run_rk54_device_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$acc serial copyout(on_device)
    on_device = acc_on_device(acc_device_nvidia)
    !$acc end serial
    if (.not. on_device) error stop "OpenACC RK54 did not run on NVIDIA"
    call run_rk54_device_test()
end program test_openacc_rk54_offload
