program test_openacc_implicit_root_offload
    use openacc, only: acc_device_nvidia, acc_on_device
    use test_gpu_implicit_root_support, only: run_implicit_root_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$acc serial copyout(on_device)
    on_device = acc_on_device(acc_device_nvidia)
    !$acc end serial
    if (.not. on_device) error stop "OpenACC implicit root did not run on NVIDIA"
    call run_implicit_root_test()
end program test_openacc_implicit_root_offload
