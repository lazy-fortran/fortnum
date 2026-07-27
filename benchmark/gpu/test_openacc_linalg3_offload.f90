program test_openacc_linalg3_offload
    use openacc, only: acc_device_nvidia, acc_on_device
    use test_gpu_linalg3_support, only: run_linalg3_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$acc serial copyout(on_device)
    on_device = acc_on_device(acc_device_nvidia)
    !$acc end serial
    if (.not. on_device) error stop "OpenACC linalg3 did not run on NVIDIA"
    call run_linalg3_test()
end program test_openacc_linalg3_offload
