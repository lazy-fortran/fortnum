program test_openacc_dawson_variants_offload
    use openacc, only: acc_device_nvidia, acc_on_device
    use test_gpu_dawson_variants_support, only: run_dawson_variants_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$acc serial copyout(on_device)
    on_device = acc_on_device(acc_device_nvidia)
    !$acc end serial
    if (.not. on_device) error stop "OpenACC variants did not run on NVIDIA"
    call run_dawson_variants_test()
end program test_openacc_dawson_variants_offload
