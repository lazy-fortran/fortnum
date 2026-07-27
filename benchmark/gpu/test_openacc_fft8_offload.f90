program test_openacc_fft8_offload
    use openacc, only: acc_device_nvidia, acc_on_device
    use test_gpu_fft8_support, only: run_fft8_test
    implicit none
    logical :: on_device

    on_device = .false.
    !$acc serial copyout(on_device)
    on_device = acc_on_device(acc_device_nvidia)
    !$acc end serial
    if (.not. on_device) error stop "OpenACC FFT8 did not run on NVIDIA"
    call run_fft8_test()
end program test_openacc_fft8_offload
