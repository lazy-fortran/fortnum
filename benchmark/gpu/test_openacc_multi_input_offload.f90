program test_openacc_multi_input_offload
    use openacc, only: acc_device_nvidia, acc_on_device
    use test_gpu_multi_input_support, only: run_multi_input_tests
    implicit none

    logical :: executed_on_nvidia

    executed_on_nvidia = .false.
    !$acc serial copyout(executed_on_nvidia)
    executed_on_nvidia = acc_on_device(acc_device_nvidia)
    !$acc end serial
    if (.not. executed_on_nvidia) then
        error stop "OpenACC region did not execute on an NVIDIA device"
    end if
    call run_multi_input_tests()
end program test_openacc_multi_input_offload
