program test_cpu_fft8_batch
    use test_gpu_fft8_support, only: run_fft8_test
    implicit none

    call run_fft8_test()
    write (*, "(a)") "PASS"
end program test_cpu_fft8_batch
