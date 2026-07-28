module openmp_enzyme_fir_kernel
    use enzyme, only: enzyme_dup, enzyme_fwddiff
    implicit none
contains
    subroutine square(x, y)
        !$omp declare target
        real, intent(in) :: x
        real, intent(out) :: y

        y = x*x
    end subroutine square

    subroutine differentiated_square(x, dx, y, dy)
        !$omp declare target
        real, intent(in) :: x, dx
        real, intent(out) :: y, dy

        call enzyme_fwddiff(square, enzyme_dup, x, dx, enzyme_dup, y, dy)
    end subroutine differentiated_square
end module openmp_enzyme_fir_kernel

program test_openmp_enzyme_fir_offload
    use omp_lib, only: omp_is_initial_device
    use openmp_enzyme_fir_kernel, only: differentiated_square
    implicit none

    real :: x, dx, y, dy
    logical :: executed_on_device

    x = 3.0
    dx = 1.0
    y = 0.0
    dy = 0.0
    executed_on_device = .false.

    !$omp target map(to: x, dx) map(from: y, dy, executed_on_device)
    call differentiated_square(x, dx, y, dy)
    executed_on_device = .not. omp_is_initial_device()
    !$omp end target

    if (.not. executed_on_device) then
        error stop "Enzyme OpenMP target ran on the initial device"
    end if
    if (abs(y - 9.0) > 1.0e-6) then
        error stop "Enzyme OpenMP primal disagrees with analytical oracle"
    end if
    if (abs(dy - 6.0) > 1.0e-6) then
        error stop "Enzyme OpenMP JVP disagrees with analytical oracle"
    end if
end program test_openmp_enzyme_fir_offload
