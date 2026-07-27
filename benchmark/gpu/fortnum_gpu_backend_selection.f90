module fortnum_gpu_backend_selection
    implicit none
    private

    integer, parameter :: selection_count = 12
    character(3), parameter :: products(selection_count) = [ &
        "jvp", "jvp", "jvp", "jvp", "jvp", "jvp", &
        "vjp", "vjp", "vjp", "vjp", "vjp", "vjp"]
    integer, parameter :: batch_sizes(selection_count) = [ &
        256, 65536, 1048576, 256, 65536, 1048576, &
        256, 65536, 1048576, 256, 65536, 1048576]
    logical, parameter :: resident(selection_count) = [ &
        .true., .true., .true., .false., .false., .false., &
        .true., .true., .true., .false., .false., .false.]
    character(7), parameter :: backends(selection_count) = [ &
        "openacc", "openacc", "openacc", "openacc", "openacc", "openacc", &
        "openacc", "openacc", "openacc", "openacc", "openacc", "openacc"]

    public :: select_multi_input_gpu_backend

contains

    pure subroutine select_multi_input_gpu_backend( &
            product, batch_size, is_resident, backend, found)
        character(*), intent(in) :: product
        integer, intent(in) :: batch_size
        logical, intent(in) :: is_resident
        character(*), intent(out) :: backend
        logical, intent(out) :: found
        integer :: i

        backend = ""
        found = .false.
        do i = 1, selection_count
            if (trim(product) /= products(i)) cycle
            if (batch_size /= batch_sizes(i)) cycle
            if (is_resident .neqv. resident(i)) cycle
            backend = backends(i)
            found = .true.
            return
        end do
    end subroutine select_multi_input_gpu_backend

end module fortnum_gpu_backend_selection
