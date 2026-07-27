program plot_gpu_report
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortplot, only: bar, figure, legend, plot, savefig_with_status, &
        set_xscale, set_yscale, set_xticks, title, xlabel, ylabel
    implicit none

    integer, parameter :: batch_rows = 3, active_rows = 4, profile_rows = 4
    real(dp) :: batch(batch_rows), batch_values(batch_rows, 10)
    real(dp) :: active(active_rows), active_values(active_rows, 6)
    real(dp) :: profile_values(profile_rows, 5)
    character(len=16) :: profile_labels(profile_rows)
    character(len=512) :: batch_path, active_path, profile_path, output_dir

    call require_arguments(batch_path, active_path, profile_path, output_dir)
    call read_batch_data(trim(batch_path), batch, batch_values)
    call read_active_data(trim(active_path), active, active_values)
    call read_profile_data(trim(profile_path), profile_labels, profile_values)
    call plot_batch_product( &
        batch, batch_values(:, 1:5), 'JVP', &
        trim(output_dir)//'/gpu_jvp_batch_scaling.png')
    call plot_batch_product( &
        batch, batch_values(:, 6:10), 'VJP', &
        trim(output_dir)//'/gpu_vjp_batch_scaling.png')
    call plot_active_product( &
        active, active_values(:, 1:3), 'JVP', &
        trim(output_dir)//'/gpu_jvp_active_scaling.png')
    call plot_active_product( &
        active, active_values(:, 4:6), 'VJP', &
        trim(output_dir)//'/gpu_vjp_active_scaling.png')
    call plot_profile_bars( &
        profile_labels, profile_values(:, 1), &
        'Achieved device-memory bandwidth', 'bandwidth (GB/s)', &
        trim(output_dir)//'/gpu_profile_bandwidth.png')
    call plot_profile_bars( &
        profile_labels, profile_values(:, 5), &
        'Peak device allocation', 'device allocation (MB)', &
        trim(output_dir)//'/gpu_peak_device_memory.png')

contains

    subroutine require_arguments(batch_file, active_file, profile_file, output)
        character(len=*), intent(out) :: batch_file, active_file, profile_file
        character(len=*), intent(out) :: output

        if (command_argument_count() /= 4) then
            error stop 'usage: plot_gpu_report BATCH.csv ACTIVE.csv PROFILE.csv OUTPUT_DIR'
        end if
        call get_command_argument(1, batch_file)
        call get_command_argument(2, active_file)
        call get_command_argument(3, profile_file)
        call get_command_argument(4, output)
    end subroutine require_arguments

    subroutine read_batch_data(path, x, values)
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: x(:), values(:, :)
        character(len=1024) :: header
        integer :: i, unit, stat

        open (newunit=unit, file=path, status='old', action='read', iostat=stat)
        if (stat /= 0) error stop 'cannot open GPU batch data'
        read (unit, '(a)', iostat=stat) header
        if (stat /= 0) error stop 'cannot read GPU batch header'
        do i = 1, size(x)
            read (unit, *, iostat=stat) x(i), values(i, :)
            if (stat /= 0) error stop 'invalid GPU batch row'
        end do
        close (unit)
    end subroutine read_batch_data

    subroutine read_active_data(path, x, values)
        character(len=*), intent(in) :: path
        real(dp), intent(out) :: x(:), values(:, :)
        character(len=1024) :: header
        integer :: i, unit, stat

        open (newunit=unit, file=path, status='old', action='read', iostat=stat)
        if (stat /= 0) error stop 'cannot open GPU active-input data'
        read (unit, '(a)', iostat=stat) header
        if (stat /= 0) error stop 'cannot read GPU active-input header'
        do i = 1, size(x)
            read (unit, *, iostat=stat) x(i), values(i, :)
            if (stat /= 0) error stop 'invalid GPU active-input row'
        end do
        close (unit)
    end subroutine read_active_data

    subroutine read_profile_data(path, labels, values)
        character(len=*), intent(in) :: path
        character(len=*), intent(out) :: labels(:)
        real(dp), intent(out) :: values(:, :)
        character(len=1024) :: header
        integer :: i, unit, stat

        open (newunit=unit, file=path, status='old', action='read', iostat=stat)
        if (stat /= 0) error stop 'cannot open GPU profile data'
        read (unit, '(a)', iostat=stat) header
        if (stat /= 0) error stop 'cannot read GPU profile header'
        do i = 1, size(labels)
            read (unit, *, iostat=stat) labels(i), values(i, :)
            if (stat /= 0) error stop 'invalid GPU profile row'
        end do
        close (unit)
    end subroutine read_profile_data

    subroutine plot_batch_product(x, values, product, path)
        real(dp), contiguous, intent(in) :: x(:), values(:, :)
        character(len=*), intent(in) :: product, path
        real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
        real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
        real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp

        call figure(figsize=[9.0_dp, 6.2_dp])
        call plot(x, values(:, 1), label='CPU', color=blue, &
            linestyle='-', marker='o', linewidth=2.2_dp)
        call plot(x, values(:, 2), label='OpenACC resident', color=orange, &
            linestyle='-', marker='s', linewidth=2.0_dp)
        call plot(x, values(:, 3), label='OpenACC transfer', color=orange, &
            linestyle='--', marker='x', linewidth=2.0_dp)
        call plot(x, values(:, 4), label='OpenMP resident', color=green, &
            linestyle='-', marker='d', linewidth=2.0_dp)
        call plot(x, values(:, 5), label='OpenMP transfer', color=green, &
            linestyle=':', marker='+', linewidth=2.0_dp)
        call set_xscale('log')
        call set_yscale('log')
        call title('Generated analytical '//trim(product)// &
            ': CPU/GPU crossover')
        call xlabel('batch size (points, log scale)')
        call ylabel('complete wall clock (ms, log scale)')
        call legend()
        call save_checked(path)
    end subroutine plot_batch_product

    subroutine plot_active_product(x, values, product, path)
        real(dp), contiguous, intent(in) :: x(:), values(:, :)
        character(len=*), intent(in) :: product, path
        real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
        real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
        real(dp), parameter :: green(3) = [0.0_dp, 158.0_dp, 115.0_dp]/255.0_dp

        call figure(figsize=[9.0_dp, 6.2_dp])
        call plot(x, values(:, 1), label='CPU', color=blue, &
            linestyle='-', marker='o', linewidth=2.2_dp)
        call plot(x, values(:, 2), label='OpenACC resident', color=orange, &
            linestyle='--', marker='s', linewidth=2.0_dp)
        call plot(x, values(:, 3), label='OpenMP resident', color=green, &
            linestyle=':', marker='d', linewidth=2.0_dp)
        call set_xscale('log', base=2.0_dp)
        call set_yscale('log')
        call title('Generated analytical '//trim(product)// &
            ': active-input scaling')
        call xlabel('active inputs (count, base-2 log scale)')
        call ylabel('wall clock at 65,536 points (ms, log scale)')
        call legend()
        call save_checked(path)
    end subroutine plot_active_product

    subroutine plot_profile_bars(labels, values, heading, y_label, path)
        character(len=*), intent(in) :: labels(:)
        real(dp), intent(in) :: values(:)
        character(len=*), intent(in) :: heading, y_label, path
        real(dp) :: x(size(values)), one_x(1), one_y(1)
        real(dp), parameter :: colors(3, 4) = reshape([ &
            230.0_dp, 159.0_dp, 0.0_dp, &
            230.0_dp, 159.0_dp, 0.0_dp, &
            0.0_dp, 158.0_dp, 115.0_dp, &
            0.0_dp, 158.0_dp, 115.0_dp], [3, 4])/255.0_dp
        integer :: i

        do i = 1, size(values)
            x(i) = real(i, dp)
        end do
        call figure(figsize=[9.0_dp, 6.2_dp])
        do i = 1, size(values)
            one_x(1) = x(i)
            one_y(1) = values(i)
            call bar(one_x, one_y, color=colors(:, i), width=0.65_dp)
        end do
        call set_xticks(x, labels)
        call title(heading)
        call ylabel(y_label)
        call save_checked(path)
    end subroutine plot_profile_bars

    subroutine save_checked(path)
        character(len=*), intent(in) :: path
        integer :: stat

        call savefig_with_status(path, stat)
        if (stat /= 0) error stop 'fortplot could not write GPU figure'
    end subroutine save_checked

end program plot_gpu_report
