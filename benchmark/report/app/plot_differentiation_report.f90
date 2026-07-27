program plot_differentiation_report
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortplot, only: bar, figure, hist, legend, plot, savefig_with_status, &
        set_xticks, title, xlabel, ylabel
    implicit none

    integer, parameter :: max_rows = 128
    integer, parameter :: mechanism_count = 4
    character(len=16), parameter :: mechanisms(mechanism_count) = [ &
        character(len=16) :: 'analytical', 'autodiff', 'hybrid', 'diagnostic']
    character(len=16) :: selected(max_rows), fastest(max_rows)
    real(dp) :: fastest_ns(max_rows), runner_up_ns(max_rows)
    integer :: counts(mechanism_count)
    integer :: n_rows
    character(len=512) :: data_path, output_dir

    call require_arguments(data_path, output_dir)
    call read_tournaments(trim(data_path), selected, fastest, fastest_ns, &
        runner_up_ns, n_rows)
    call count_mechanisms(selected(:n_rows), mechanisms, counts)
    call plot_wins(counts, mechanisms, trim(output_dir)//'/mechanism_wins.png')
    call plot_slowdown(fastest_ns(:n_rows), runner_up_ns(:n_rows), &
        trim(output_dir)//'/runner_up_slowdown.png')
    call plot_direct_solver_scaling(trim(output_dir))

contains

    subroutine require_arguments(data_file, output_path)
        character(len=*), intent(out) :: data_file, output_path

        if (command_argument_count() /= 2) then
            error stop 'usage: plot_differentiation_report DATA.csv OUTPUT_DIR'
        end if
        call get_command_argument(1, data_file)
        call get_command_argument(2, output_path)
    end subroutine require_arguments

    subroutine read_tournaments(path, selected_values, fastest_values, &
            fastest_times, runner_up_times, n)
        character(len=*), intent(in) :: path
        character(len=*), intent(out) :: selected_values(:), fastest_values(:)
        real(dp), intent(out) :: fastest_times(:), runner_up_times(:)
        integer, intent(out) :: n
        character(len=512) :: line
        character(len=96) :: workload, source_record
        integer :: unit, stat

        open (newunit=unit, file=path, status='old', action='read', iostat=stat)
        if (stat /= 0) error stop 'cannot open tournament data'
        read (unit, '(a)', iostat=stat) line
        if (stat /= 0) error stop 'cannot read tournament header'

        n = 0
        do
            read (unit, '(a)', iostat=stat) line
            if (stat < 0) exit
            if (stat > 0) error stop 'cannot read tournament row'
            if (len_trim(line) == 0) cycle
            n = n + 1
            if (n > size(selected_values)) error stop 'too many tournament rows'
            read (line, *, iostat=stat) workload, selected_values(n), &
                fastest_values(n), fastest_times(n), runner_up_times(n), source_record
            if (stat /= 0) error stop 'invalid tournament row'
            if (fastest_times(n) <= 0.0_dp) error stop 'non-positive fastest time'
            if (runner_up_times(n) < fastest_times(n)) then
                error stop 'runner-up time is smaller than fastest time'
            end if
        end do
        close (unit)
        if (n == 0) error stop 'tournament data is empty'
    end subroutine read_tournaments

    subroutine count_mechanisms(values, names, result)
        character(len=*), intent(in) :: values(:), names(:)
        integer, intent(out) :: result(:)
        integer :: i, j
        logical :: matched

        result = 0
        do i = 1, size(values)
            matched = .false.
            do j = 1, size(names)
                if (trim(values(i)) == trim(names(j))) then
                    result(j) = result(j) + 1
                    matched = .true.
                    exit
                end if
            end do
            if (.not. matched) error stop 'unknown selected mechanism'
        end do
    end subroutine count_mechanisms

    subroutine plot_wins(values, labels, path)
        integer, intent(in) :: values(:)
        character(len=*), intent(in) :: labels(:)
        character(len=*), intent(in) :: path
        real(dp) :: x(size(values)), y(1), one_x(1)
        integer :: i
        real(dp), parameter :: colors(3, 4) = reshape([ &
            0.0_dp, 114.0_dp, 178.0_dp, &
            230.0_dp, 159.0_dp, 0.0_dp, &
            0.0_dp, 158.0_dp, 115.0_dp, &
            204.0_dp, 121.0_dp, 167.0_dp], [3, 4])/255.0_dp

        do i = 1, size(values)
            x(i) = real(i, dp)
        end do
        call figure()
        do i = 1, size(values)
            one_x(1) = x(i)
            y(1) = real(values(i), dp)
            call bar(one_x, y, color=colors(:, i), width=0.7_dp, &
                label=trim(labels(i)))
        end do
        call set_xticks(x, labels)
        call title('Selected mechanism across measured tournaments')
        call ylabel('workloads selected (count)')
        call save_checked(path)
    end subroutine plot_wins

    subroutine plot_slowdown(best, second, path)
        real(dp), intent(in) :: best(:), second(:)
        character(len=*), intent(in) :: path
        real(dp) :: slowdown(size(best))

        slowdown = log10(second/best)
        call figure()
        call hist(slowdown, bins=8, color='#0072B2', alpha=0.85_dp)
        call title('Separation between fastest and second-fastest candidates')
        call xlabel('log10(second-fastest / fastest wall clock)')
        call ylabel('workloads (count)')
        call save_checked(path)
    end subroutine plot_slowdown

    subroutine plot_direct_solver_scaling(output_path)
        character(len=*), intent(in) :: output_path
        real(dp), parameter :: products(3) = [1.0_dp, 4.0_dp, 16.0_dp]
        real(dp), parameter :: jvp_analytical(3) = [55.3112_dp, 207.7022_dp, &
            773.70575_dp]
        real(dp), parameter :: jvp_autodiff(3) = [189.42425_dp, 720.94535_dp, &
            2719.0076_dp]
        real(dp), parameter :: jvp_diagnostic(3) = [72.16495_dp, 257.4344_dp, &
            982.2762_dp]
        real(dp), parameter :: vjp_analytical(3) = [59.397_dp, 227.0375_dp, &
            852.922_dp]
        real(dp), parameter :: vjp_autodiff(3) = [73.5035_dp, 307.519_dp, &
            1089.863_dp]
        real(dp), parameter :: vjp_diagnostic(3) = [1025.3865_dp, 3972.388_dp, &
            15454.5425_dp]

        call scaling_figure(products, jvp_analytical, jvp_autodiff, &
            jvp_diagnostic, 'Direct-solver forward JVP scaling', &
            'JVP directions (count)', &
            trim(output_path)//'/direct_solver_jvp_scaling.png')
        call scaling_figure(products, vjp_analytical, vjp_autodiff, &
            vjp_diagnostic, 'Direct-solver reverse VJP scaling', &
            'VJP cotangents (count)', &
            trim(output_path)//'/direct_solver_vjp_scaling.png')
    end subroutine plot_direct_solver_scaling

    subroutine scaling_figure(x, analytical, autodiff, diagnostic, heading, &
            x_label, path)
        real(dp), contiguous, intent(in) :: x(:)
        real(dp), intent(in) :: analytical(:), autodiff(:), diagnostic(:)
        character(len=*), intent(in) :: heading, x_label, path
        real(dp), parameter :: blue(3) = [0.0_dp, 114.0_dp, 178.0_dp]/255.0_dp
        real(dp), parameter :: orange(3) = [230.0_dp, 159.0_dp, 0.0_dp]/255.0_dp
        real(dp), parameter :: purple(3) = [204.0_dp, 121.0_dp, 167.0_dp]/255.0_dp
        real(dp) :: analytical_us(size(x)), autodiff_us(size(x))
        real(dp) :: diagnostic_us(size(x))

        analytical_us = analytical/1000.0_dp
        autodiff_us = autodiff/1000.0_dp
        diagnostic_us = diagnostic/1000.0_dp
        call figure()
        call plot(x, analytical_us, label='analytical', color=blue, &
            linestyle='-', marker='o', linewidth=2.0_dp)
        call plot(x, autodiff_us, label='autodiff', color=orange, &
            linestyle='--', marker='s', linewidth=2.0_dp)
        call plot(x, diagnostic_us, label='finite-difference diagnostic', &
            color=purple, linestyle=':', marker='x', linewidth=2.0_dp)
        call title(heading)
        call xlabel(x_label)
        call ylabel('complete product wall clock (microseconds)')
        call legend()
        call save_checked(path)
    end subroutine scaling_figure

    subroutine save_checked(path)
        character(len=*), intent(in) :: path
        integer :: stat

        call savefig_with_status(path, stat)
        if (stat /= 0) error stop 'fortplot could not write figure'
    end subroutine save_checked

end program plot_differentiation_report
