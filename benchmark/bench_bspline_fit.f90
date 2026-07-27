program bench_bspline_fit
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_benchmark_memory, only: peak_rss_bytes
    use fortnum_bspline, only: bspline_workspace_t, bspline_init, &
        bspline_set_knots, bspline_eval_basis, bspline_fit_jvp_factored, &
        bspline_fit_vjp_factored
    use fortnum_linalg, only: lu_factor, lu_solve, LINALG_OK
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

    integer, parameter :: samples = 15
    real(dp), parameter :: h = 1.0e-5_dp
    real(dp), allocatable :: basis(:, :), factors(:, :), transpose_factors(:, :)
    real(dp), allocatable :: dbasis(:, :), basis_bar(:, :)
    real(dp), allocatable :: values(:), dvalues(:), coef(:), dcoef(:)
    real(dp), allocatable :: cbar(:), values_bar(:)
    integer, allocatable :: pivots(:), transpose_pivots(:)
    integer :: n, info, sample
    integer(int64) :: reps
    character(16) :: candidate, product, size_arg, mode
    logical :: memory_only
    real(dp) :: warmup

    call get_command_argument(1, candidate)
    call get_command_argument(2, product)
    call get_command_argument(3, size_arg)
    call get_command_argument(4, mode)
    read (size_arg, *) n
    memory_only = trim(mode) == "--peak-rss"
    if ((trim(candidate) /= "analytical") .and. &
        (trim(candidate) /= "diagnostic")) error stop "invalid candidate"
    if ((trim(product) /= "jvp") .and. &
        (trim(product) /= "vjp")) error stop "invalid product"
    if ((n /= 4) .and. (n /= 8) .and. (n /= 16)) error stop "n must be 4, 8, or 16"

    call initialize_inputs()
    if (trim(candidate) == "analytical") then
        reps = merge(200000_int64, 100000_int64, trim(product) == "jvp")
    else
        reps = merge(20000_int64, 200_int64, trim(product) == "jvp")
    end if
    if (memory_only) then
        warmup = run_sample(max(1_int64, reps/10_int64))
        write (*, "(i0)") peak_rss_bytes()
        stop
    end if
    do sample = 1, 3
        warmup = run_sample(max(1_int64, reps/10_int64))
    end do
    do sample = 1, samples
        write (*, "(f0.4)") run_sample(reps)
    end do

contains

    subroutine initialize_inputs()
        type(bspline_workspace_t) :: ws
        type(fortnum_status_t) :: status
        real(dp), allocatable :: breakpoints(:), row(:)
        integer :: i, j, nbreak

        nbreak = n - 2
        allocate (basis(n, n), factors(n, n), transpose_factors(n, n))
        allocate (dbasis(n, n), basis_bar(n, n))
        allocate (values(n), dvalues(n), coef(n), dcoef(n), cbar(n), values_bar(n))
        allocate (pivots(n), transpose_pivots(n), breakpoints(nbreak), row(n))
        do i = 1, nbreak
            breakpoints(i) = real(i - 1, dp)/real(nbreak - 1, dp)
        end do
        call bspline_init(ws, 4, nbreak, status)
        call bspline_set_knots(ws, breakpoints, status)
        if (status%code /= FORTNUM_OK) error stop "spline setup failed"
        do i = 1, n
            call bspline_eval_basis(ws, real(i - 1, dp)/real(n - 1, dp), &
                row, status)
            basis(i, :) = row
            coef(i) = sin(0.4_dp*real(i, dp))
            dvalues(i) = 0.03_dp*real(mod(5*i, 7) - 3, dp)
            cbar(i) = 0.05_dp*real(mod(3*i, 11) - 5, dp)
            do j = 1, n
                dbasis(i, j) = 0.002_dp*real(mod(3*i + 5*j, 13) - 6, dp)
            end do
        end do
        values = matmul(basis, coef)
        factors = basis
        call lu_factor(n, factors, pivots, info)
        transpose_factors = transpose(basis)
        call lu_factor(n, transpose_factors, transpose_pivots, info)
        if (info /= LINALG_OK) error stop "collocation factorization failed"
    end subroutine initialize_inputs

    function run_sample(count) result(ns_per_call)
        integer(int64), intent(in) :: count
        real(dp) :: ns_per_call, sink
        integer(int64) :: iteration, tick0, tick1, rate

        sink = 0.0_dp
        call system_clock(tick0, rate)
        do iteration = 1_int64, count
            dvalues(1) = 0.01_dp*real(mod(iteration, 17_int64) - 8_int64, dp)
            if (trim(product) == "jvp") then
                call run_jvp()
                sink = sink + dcoef(1)
            else
                call run_vjp()
                sink = sink + basis_bar(1, 1) + values_bar(1)
            end if
        end do
        call system_clock(tick1)
        if (sink /= sink) error stop "benchmark failed"
        ns_per_call = real(tick1 - tick0, dp)*1.0e9_dp/ &
            (real(rate, dp)*real(count, dp))
    end function run_sample

    subroutine run_jvp()
        type(fortnum_status_t) :: status
        real(dp) :: plus_a(n, n), minus_a(n, n)
        real(dp) :: plus_b(n), minus_b(n)

        if (trim(candidate) == "analytical") then
            call bspline_fit_jvp_factored(factors, pivots, coef, dbasis, &
                dvalues, dcoef, status)
            return
        end if
        plus_a = basis + h*dbasis
        minus_a = basis - h*dbasis
        plus_b = values + h*dvalues
        minus_b = values - h*dvalues
        call lu_solve(n, plus_a, plus_b, info)
        call lu_solve(n, minus_a, minus_b, info)
        dcoef = (plus_b - minus_b)/(2.0_dp*h)
    end subroutine run_jvp

    subroutine run_vjp()
        type(fortnum_status_t) :: status
        real(dp) :: work_a(n, n), work_b(n), plus, minus
        integer :: i, j

        if (trim(candidate) == "analytical") then
            call bspline_fit_vjp_factored(transpose_factors, transpose_pivots, &
                coef, cbar, basis_bar, values_bar, status)
            return
        end if
        do j = 1, n
            do i = 1, n
                work_a = basis
                work_a(i, j) = work_a(i, j) + h
                work_b = values
                call lu_solve(n, work_a, work_b, info)
                plus = dot_product(cbar, work_b)
                work_a = basis
                work_a(i, j) = work_a(i, j) - h
                work_b = values
                call lu_solve(n, work_a, work_b, info)
                minus = dot_product(cbar, work_b)
                basis_bar(i, j) = (plus - minus)/(2.0_dp*h)
            end do
            work_a = basis
            work_b = values
            work_b(j) = work_b(j) + h
            call lu_solve(n, work_a, work_b, info)
            plus = dot_product(cbar, work_b)
            work_a = basis
            work_b = values
            work_b(j) = work_b(j) - h
            call lu_solve(n, work_a, work_b, info)
            minus = dot_product(cbar, work_b)
            values_bar(j) = (plus - minus)/(2.0_dp*h)
        end do
    end subroutine run_vjp

end program bench_bspline_fit
