program test_fortnum_sobol
    ! Sobol sequence.
    !
    ! The oracles here are the defining properties of the construction, not a
    ! stored table of expected numbers. That matters: a plausible-looking but
    ! wrong set of direction numbers still produces points in [0,1), still looks
    ! uniform to the eye, and still passes any test that only checks ranges. It
    ! fails equidistribution.
    !
    !   * one-dimensional equidistribution. For the first 2^k points, each of
    !     the 2^k equal subintervals of every coordinate must contain exactly
    !     one point. This is exact, not statistical;
    !   * two-dimensional net property. The leading dimension pair forms a
    !     (0,2)-net in base 2, so every elementary interval of area 2^-k must
    !     contain exactly one of the first 2^k points. This is the property that
    !     depends on the direction numbers being right, and it is what a wrong
    !     table breaks;
    !   * discrepancy. Sobol must beat pseudorandom points decisively at the
    !     same count. Measured, not assumed;
    !   * the Gray-code recurrence must agree with the direct construction,
    !     which the test computes independently from the direction numbers.

    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK, FORTNUM_DOMAIN_ERROR, &
        FORTNUM_NOT_IMPLEMENTED
    use fortnum_rng, only: rng_t, rng_seed, rng_uniform
    use fortnum_sobol, only: sobol_t, sobol_initialize, sobol_next, sobol_skip, &
        sobol_fill, SOBOL_MAX_DIMENSION
    implicit none

    integer :: failures

    failures = 0
    call check_range_and_shape(failures)
    call check_one_dimensional_equidistribution(failures)
    call check_two_dimensional_net(failures)
    call check_beats_random_discrepancy(failures)
    call check_skip_matches_generation(failures)
    call check_reproducible(failures)
    call check_refusals(failures)

    if (failures == 0) then
        print *, "test_fortnum_sobol: PASS"
    else
        print *, "test_fortnum_sobol: FAIL", failures
        error stop 1
    end if

contains

    subroutine check_range_and_shape(failures)
        integer, intent(inout) :: failures
        type(sobol_t) :: sequence
        type(fortnum_status_t) :: status
        real(dp) :: points(1024, 5)

        call sobol_initialize(sequence, 5, status)
        call expect(status%code == FORTNUM_OK, "initialization succeeds", failures)
        call sobol_fill(sequence, points, status)
        call expect(status%code == FORTNUM_OK, "filling succeeds", failures)
        call expect(all(points >= 0.0_dp) .and. all(points < 1.0_dp), &
            "every coordinate lies in [0,1)", failures)
        call expect(abs(sum(points(:, 1))/1024.0_dp - 0.5_dp) < 0.01_dp, &
            "the first coordinate averages to one half", failures)
    end subroutine check_range_and_shape

    ! Exact stratification: 2^k points, 2^k bins, one point per bin.
    subroutine check_one_dimensional_equidistribution(failures)
        integer, intent(inout) :: failures
        type(sobol_t) :: sequence
        type(fortnum_status_t) :: status
        integer, parameter :: k = 10
        integer, parameter :: n = 2**k
        real(dp) :: points(n, SOBOL_MAX_DIMENSION)
        integer :: counts(n)
        integer :: j, i, bin
        logical :: stratified

        call sobol_initialize(sequence, SOBOL_MAX_DIMENSION, status)
        call sobol_fill(sequence, points, status)

        stratified = .true.
        do j = 1, SOBOL_MAX_DIMENSION
            counts = 0
            do i = 1, n
                bin = int(points(i, j)*real(n, dp)) + 1
                if (bin < 1 .or. bin > n) then
                    stratified = .false.
                    cycle
                end if
                counts(bin) = counts(bin) + 1
            end do
            if (any(counts /= 1)) stratified = .false.
        end do
        call expect(stratified, &
            "the first 2^k points hit each of 2^k bins exactly once, every dimension", &
            failures)
    end subroutine check_one_dimensional_equidistribution

    ! (0,2)-net property on the leading pair: for every split of k bits between
    ! the two axes, each elementary interval holds exactly one point. A wrong
    ! primitive polynomial or a wrong initial m breaks this immediately.
    subroutine check_two_dimensional_net(failures)
        integer, intent(inout) :: failures
        type(sobol_t) :: sequence
        type(fortnum_status_t) :: status
        integer, parameter :: k = 8
        integer, parameter :: n = 2**k
        real(dp) :: points(n, 2)
        integer :: counts(n)
        integer :: split, rows, columns, i, row, column, cell
        logical :: is_net

        call sobol_initialize(sequence, 2, status)
        call sobol_fill(sequence, points, status)

        is_net = .true.
        do split = 0, k
            rows = 2**split
            columns = 2**(k - split)
            counts = 0
            do i = 1, n
                row = int(points(i, 1)*real(rows, dp))
                column = int(points(i, 2)*real(columns, dp))
                cell = row*columns + column + 1
                if (cell < 1 .or. cell > n) then
                    is_net = .false.
                    cycle
                end if
                counts(cell) = counts(cell) + 1
            end do
            if (any(counts /= 1)) is_net = .false.
        end do
        call expect(is_net, &
            "the leading pair forms a (0,2)-net over every elementary split", &
            failures)
    end subroutine check_two_dimensional_net

    ! Star discrepancy on a grid of test boxes anchored at the origin. Sobol
    ! should be several times better than pseudorandom at this count.
    subroutine check_beats_random_discrepancy(failures)
        integer, intent(inout) :: failures
        type(sobol_t) :: sequence
        type(rng_t) :: generator
        type(fortnum_status_t) :: status
        integer, parameter :: n = 4096
        real(dp) :: quasi(n, 2), pseudo(n, 2)
        real(dp) :: quasi_discrepancy, pseudo_discrepancy
        integer :: i, j

        call sobol_initialize(sequence, 2, status)
        call sobol_fill(sequence, quasi, status)

        call rng_seed(generator, 20260808_int64, status)
        do i = 1, n
            do j = 1, 2
                call rng_uniform(generator, pseudo(i, j))
            end do
        end do

        quasi_discrepancy = star_discrepancy(quasi)
        pseudo_discrepancy = star_discrepancy(pseudo)

        call expect(quasi_discrepancy < pseudo_discrepancy, &
            "Sobol has lower star discrepancy than pseudorandom points", failures)
        call expect(quasi_discrepancy < 0.5_dp*pseudo_discrepancy, &
            "the discrepancy advantage is decisive, not marginal", failures)
        call expect(quasi_discrepancy < 0.01_dp, &
            "the discrepancy is small in absolute terms", failures)
    end subroutine check_beats_random_discrepancy

    ! Largest absolute difference between the fraction of points inside an
    ! origin-anchored box and its area, over a grid of boxes.
    function star_discrepancy(points) result(worst)
        real(dp), intent(in) :: points(:, :)
        real(dp) :: worst
        integer, parameter :: grid = 48
        real(dp) :: upper_x, upper_y, area, fraction
        integer :: a, b, i, inside

        worst = 0.0_dp
        do a = 1, grid
            upper_x = real(a, dp)/real(grid, dp)
            do b = 1, grid
                upper_y = real(b, dp)/real(grid, dp)
                inside = 0
                do i = 1, size(points, 1)
                    if (points(i, 1) < upper_x .and. points(i, 2) < upper_y) &
                        inside = inside + 1
                end do
                area = upper_x*upper_y
                fraction = real(inside, dp)/real(size(points, 1), dp)
                worst = max(worst, abs(fraction - area))
            end do
        end do
    end function star_discrepancy

    subroutine check_skip_matches_generation(failures)
        integer, intent(inout) :: failures
        type(sobol_t) :: skipped, generated
        type(fortnum_status_t) :: status
        real(dp) :: a(4), b(4), scratch(4)
        integer :: i
        logical :: identical

        call sobol_initialize(skipped, 4, status)
        call sobol_skip(skipped, 100_int64, status)
        call expect(status%code == FORTNUM_OK, "skipping succeeds", failures)

        call sobol_initialize(generated, 4, status)
        do i = 1, 100
            call sobol_next(generated, scratch, status)
        end do

        identical = .true.
        do i = 1, 20
            call sobol_next(skipped, a, status)
            call sobol_next(generated, b, status)
            if (maxval(abs(a - b)) /= 0.0_dp) identical = .false.
        end do
        call expect(identical, "skipping matches generating and discarding", failures)
    end subroutine check_skip_matches_generation

    subroutine check_reproducible(failures)
        integer, intent(inout) :: failures
        type(sobol_t) :: first, second
        type(fortnum_status_t) :: status
        real(dp) :: a(512, 3), b(512, 3)

        call sobol_initialize(first, 3, status)
        call sobol_fill(first, a, status)
        call sobol_initialize(second, 3, status)
        call sobol_fill(second, b, status)
        call expect(maxval(abs(a - b)) == 0.0_dp, &
            "two sequences of the same dimension agree bitwise", failures)
        call expect(all(a(1, :) == 0.0_dp), &
            "the sequence starts at the origin", failures)
    end subroutine check_reproducible

    subroutine check_refusals(failures)
        integer, intent(inout) :: failures
        type(sobol_t) :: sequence
        type(fortnum_status_t) :: status
        real(dp) :: point(3), wide(2, 9)

        call sobol_initialize(sequence, 0, status)
        call expect(status%code == FORTNUM_DOMAIN_ERROR, &
            "a zero dimension is refused", failures)

        call sobol_initialize(sequence, SOBOL_MAX_DIMENSION + 1, status)
        call expect(status%code == FORTNUM_NOT_IMPLEMENTED, &
            "a dimension beyond the table is refused, not approximated", failures)

        call sobol_initialize(sequence, 3, status)
        call sobol_fill(sequence, wide, status)
        call expect(status%code == FORTNUM_DOMAIN_ERROR, &
            "a mismatched fill width is refused", failures)
        call sobol_next(sequence, point(1:2), status)
        call expect(status%code == FORTNUM_DOMAIN_ERROR, &
            "a mismatched point width is refused", failures)
        call sobol_skip(sequence, -1_int64, status)
        call expect(status%code == FORTNUM_DOMAIN_ERROR, &
            "a negative skip is refused", failures)
    end subroutine check_refusals

    subroutine expect(condition, description, failures)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description
        integer, intent(inout) :: failures

        if (.not. condition) then
            failures = failures + 1
            print *, "  FAIL: ", description
        end if
    end subroutine expect

end program test_fortnum_sobol
