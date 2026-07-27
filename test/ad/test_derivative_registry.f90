program test_derivative_registry
    use, intrinsic :: iso_fortran_env, only: error_unit, int64, real64
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_ad_interfaces, only: fortnum_ad_status_t, ad_status_set, &
        ad_status_merge, FORTNUM_AD_BACKEND_AUTODIFF, &
        FORTNUM_AD_BACKEND_ANALYTICAL, FORTNUM_AD_BACKEND_HYBRID, &
        FORTNUM_AD_QUALITY_EXACT, FORTNUM_AD_QUALITY_APPROXIMATE, &
        FORTNUM_AD_QUALITY_NONSMOOTH
    use fortnum_derivative_registry, only: derivative_candidate_t, &
        select_derivative_candidate
    use fortnum_build_selection, only: fortnum_lookup_build_selection
    implicit none

    integer, parameter :: dp = real64
    integer :: nfail

    nfail = 0
    call test_fastest_validated_wins(nfail)
    call test_timing_tie_uses_memory(nfail)
    call test_stable_id_breaks_complete_tie(nfail)
    call test_no_validated_candidate(nfail)
    call test_hybrid_status_merge(nfail)
    call test_static_build_selection(nfail)
    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        error stop 1
    end if
    print *, "PASS"

contains

    subroutine test_fastest_validated_wins(nfail)
        integer, intent(inout) :: nfail
        type(derivative_candidate_t) :: c(3)
        type(fortnum_status_t) :: status
        integer :: selected

        c = candidate([character(7) :: "slow", "invalid", "fast"], &
            [12.0_dp, 1.0_dp, 8.0_dp], &
            [logical :: .true., .false., .true.], [100_int64, 1_int64, 200_int64])
        call select_derivative_candidate(c, selected, status)
        call check(nfail, status_ok(status) .and. selected == 3, &
            "fastest validated candidate")
    end subroutine test_fastest_validated_wins

    subroutine test_timing_tie_uses_memory(nfail)
        integer, intent(inout) :: nfail
        type(derivative_candidate_t) :: c(2)
        type(fortnum_status_t) :: status
        integer :: selected

        c = candidate([character(12) :: "lower-time", "lower-memory"], &
            [100.0_dp, 102.0_dp], &
            [logical :: .true., .true.], [1000_int64, 100_int64])
        call select_derivative_candidate(c, selected, status)
        call check(nfail, selected == 2, "three-percent tie uses memory")
    end subroutine test_timing_tie_uses_memory

    subroutine test_stable_id_breaks_complete_tie(nfail)
        integer, intent(inout) :: nfail
        type(derivative_candidate_t) :: c(2)
        type(fortnum_status_t) :: status
        integer :: selected

        c = candidate([character(5) :: "zeta", "alpha"], [10.0_dp, 10.0_dp], &
            [logical :: .true., .true.], [4_int64, 4_int64])
        call select_derivative_candidate(c, selected, status)
        call check(nfail, selected == 2, "stable candidate id tie break")
    end subroutine test_stable_id_breaks_complete_tie

    subroutine test_no_validated_candidate(nfail)
        integer, intent(inout) :: nfail
        type(derivative_candidate_t) :: c(1)
        type(fortnum_status_t) :: status
        integer :: selected

        c = candidate(["bad"], [1.0_dp], [logical :: .false.], [1_int64])
        call select_derivative_candidate(c, selected, status)
        call check(nfail, selected == 0 .and. .not. status_ok(status), &
            "no validated candidate is an error")
    end subroutine test_no_validated_candidate

    subroutine test_hybrid_status_merge(nfail)
        integer, intent(inout) :: nfail
        type(fortnum_ad_status_t) :: outer, inner, merged

        call ad_status_set(outer, 0, "", FORTNUM_AD_BACKEND_AUTODIFF, &
            FORTNUM_AD_QUALITY_EXACT)
        call ad_status_set(inner, 0, "", FORTNUM_AD_BACKEND_ANALYTICAL, &
            FORTNUM_AD_QUALITY_APPROXIMATE)
        call ad_status_merge(outer, inner, merged)
        call check(nfail, merged%backend == FORTNUM_AD_BACKEND_HYBRID .and. &
            merged%quality == FORTNUM_AD_QUALITY_APPROXIMATE, &
            "mixed mechanisms merge to hybrid")

        inner%quality = FORTNUM_AD_QUALITY_NONSMOOTH
        call ad_status_merge(outer, inner, merged)
        call check(nfail, merged%quality == FORTNUM_AD_QUALITY_NONSMOOTH, &
            "worst derivative quality propagates")
    end subroutine test_hybrid_status_merge

    subroutine test_static_build_selection(nfail)
        integer, intent(inout) :: nfail
        character(64) :: selected
        logical :: found

        call fortnum_lookup_build_selection("linear_solve", "vjp", &
            "16x16 dense system, 200000 cotangents per sample, "// &
            "transposed primal LU reusable", selected, found)
        call check(nfail, found .and. trim(selected) == &
            "reuse_transposed_primal_lu", "static workload selection")

        call fortnum_lookup_build_selection("unknown", "jvp", "none", &
            selected, found)
        call check(nfail,.not. found .and. trim(selected) == "", &
            "unknown static workload")
    end subroutine test_static_build_selection

    function candidate(ids, times, valid, memory) result(c)
        character(*), intent(in) :: ids(:)
        real(dp), intent(in) :: times(:)
        logical, intent(in) :: valid(:)
        integer(int64), intent(in) :: memory(:)
        type(derivative_candidate_t) :: c(size(ids))
        integer :: i

        do i = 1, size(ids)
            c(i)%candidate_id = ids(i)
            c(i)%median_ns = times(i)
            c(i)%validated = valid(i)
            c(i)%peak_bytes = memory(i)
            c(i)%code_bytes = 4_int64
        end do
    end function candidate

    subroutine check(nfail, condition, label)
        integer, intent(inout) :: nfail
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL "//label
        end if
    end subroutine check

end program test_derivative_registry
