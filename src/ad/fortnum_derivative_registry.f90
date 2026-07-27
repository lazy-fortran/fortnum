module fortnum_derivative_registry
    use, intrinsic :: iso_fortran_env, only: int64, real64
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR
    use fortnum_ad_interfaces, only: FORTNUM_AD_BACKEND_AUTODIFF, &
        FORTNUM_AD_BACKEND_ANALYTICAL, FORTNUM_AD_BACKEND_HYBRID
    implicit none
    private

    integer, parameter :: dp = real64
    integer, parameter, public :: FORTNUM_PRODUCT_JVP = 1
    integer, parameter, public :: FORTNUM_PRODUCT_VJP = 2
    integer, parameter, public :: FORTNUM_PRODUCT_GRADIENT = 3
    integer, parameter, public :: FORTNUM_PRODUCT_HVP = 4
    real(dp), parameter, public :: FORTNUM_TIMING_TIE_FRACTION = 0.03_dp

    type, public :: derivative_candidate_t
        character(48) :: operator = ""
        character(48) :: workload = ""
        character(64) :: candidate_id = ""
        integer :: product = 0
        integer :: mechanism = 0
        logical :: validated = .false.
        real(dp) :: validation_error = huge(0.0_dp)
        real(dp) :: median_ns = huge(0.0_dp)
        real(dp) :: runtime_mad_ns = huge(0.0_dp)
        integer(int64) :: peak_bytes = huge(0_int64)
        integer(int64) :: code_bytes = huge(0_int64)
    end type derivative_candidate_t

    public :: select_derivative_candidate
    public :: mechanism_name

contains

    pure subroutine select_derivative_candidate(candidates, selected, status)
        type(derivative_candidate_t), intent(in) :: candidates(:)
        integer, intent(out) :: selected
        type(fortnum_status_t), intent(out) :: status
        integer :: i

        selected = 0
        call status_set(status, FORTNUM_DOMAIN_ERROR, &
            "select_derivative_candidate: no validated candidate")
        do i = 1, size(candidates)
            if (.not. usable(candidates(i))) cycle
            if (selected == 0) then
                selected = i
            else if (prefer(candidates(i), candidates(selected))) then
                selected = i
            end if
        end do
        if (selected > 0) call status_set(status, FORTNUM_OK, "")
    end subroutine select_derivative_candidate

    pure logical function usable(candidate)
        type(derivative_candidate_t), intent(in) :: candidate
        usable = candidate%validated .and. candidate%median_ns > 0.0_dp .and. &
            candidate%median_ns < huge(0.0_dp)
    end function usable

    pure logical function prefer(candidate, incumbent)
        type(derivative_candidate_t), intent(in) :: candidate, incumbent
        real(dp) :: scale

        scale = max(candidate%median_ns, incumbent%median_ns)
        if (abs(candidate%median_ns - incumbent%median_ns) > &
            FORTNUM_TIMING_TIE_FRACTION*scale) then
            prefer = candidate%median_ns < incumbent%median_ns
        else if (candidate%peak_bytes /= incumbent%peak_bytes) then
            prefer = candidate%peak_bytes < incumbent%peak_bytes
        else if (candidate%code_bytes /= incumbent%code_bytes) then
            prefer = candidate%code_bytes < incumbent%code_bytes
        else
            prefer = candidate%candidate_id < incumbent%candidate_id
        end if
    end function prefer

    pure function mechanism_name(mechanism) result(name)
        integer, intent(in) :: mechanism
        character(10) :: name
        select case (mechanism)
        case (FORTNUM_AD_BACKEND_AUTODIFF)
            name = "autodiff"
        case (FORTNUM_AD_BACKEND_ANALYTICAL)
            name = "analytical"
        case (FORTNUM_AD_BACKEND_HYBRID)
            name = "hybrid"
        case default
            name = "unknown"
        end select
    end function mechanism_name

end module fortnum_derivative_registry
