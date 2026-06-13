program test_fortnum_status
    use fortnum_status, only: fortnum_status_t, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR, &
        FORTNUM_CONVERGENCE_ERROR, FORTNUM_NOT_IMPLEMENTED, &
        status_ok, status_set
    implicit none

    type(fortnum_status_t) :: s

    ! Default-constructed status must be OK.
    if (.not. status_ok(s))          stop 1
    if (s%code /= FORTNUM_OK)        stop 1

    ! status_set with a domain error.
    call status_set(s, FORTNUM_DOMAIN_ERROR, "x < 0 not allowed")
    if (status_ok(s))                stop 1
    if (s%code /= FORTNUM_DOMAIN_ERROR) stop 1

    ! status_set with convergence error.
    call status_set(s, FORTNUM_CONVERGENCE_ERROR, "did not converge")
    if (s%code /= FORTNUM_CONVERGENCE_ERROR) stop 1

    ! status_set with not-implemented.
    call status_set(s, FORTNUM_NOT_IMPLEMENTED, "feature pending")
    if (s%code /= FORTNUM_NOT_IMPLEMENTED) stop 1

    ! Reset to OK.
    call status_set(s, FORTNUM_OK, "")
    if (.not. status_ok(s))          stop 1

    ! Named constants are distinct and non-negative.
    if (FORTNUM_OK < 0)                  stop 1
    if (FORTNUM_DOMAIN_ERROR == FORTNUM_OK) stop 1
    if (FORTNUM_CONVERGENCE_ERROR == FORTNUM_DOMAIN_ERROR) stop 1
    if (FORTNUM_NOT_IMPLEMENTED == FORTNUM_CONVERGENCE_ERROR) stop 1

    stop 0
end program test_fortnum_status
