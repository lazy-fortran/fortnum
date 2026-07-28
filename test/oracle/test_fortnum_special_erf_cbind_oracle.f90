program test_fortnum_special_erf_cbind_oracle
    ! Load the scipy.special.erf / scipy.special.erfc reference tables and
    ! verify fortnum_erf(x) / fortnum_erfc(x) agree within tolerance over the
    ! full real line (KAMEL/KIM domain: any sign, dense near 0 out to the
    ! saturated/underflow tails).
    !
    ! The wrappers forward to the F2008 intrinsics, so the only error is the
    ! intrinsic vs SciPy last-ULP scatter. atol = 1e-14, rtol = 1e-13.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_oracle, only: oracle_table_t, oracle_read, oracle_check
    use fortnum_special_erf_cbind, only: fortnum_erf, fortnum_erfc, &
        fortnum_erf_jvp, fortnum_erf_grad, fortnum_erfc_jvp, fortnum_erfc_grad
    implicit none

    type(oracle_table_t)   :: table_erf, table_erfc
    type(fortnum_status_t) :: status
    character(len=4096)    :: path_erf, path_erfc
    integer                :: len1, len2, stat1, stat2

    abstract interface
        subroutine product_fn(x, seed, product)
            import :: dp
            real(dp), intent(in) :: x(:), seed(:)
            real(dp), intent(out) :: product(:)
        end subroutine product_fn
    end interface

    call get_command_argument(1, path_erf, len1, stat1)
    call get_command_argument(2, path_erfc, len2, stat2)
    if (stat1 /= 0 .or. len1 == 0 .or. stat2 /= 0 .or. len2 == 0) then
        write (error_unit, "(a)") &
            "usage: test_fortnum_special_erf_cbind_oracle <erf.csv> <erfc.csv>"
        stop 1
    end if

    call oracle_read(path_erf(1:len1), table_erf, status)
    if (.not. status_ok(status)) then
        write (error_unit, "(a)") "oracle_read erf failed: "//trim(status%msg)
        stop 1
    end if
    call oracle_read(path_erfc(1:len2), table_erfc, status)
    if (.not. status_ok(status)) then
        write (error_unit, "(a)") "oracle_read erfc failed: "//trim(status%msg)
        stop 1
    end if

    if (size(table_erf%x) < 1 .or. size(table_erfc%x) < 1) then
        write (error_unit, "(a)") "empty oracle table"
        stop 1
    end if

    call oracle_check(table_erf, erf_wrapper, 1.0e-14_dp, 1.0e-13_dp, status)
    if (.not. status_ok(status)) then
        write (error_unit, "(a)") "erf failed oracle: "//trim(status%msg)
        stop 1
    end if

    call oracle_check(table_erfc, erfc_wrapper, 1.0e-14_dp, 1.0e-13_dp, status)
    if (.not. status_ok(status)) then
        write (error_unit, "(a)") "erfc failed oracle: "//trim(status%msg)
        stop 1
    end if
    call check_products(table_erf, fortnum_erf_jvp, fortnum_erf_grad, "erf")
    call check_products(table_erfc, fortnum_erfc_jvp, fortnum_erfc_grad, "erfc")

    write (*, "(a,i0,a,i0,a)") "erf/erfc oracle passed: ", &
        size(table_erf%x), " + ", size(table_erfc%x), " rows"
    stop 0

contains

    subroutine check_products(table, jvp, vjp, label)
        type(oracle_table_t), intent(in) :: table
        procedure(product_fn) :: jvp, vjp
        character(*), intent(in) :: label
        real(dp), parameter :: tangent = 0.375_dp
        real(dp), parameter :: cotangent = -1.25_dp
        real(dp) :: x(1), seed(1), got(1), expected, tolerance
        integer :: i

        if (.not. table%has_derivative) then
            write (error_unit, "(a)") trim(label)//" derivative oracle absent"
            stop 1
        end if
        do i = 1, size(table%x)
            x(1) = table%x(i)
            seed(1) = tangent
            call jvp(x, seed, got)
            expected = table%derivative(i)*tangent
            tolerance = 2.0e-15_dp + 3.0e-14_dp*abs(expected)
            if (.not. (abs(got(1) - expected) <= tolerance)) then
                write (error_unit, "(a,i0)") &
                    trim(label)//" JVP oracle failed at row ", i
                stop 1
            end if

            seed(1) = cotangent
            call vjp(x, seed, got)
            expected = table%derivative(i)*cotangent
            tolerance = 2.0e-15_dp + 3.0e-14_dp*abs(expected)
            if (.not. (abs(got(1) - expected) <= tolerance)) then
                write (error_unit, "(a,i0)") &
                    trim(label)//" VJP oracle failed at row ", i
                stop 1
            end if
        end do
    end subroutine check_products

    pure function erf_wrapper(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = fortnum_erf(x)
    end function erf_wrapper

    pure function erfc_wrapper(x) result(y)
        real(dp), intent(in) :: x
        real(dp)             :: y
        y = fortnum_erfc(x)
    end function erfc_wrapper

end program test_fortnum_special_erf_cbind_oracle
