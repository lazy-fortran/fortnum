module fortnum_build_selection
    ! Fallback for build systems that do not consume CMake benchmark records.
    ! CMake compiles its generated module of the same name instead of this file.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_generated_dawson_outer, only: fortnum_dawson_outer_kernel
    implicit none
    private

    character(*), parameter, public :: FORTNUM_DAWSON_OUTER_JVP_CANDIDATE = &
        "analytical"

    public :: fortnum_selected_dawson_outer_jvp
    public :: fortnum_lookup_build_selection

    integer, parameter :: registry_size = 4
    character(48), parameter :: registry_operators(registry_size) = [character(48) :: &
        "dawson_outer", "linear_solve", "linear_solve", "multiroot_implicit"]
    character(16), parameter :: registry_products(registry_size) = [character(16) :: &
        "jvp", "jvp", "vjp", "jvp"]
    character(128), parameter :: registry_workloads(registry_size) = [character(128) :: &
        "scalar calls, x cycled over [0.65, 0.75], one direction", &
        "16x16 dense system, 200000 directions per sample, primal LU reusable", &
        "16x16 dense system, 200000 cotangents per sample, transposed primal LU reusable", &
        "16x16 dense residual Jacobian, 200000 parameter directions per sample"]
    character(64), parameter :: registry_candidates(registry_size) = [character(64) :: &
        "analytical", "reuse_primal_lu", "reuse_transposed_primal_lu", &
        "default_dense_solve"]

contains

    subroutine fortnum_selected_dawson_outer_jvp(x, f, v, value, jvp)
        real(dp), intent(in) :: x, f, v
        real(dp), intent(out) :: value, jvp

        call fortnum_dawson_outer_kernel(x, f, v, value, jvp)
    end subroutine fortnum_selected_dawson_outer_jvp

    pure subroutine fortnum_lookup_build_selection(operator, product, workload, &
            candidate, found)
        character(*), intent(in) :: operator, product, workload
        character(*), intent(out) :: candidate
        logical, intent(out) :: found
        integer :: i

        candidate = ""
        found = .false.
        do i = 1, registry_size
            if (trim(operator) /= trim(registry_operators(i))) cycle
            if (trim(product) /= trim(registry_products(i))) cycle
            if (trim(workload) /= trim(registry_workloads(i))) cycle
            candidate = trim(registry_candidates(i))
            found = .true.
            return
        end do
    end subroutine fortnum_lookup_build_selection

end module fortnum_build_selection
