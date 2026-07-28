program gen_enzyme_scalar_wrappers
    use fortsym_string, only: str, chars
    use fortsym_enzyme, only: enzyme_fixed_array_wrapper_spec_t, &
        enzyme_scalar_vector_wrapper_spec_t, enzyme_scalar_wrapper_spec_t, &
        emit_enzyme_fixed_array_wrapper, emit_enzyme_scalar_vector_wrapper, &
        emit_enzyme_scalar_wrapper
    use fortnum_codegen_provenance, only: codegen_log, fortsym_revision
    implicit none

    character(:), allocatable :: output_directory
    integer :: active_inputs

    call read_output_directory(output_directory)
    do active_inputs = 1, 4
        call write_wrapper(output_directory, active_inputs)
    end do
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_bessel_outer", &
        "fortnum_enzyme_bessel_outer", "fortnum_bessel_outer", &
        "fortnum_enzyme_bessel_outer.f90")
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_bessel_outer_autodiff", &
        "fortnum_enzyme_bessel_outer_autodiff", &
        "fortnum_bessel_outer_autodiff", &
        "fortnum_enzyme_bessel_outer_autodiff.f90")
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_gamma_outer", &
        "fortnum_enzyme_gamma_outer", "fortnum_gamma_outer", &
        "fortnum_enzyme_gamma_outer.f90")
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_gamma_outer_autodiff", &
        "fortnum_enzyme_gamma_outer_autodiff", &
        "fortnum_gamma_outer_autodiff", &
        "fortnum_enzyme_gamma_outer_autodiff.f90")
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_erf_outer", &
        "fortnum_enzyme_erf_outer", "fortnum_erf_outer", &
        "fortnum_enzyme_erf_outer.f90")
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_erf_outer_autodiff", &
        "fortnum_enzyme_erf_outer_autodiff", &
        "fortnum_erf_outer_autodiff", &
        "fortnum_enzyme_erf_outer_autodiff.f90")
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_dawson_outer", &
        "fortnum_enzyme_dawson_outer", "fortnum_dawson_outer", &
        "fortnum_enzyme_dawson_outer.f90")
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_dawson_outer_autodiff", &
        "fortnum_enzyme_dawson_outer_autodiff", &
        "fortnum_dawson_outer_autodiff", &
        "fortnum_enzyme_dawson_outer_autodiff.f90")
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_square", "fortnum_enzyme_square", &
        "fortnum_smoke_square", "fortnum_enzyme_square.f90")
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_scalar_root_residual", &
        "fortnum_enzyme_scalar_root_residual", "fortnum_scalar_root_residual", &
        "fortnum_enzyme_scalar_root_residual.f90", 3)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_scalar_root_newton", &
        "fortnum_enzyme_scalar_root_newton", "fortnum_scalar_root_newton_solve", &
        "fortnum_enzyme_scalar_root_newton.f90", 3)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_scalar_root_vjp_residual", &
        "fortnum_enzyme_scalar_root_vjp_residual", &
        "fortnum_scalar_root_vjp_residual", &
        "fortnum_enzyme_scalar_root_vjp_residual.f90", 3)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_scalar_root_vjp_newton", &
        "fortnum_enzyme_scalar_root_vjp_newton", &
        "fortnum_scalar_root_vjp_newton_solve", &
        "fortnum_enzyme_scalar_root_vjp_newton.f90", 3)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_vector_root_residual_one", &
        "fortnum_enzyme_vector_root_residual_one", &
        "fortnum_vector_root_residual_one", &
        "fortnum_enzyme_vector_root_residual_one.f90", 4)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_vector_root_residual_two", &
        "fortnum_enzyme_vector_root_residual_two", &
        "fortnum_vector_root_residual_two", &
        "fortnum_enzyme_vector_root_residual_two.f90", 4)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_vector_root_newton_one", &
        "fortnum_enzyme_vector_root_newton_one", &
        "fortnum_vector_newton_root_one", &
        "fortnum_enzyme_vector_root_newton_one.f90", 4)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_vector_root_newton_two", &
        "fortnum_enzyme_vector_root_newton_two", &
        "fortnum_vector_newton_root_two", &
        "fortnum_enzyme_vector_root_newton_two.f90", 4)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_vector_root_vjp_residual_one", &
        "fortnum_enzyme_vector_root_vjp_residual_one", &
        "fortnum_vector_root_vjp_residual_one", &
        "fortnum_enzyme_vector_root_vjp_residual_one.f90", 4)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_vector_root_vjp_residual_two", &
        "fortnum_enzyme_vector_root_vjp_residual_two", &
        "fortnum_vector_root_vjp_residual_two", &
        "fortnum_enzyme_vector_root_vjp_residual_two.f90", 4)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_vector_root_objective", &
        "fortnum_enzyme_vector_root_objective", &
        "fortnum_vector_newton_objective", &
        "fortnum_enzyme_vector_root_objective.f90", 4)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_fixed_quadrature_integrand", &
        "fortnum_enzyme_fixed_quadrature_integrand", &
        "fortnum_fixed_quadrature_integrand", &
        "fortnum_enzyme_fixed_quadrature_integrand.f90", 5)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_fixed_quadrature_kernel", &
        "fortnum_enzyme_fixed_quadrature_kernel", &
        "fortnum_fixed_quadrature_kernel", &
        "fortnum_enzyme_fixed_quadrature_kernel.f90", 4)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_adaptive_integrand", &
        "fortnum_enzyme_adaptive_integrand", &
        "fortnum_adaptive_trace_integrand", &
        "fortnum_enzyme_adaptive_integrand.f90", 2)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_singular_integrand", &
        "fortnum_enzyme_singular_integrand", &
        "fortnum_singular_trace_integrand", &
        "fortnum_enzyme_singular_integrand.f90", 2)
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_adaptive_frozen_trace", &
        "fortnum_enzyme_adaptive_frozen_trace", &
        "fortnum_adaptive_frozen_trace_value", &
        "fortnum_enzyme_adaptive_frozen_trace.f90")
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_singular_frozen_trace", &
        "fortnum_enzyme_singular_frozen_trace", &
        "fortnum_singular_frozen_trace_value", &
        "fortnum_enzyme_singular_frozen_trace.f90")
    call write_named_wrapper(output_directory, &
        "fortnum_generated_enzyme_ode_scalar_rhs", &
        "fortnum_enzyme_ode_scalar_rhs", "fortnum_ode_scalar_rhs", &
        "fortnum_enzyme_ode_scalar_rhs.f90", 3)
    call write_scalar_vector_wrapper(output_directory)
    call write_fixed_array_wrapper(output_directory, &
        "fortnum_generated_enzyme_direct_solver_component", &
        "fortnum_enzyme_direct_solver_component", &
        "fortnum_direct_solve_component", [16, 4], 1, &
        "fortnum_enzyme_direct_solver_component.f90")
    call write_fixed_array_wrapper(output_directory, &
        "fortnum_generated_enzyme_direct_solver_objective", &
        "fortnum_enzyme_direct_solver_objective", &
        "fortnum_direct_solve_objective", [16, 4, 4], 0, &
        "fortnum_enzyme_direct_solver_objective.f90")
    call write_fixed_array_wrapper(output_directory, &
        "fortnum_generated_enzyme_iterative_solver_component", &
        "fortnum_enzyme_iterative_solver_component", &
        "fortnum_iterative_solve_component", [16, 4], 2, &
        "fortnum_enzyme_iterative_solver_component.f90")

contains

    subroutine read_output_directory(directory)
        character(:), allocatable, intent(out) :: directory
        character(4096) :: buffer
        integer :: length, status

        call get_environment_variable( &
            "FORTNUM_ENZYME_WRAPPER_OUTPUT_DIR", buffer, length, status)
        if (status /= 0 .or. length < 1) then
            error stop "FORTNUM_ENZYME_WRAPPER_OUTPUT_DIR is required"
        end if
        directory = buffer(:length)
    end subroutine read_output_directory

    subroutine write_wrapper(directory, count)
        character(*), intent(in) :: directory
        integer, intent(in) :: count
        type(enzyme_scalar_wrapper_spec_t) :: spec
        character(:), allocatable :: code, path, suffix
        character(16) :: count_text
        integer :: unit, ios

        write (count_text, "(i0)") count
        suffix = "p"//trim(count_text)
        spec%module_name = str("fortnum_generated_enzyme_scalar_"//suffix)
        spec%wrapper_prefix = str("fortnum_enzyme_scalar_"//suffix)
        spec%primal_symbol = str("fortnum_scalar_primal_"//suffix)
        spec%active_inputs = count
        spec%generator = str("gen_enzyme_scalar_wrappers")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && FORTNUM_ENZYME_WRAPPER_OUTPUT_DIR=<dir> "// &
            "fo exec gen_enzyme_scalar_wrappers")
        if (count == 1) then
            spec%analytical_jvp_symbol = str( &
                "fortnum_scalar_analytical_p1_jvp")
            spec%custom_forward_symbol = str( &
                "fortnum_scalar_custom_p1_forward")
            spec%custom_forward_counter_symbol = str( &
                "fortnum_enzyme_rule_counter_record")
        end if

        path = directory//"/fortnum_enzyme_scalar_"//suffix//".f90"
        open (newunit=unit, file=path, status="replace", action="write", &
            iostat=ios)
        if (ios /= 0) error stop "cannot write "//path
        code = chars(emit_enzyme_scalar_wrapper(spec))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)
        call codegen_log("wrote "//path)
    end subroutine write_wrapper

    subroutine write_named_wrapper(directory, module_name, wrapper_prefix, &
            primal_symbol, filename, active_inputs)
        character(*), intent(in) :: directory, module_name, wrapper_prefix
        character(*), intent(in) :: primal_symbol, filename
        integer, intent(in), optional :: active_inputs
        type(enzyme_scalar_wrapper_spec_t) :: spec
        character(:), allocatable :: code, path
        integer :: unit, ios

        spec%module_name = str(module_name)
        spec%wrapper_prefix = str(wrapper_prefix)
        spec%primal_symbol = str(primal_symbol)
        spec%active_inputs = 1
        if (present(active_inputs)) spec%active_inputs = active_inputs
        spec%generator = str("gen_enzyme_scalar_wrappers")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && FORTNUM_ENZYME_WRAPPER_OUTPUT_DIR=<dir> "// &
            "fo exec gen_enzyme_scalar_wrappers")

        path = directory//"/"//filename
        open (newunit=unit, file=path, status="replace", action="write", &
            iostat=ios)
        if (ios /= 0) error stop "cannot write "//path
        code = chars(emit_enzyme_scalar_wrapper(spec))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)
        call codegen_log("wrote "//path)
    end subroutine write_named_wrapper

    subroutine write_scalar_vector_wrapper(directory)
        character(*), intent(in) :: directory
        type(enzyme_scalar_vector_wrapper_spec_t) :: spec
        character(:), allocatable :: code, path
        integer :: unit, ios

        spec%module_name = str("fortnum_generated_enzyme_bspline_fixed_span")
        spec%wrapper_prefix = str("fortnum_enzyme_bspline_fixed_span")
        spec%primal_symbol = str("fortnum_bspline_fixed_span_value")
        spec%vector_size = 4
        spec%generator = str("gen_enzyme_scalar_wrappers")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && FORTNUM_ENZYME_WRAPPER_OUTPUT_DIR=<dir> "// &
            "fo exec gen_enzyme_scalar_wrappers")

        path = directory//"/fortnum_enzyme_bspline_fixed_span.f90"
        open (newunit=unit, file=path, status="replace", action="write", &
            iostat=ios)
        if (ios /= 0) error stop "cannot write "//path
        code = chars(emit_enzyme_scalar_vector_wrapper(spec))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)
        call codegen_log("wrote "//path)
    end subroutine write_scalar_vector_wrapper

    subroutine write_fixed_array_wrapper(directory, module_name, wrapper_prefix, &
            primal_symbol, array_sizes, inactive_integer_count, filename)
        character(*), intent(in) :: directory, module_name, wrapper_prefix
        character(*), intent(in) :: primal_symbol, filename
        integer, intent(in) :: array_sizes(:)
        integer, intent(in) :: inactive_integer_count
        type(enzyme_fixed_array_wrapper_spec_t) :: spec
        character(:), allocatable :: code, path
        integer :: unit, ios

        spec%module_name = str(module_name)
        spec%wrapper_prefix = str(wrapper_prefix)
        spec%primal_symbol = str(primal_symbol)
        spec%array_sizes = array_sizes
        spec%inactive_integer_count = inactive_integer_count
        spec%generator = str("gen_enzyme_scalar_wrappers")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && FORTNUM_ENZYME_WRAPPER_OUTPUT_DIR=<dir> "// &
            "fo exec gen_enzyme_scalar_wrappers")

        path = directory//"/"//filename
        open (newunit=unit, file=path, status="replace", action="write", &
            iostat=ios)
        if (ios /= 0) error stop "cannot write "//path
        code = chars(emit_enzyme_fixed_array_wrapper(spec))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)
        call codegen_log("wrote "//path)
    end subroutine write_fixed_array_wrapper

end program gen_enzyme_scalar_wrappers
