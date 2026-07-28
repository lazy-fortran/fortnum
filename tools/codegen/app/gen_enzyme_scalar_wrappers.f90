program gen_enzyme_scalar_wrappers
    use fortsym_string, only: str, chars
    use fortsym_enzyme, only: enzyme_scalar_wrapper_spec_t, &
        emit_enzyme_scalar_wrapper
    use fortnum_codegen_provenance, only: fortsym_revision
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
        print "(a)", "wrote "//path
    end subroutine write_wrapper

    subroutine write_named_wrapper(directory, module_name, wrapper_prefix, &
            primal_symbol, filename)
        character(*), intent(in) :: directory, module_name, wrapper_prefix
        character(*), intent(in) :: primal_symbol, filename
        type(enzyme_scalar_wrapper_spec_t) :: spec
        character(:), allocatable :: code, path
        integer :: unit, ios

        spec%module_name = str(module_name)
        spec%wrapper_prefix = str(wrapper_prefix)
        spec%primal_symbol = str(primal_symbol)
        spec%active_inputs = 1
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
        print "(a)", "wrote "//path
    end subroutine write_named_wrapper

end program gen_enzyme_scalar_wrappers
