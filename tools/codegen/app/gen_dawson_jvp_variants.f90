program gen_dawson_jvp_variants
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, exp, operator(+), operator(-), &
        operator(*)
    use fortsym_products, only: directional_derivative
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, count_operations, &
        operation_count_t, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortnum_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: x, d, tx, td, value
    type(expr_t) :: variables(2), tangents(2), values(1)
    type(expr_t) :: raw, simplified, factored
    type(engine_result_t) :: result
    character(:), allocatable :: output_directory
    logical :: emit_temporary_variants

    call read_output_directory(output_directory, emit_temporary_variants)
    call arena%init()
    engine = make_symengine_engine(arena)
    x = sym(arena, "x")
    d = sym(arena, "d")
    tx = sym(arena, "tx")
    td = sym(arena, "td")
    value = exp(-d*d) + x*d
    variables = [x, d]
    tangents = [tx, td]
    values(1) = value
    raw = directional_derivative(value, variables, tangents)

    result = engine%simplify(raw)
    if (.not. result%ok) error stop "SymEngine did not simplify raw variant"
    simplified = result%value
    factored = d*tx + td*(x - num(arena, 2)*d*exp(-d*d))

    call prove_equivalent(engine, raw, simplified, "simplified")
    call prove_equivalent(engine, raw, factored, "factored")
    call write_variant( &
        generated_path("fortnum_dawson_identity_jvp_kernel.f90"), &
        "selected simplified", "fortnum_generated_dawson_identity_jvp", &
        "fortnum_dawson_identity_jvp_kernel", &
        "cd tools/codegen && fo exec gen_dawson_jvp_variants", simplified)
    if (emit_temporary_variants) then
        call write_variant( &
            output_directory//"/fortnum_dawson_jvp_raw.f90", "raw", &
            "fortnum_dawson_variant_raw", "fortnum_dawson_jvp_raw", &
            temporary_regenerate_command(), raw)
        call write_variant( &
            output_directory//"/fortnum_dawson_jvp_factored.f90", "factored", &
            "fortnum_dawson_variant_factored", "fortnum_dawson_jvp_factored", &
            temporary_regenerate_command(), factored)
    end if

contains

    subroutine read_output_directory(directory, present)
        character(:), allocatable, intent(out) :: directory
        logical, intent(out) :: present
        character(1024) :: buffer
        integer :: length, status

        call get_environment_variable( &
            "FORTNUM_VARIANT_OUTPUT_DIR", buffer, length, status)
        present = status == 0 .and. length > 0
        if (present) then
            directory = buffer(:length)
        else
            directory = ""
        end if
    end subroutine read_output_directory

    function temporary_regenerate_command() result(command)
        character(:), allocatable :: command

        command = "cd tools/codegen && FORTNUM_VARIANT_OUTPUT_DIR=<dir> "// &
            "fo exec gen_dawson_jvp_variants"
    end function temporary_regenerate_command

    subroutine prove_equivalent(symbolic_engine, reference, candidate, label)
        type(symengine_engine_t), intent(inout) :: symbolic_engine
        type(expr_t), intent(in) :: reference, candidate
        character(*), intent(in) :: label
        type(engine_result_t) :: proof

        proof = symbolic_engine%zero_test(candidate - reference)
        if (proof%verdict /= VERDICT_TRUE) then
            error stop label//" variant was not proved equivalent"
        end if
        write (output_unit, "(a)") "proved equivalent: "//label
    end subroutine prove_equivalent

    subroutine write_variant(path, label, module_name, procedure_name, &
            regenerate_command, expression)
        character(*), intent(in) :: path, label, module_name, procedure_name
        character(*), intent(in) :: regenerate_command
        type(expr_t), intent(in) :: expression
        type(kernel_spec_t) :: spec
        type(operation_count_t) :: operations
        character(:), allocatable :: code
        integer :: unit, ios

        spec%name = str(procedure_name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str("t")
        spec%generator = str("gen_dawson_jvp_variants")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str(regenerate_command)
        spec%pure_procedure = .true.
        spec%openmp_declare_target = .true.
        spec%openacc_routine_seq = .true.
        spec%args = [str("x"), str("d"), str("tx"), str("td")]
        spec%outputs = [str("jvp")]

        open (newunit=unit, file=path, status="replace", action="write", iostat=ios)
        if (ios /= 0) error stop "cannot write "//path
        code = chars(emit_kernel([expression], spec))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)

        operations = count_operations([expression])
        write (output_unit, "(a)") "wrote "//path
        write (output_unit, "(a,i0)") label//" post-CSE operations: ", &
            operations%total
    end subroutine write_variant

end program gen_dawson_jvp_variants
