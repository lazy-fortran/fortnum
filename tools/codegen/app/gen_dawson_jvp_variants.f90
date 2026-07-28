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
    use fortnum_codegen_provenance, only: fortsym_revision
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: x, d, tx, td, value
    type(expr_t) :: variables(2), tangents(2), values(1)
    type(expr_t) :: raw, simplified, factored
    type(engine_result_t) :: result
    character(:), allocatable :: output_directory

    call read_output_directory(output_directory)
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
    call write_variant(output_directory, "raw", raw)
    call write_variant(output_directory, "simplified", simplified)
    call write_variant(output_directory, "factored", factored)

contains

    subroutine read_output_directory(directory)
        character(:), allocatable, intent(out) :: directory
        character(1024) :: buffer
        integer :: length, status

        call get_environment_variable( &
            "FORTNUM_VARIANT_OUTPUT_DIR", buffer, length, status)
        if (status /= 0 .or. length < 1) then
            error stop "FORTNUM_VARIANT_OUTPUT_DIR is required"
        end if
        directory = buffer(:length)
    end subroutine read_output_directory

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

    subroutine write_variant(directory, label, expression)
        character(*), intent(in) :: directory, label
        type(expr_t), intent(in) :: expression
        type(kernel_spec_t) :: spec
        type(operation_count_t) :: operations
        character(:), allocatable :: path, code
        integer :: unit, ios

        spec%name = str("fortnum_dawson_jvp_"//label)
        spec%module_name = str("fortnum_dawson_variant_"//label)
        spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str("t")
        spec%generator = str("gen_dawson_jvp_variants")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "FORTNUM_VARIANT_OUTPUT_DIR=<dir> fo exec gen_dawson_jvp_variants")
        spec%pure_procedure = .true.
        spec%openmp_declare_target = .true.
        spec%openacc_routine_seq = .true.
        spec%args = [str("x"), str("d"), str("tx"), str("td")]
        spec%outputs = [str("jvp")]

        path = directory//"/fortnum_dawson_jvp_"//label//".f90"
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
