program gen_dawson_outer
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(**), sin
    use fortsym_products, only: jvp, vjp
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, count_operations, &
        operation_count_t, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortnum_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: x, f, v, u, value
    type(expr_t) :: values(1), variables(1), tangents(1), cotangents(1)
    type(expr_t) :: jvp_products(1), vjp_products(1)
    type(expr_t) :: value_root(1), fused_jvp(2), fused_vjp(2)
    character(len=1), parameter :: value_args(1) = ["f"]
    character(len=1), parameter :: jvp_args(3) = ["x", "f", "v"]
    character(len=1), parameter :: vjp_args(3) = ["x", "f", "u"]
    character(len=5), parameter :: value_output(1) = ["value"]
    character(len=5), parameter :: fused_jvp_outputs(2) = ["value", "jvp  "]
    character(len=3), parameter :: jvp_output(1) = ["jvp"]
    character(len=5), parameter :: fused_vjp_outputs(2) = ["value", "vjp  "]
    character(len=3), parameter :: vjp_output(1) = ["vjp"]

    call arena%init()
    engine = make_symengine_engine(arena)
    x = sym(arena, "x")
    f = sym(arena, "f")
    v = sym(arena, "v")
    u = sym(arena, "u")

    ! f is the already-computed Dawson value. This preserves the numerical
    ! approximation as an operator boundary. Every emitted candidate below is
    ! derived from this one outer-expression DAG.
    value = sin(f) + f**2
    values(1) = value
    variables(1) = f
    tangents(1) = (1 - 2*x*f)*v
    cotangents(1) = u
    jvp_products = jvp(values, variables, tangents)
    vjp_products = vjp(values, variables, cotangents)
    vjp_products(1) = vjp_products(1)*(1 - 2*x*f)

    value_root(1) = value
    fused_jvp(1) = value
    fused_jvp(2) = jvp_products(1)
    fused_vjp(1) = value
    fused_vjp(2) = vjp_products(1)
    call write_variant( &
        "fortnum_dawson_outer_value_kernel.f90", &
        "fortnum_dawson_outer_value_kernel", &
        "fortnum_generated_dawson_outer_value", &
        value_args, value_output, value_root)
    call write_variant( &
        "fortnum_dawson_outer_kernel.f90", &
        "fortnum_dawson_outer_kernel", &
        "fortnum_generated_dawson_outer", &
        jvp_args, fused_jvp_outputs, fused_jvp)
    call write_variant( &
        "fortnum_dawson_outer_jvp_kernel.f90", &
        "fortnum_dawson_outer_jvp_kernel", &
        "fortnum_generated_dawson_outer_jvp", &
        jvp_args, jvp_output, jvp_products)
    call write_variant( &
        "fortnum_dawson_outer_value_vjp_kernel.f90", &
        "fortnum_dawson_outer_value_vjp_kernel", &
        "fortnum_generated_dawson_outer_value_vjp", &
        vjp_args, fused_vjp_outputs, fused_vjp)
    call write_variant( &
        "fortnum_dawson_outer_vjp_kernel.f90", &
        "fortnum_dawson_outer_vjp_kernel", &
        "fortnum_generated_dawson_outer_vjp", &
        vjp_args, vjp_output, vjp_products)

contains

    subroutine write_variant(filename, name, module_name, arguments, outputs, &
            expressions)
        character(*), intent(in) :: filename, name, module_name
        character(*), intent(in) :: arguments(:), outputs(:)
        type(expr_t), intent(in) :: expressions(:)
        character(:), allocatable :: output, code
        type(expr_t) :: roots(size(expressions)), candidate(size(expressions))
        type(kernel_spec_t) :: spec
        type(operation_count_t) :: operations, candidate_operations
        type(engine_result_t) :: simplified
        integer :: unit, ios, k

        output = generated_path(filename)
        roots = expressions
        operations = count_operations(roots)
        candidate = roots
        do k = 1, size(candidate)
            simplified = engine%simplify(candidate(k))
            if (simplified%ok) candidate(k) = simplified%value
        end do
        candidate_operations = count_operations(candidate)
        if (candidate_operations%total < operations%total) then
            roots = candidate
            operations = candidate_operations
        end if

        spec%name = str(name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str("t")
        spec%generator = str("gen_dawson_outer")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && fo exec gen_dawson_outer")
        spec%pure_procedure = .true.
        allocate (spec%args(size(arguments)), spec%outputs(size(outputs)))
        do k = 1, size(arguments)
            spec%args(k) = str(trim(arguments(k)))
        end do
        do k = 1, size(outputs)
            spec%outputs(k) = str(trim(outputs(k)))
        end do

        open (newunit=unit, file=output, status="replace", action="write", &
            iostat=ios)
        if (ios /= 0) error stop "cannot write "//output
        code = chars(emit_kernel(roots, spec))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)

        write (output_unit, "(a)") "wrote "//output
        write (output_unit, "(a,i0)") "post-CSE structural operations: ", &
            operations%total
    end subroutine write_variant

end program gen_dawson_outer
