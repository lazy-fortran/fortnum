program gen_special_outer
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(**), sin
    use fortsym_products, only: jvp, vjp
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, count_operations, &
        operation_count_t, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortnum_codegen_provenance, only: codegen_log, codegen_log_count, &
        fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: x, f, df, v, u, value
    type(expr_t) :: values(1), variables(1), tangents(1), cotangents(1)
    type(expr_t) :: jvp_products(1), vjp_products(1)
    type(expr_t) :: jvp_partial(1), vjp_partial(1)
    type(expr_t) :: bessel_tangent(1), bessel_jvp(1), bessel_partial(1)
    type(expr_t) :: value_root(1), fused_jvp(2), fused_vjp(2)
    character(len=1), parameter :: value_args(1) = ["f"]
    character(len=2), parameter :: bessel_jvp_args(3) = ["f ", "df", "v "]
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
    df = sym(arena, "df")
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
    bessel_tangent(1) = df*v
    bessel_jvp = jvp(values, variables, bessel_tangent)
    bessel_partial(1) = diff(value, f)*df*v
    bessel_jvp(1) = lower_operation_variant( &
        bessel_jvp(1), bessel_partial(1), engine)
    jvp_partial(1) = diff(value, f)*tangents(1)
    vjp_partial(1) = u*diff(value, f)*(1 - 2*x*f)
    jvp_products(1) = lower_operation_variant( &
        jvp_products(1), jvp_partial(1), engine)
    vjp_products(1) = lower_operation_variant( &
        vjp_products(1), vjp_partial(1), engine)

    value_root(1) = value
    fused_jvp(1) = value
    fused_jvp(2) = jvp_products(1)
    fused_vjp(1) = value
    fused_vjp(2) = vjp_products(1)
    call write_variant( &
        "fortnum_bessel_outer_jvp_kernel.f90", &
        "fortnum_bessel_outer_jvp_kernel", &
        "fortnum_generated_bessel_outer_jvp", &
        bessel_jvp_args, jvp_output, bessel_jvp)
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

    function lower_operation_variant(first, second, engine) result(best)
        type(expr_t), intent(in) :: first, second
        type(symengine_engine_t), intent(inout) :: engine
        type(expr_t) :: best
        type(expr_t) :: candidate(1)
        type(operation_count_t) :: best_operations, candidate_operations
        type(engine_result_t) :: simplified

        best = first
        candidate(1) = best
        best_operations = count_operations(candidate)
        candidate(1) = second
        simplified = engine%simplify(candidate(1))
        if (simplified%ok) candidate(1) = simplified%value
        candidate_operations = count_operations(candidate)
        if (candidate_operations%total < best_operations%total) then
            best = candidate(1)
        end if
    end function lower_operation_variant

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
        spec%generator = str("gen_special_outer")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && fo exec gen_special_outer")
        spec%pure_procedure = .true.
        spec%openmp_declare_target = .true.
        spec%openacc_routine_seq = .true.
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

        call codegen_log("wrote "//output)
        call codegen_log_count("post-CSE structural operations: ", &
            operations%total)
    end subroutine write_variant

end program gen_special_outer
