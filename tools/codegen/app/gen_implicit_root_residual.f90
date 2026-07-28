program gen_implicit_root_residual
    use fortsym_string, only: str_t, str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(**)
    use fortsym_products, only: jvp, vjp, implicit_tangent_rhs
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, count_operations, &
        operation_count_t, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, &
        make_symengine_engine
    use fortnum_codegen_provenance, only: codegen_log, codegen_log_count, &
        fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(expr_t) :: x, p, dp, residual
    type(expr_t) :: p1, p2, tp1, tp2, u, scalar_residual
    type(expr_t) :: x1, x2, u1, u2
    type(expr_t) :: x_variables(1), p_variables(1), parameter_tangent(1)
    type(expr_t) :: residuals(1), residual_x(1)
    type(expr_t) :: tangent_rhs(1), roots(3)
    type(expr_t) :: scalar_residuals(1), scalar_x_variables(1)
    type(expr_t) :: scalar_p_variables(2), scalar_p_tangents(2)
    type(expr_t) :: scalar_cotangents(1), scalar_f_x(1)
    type(expr_t) :: scalar_f_p_tp(1), scalar_f_p_t_u(2)
    type(expr_t) :: scalar_jvp_roots(2), scalar_vjp_roots(3)
    type(expr_t) :: unit_tangent(1)
    type(expr_t) :: vector_residuals(2), vector_x_variables(2)
    type(expr_t) :: vector_p_variables(2), vector_p_tangents(2)
    type(expr_t) :: vector_state_direction_one(2)
    type(expr_t) :: vector_state_direction_two(2)
    type(expr_t) :: vector_cotangents(2)
    type(expr_t) :: vector_jacobian_column_one(2)
    type(expr_t) :: vector_jacobian_column_two(2)
    type(expr_t) :: vector_parameter_jvp(2), vector_parameter_vjp(2)
    type(expr_t) :: vector_jvp_roots(6), vector_jacobian_roots(4)
    type(expr_t) :: vector_vjp_roots(2)
    type(str_t) :: base_args(3), base_outputs(3)
    type(str_t) :: scalar_jvp_args(4), scalar_jvp_outputs(2)
    type(str_t) :: scalar_vjp_args(3), scalar_vjp_outputs(3)
    type(str_t) :: vector_jvp_args(2), vector_jvp_arg_shapes(2)
    type(str_t) :: vector_jvp_outputs(2), vector_jvp_output_shapes(2)
    type(str_t) :: vector_jvp_references(6)
    type(str_t) :: vector_jacobian_args(1), vector_jacobian_arg_shapes(1)
    type(str_t) :: vector_jacobian_outputs(1)
    type(str_t) :: vector_jacobian_output_shapes(1)
    type(str_t) :: vector_jacobian_references(4)
    type(str_t) :: vector_vjp_args(1), vector_vjp_arg_shapes(1)
    type(str_t) :: vector_vjp_outputs(1), vector_vjp_output_shapes(1)
    type(str_t) :: vector_vjp_references(2)
    type(symengine_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(kernel_spec_t) :: spec
    type(operation_count_t) :: operations
    character(:), allocatable :: output, code
    integer :: unit, ios, i

    call arena%init()
    engine = make_symengine_engine(arena)
    x = sym(arena, "x")
    p = sym(arena, "p")
    dp = sym(arena, "tp")
    residual = x*x - p
    x_variables(1) = x
    p_variables(1) = p
    parameter_tangent(1) = dp
    residuals(1) = residual
    unit_tangent(1) = num(arena, 1)
    residual_x = jvp(residuals, x_variables, unit_tangent)
    tangent_rhs = implicit_tangent_rhs( &
        residuals, p_variables, parameter_tangent)
    roots(1) = residual
    roots(2) = residual_x(1)
    roots(3) = -tangent_rhs(1)
    do i = 1, size(roots)
        simplified = engine%simplify(roots(i))
        if (simplified%ok) roots(i) = simplified%value
    end do

    spec%name = str("fortnum_implicit_root_residual_kernel")
    spec%module_name = str("fortnum_generated_implicit_root_residual")
    spec%mode = KERNEL_SUBROUTINE
    spec%temp_prefix = str("t")
    spec%generator = str("gen_implicit_root_residual")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str( &
        "cd tools/codegen && fo exec gen_implicit_root_residual")
    spec%pure_procedure = .true.
    spec%openmp_declare_target = .true.
    spec%openacc_routine_seq = .true.
    base_args(1) = str("x")
    base_args(2) = str("p")
    base_args(3) = str("tp")
    base_outputs(1) = str("residual")
    base_outputs(2) = str("f_x")
    base_outputs(3) = str("f_p_tp")
    spec%args = base_args
    spec%outputs = base_outputs

    output = generated_path("fortnum_implicit_root_residual_kernel.f90")
    open (newunit=unit, file=output, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write "//output
    code = chars(emit_kernel(roots, spec))
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)

    operations = count_operations(roots)
    call codegen_log("wrote "//output)
    call codegen_log_count("post-CSE structural operations: ", operations%total)

    p1 = sym(arena, "p1")
    p2 = sym(arena, "p2")
    tp1 = sym(arena, "tp1")
    tp2 = sym(arena, "tp2")
    u = sym(arena, "u")
    scalar_residual = x**3 + p1*x - p2
    scalar_residuals(1) = scalar_residual
    scalar_x_variables(1) = x
    scalar_p_variables(1) = p1
    scalar_p_variables(2) = p2
    scalar_p_tangents(1) = tp1
    scalar_p_tangents(2) = tp2
    scalar_cotangents(1) = u
    scalar_f_x = jvp( &
        scalar_residuals, scalar_x_variables, unit_tangent)
    scalar_f_p_tp = jvp( &
        scalar_residuals, scalar_p_variables, scalar_p_tangents)
    scalar_f_p_t_u = vjp( &
        scalar_residuals, scalar_p_variables, scalar_cotangents)
    scalar_jvp_roots(1) = scalar_f_x(1)
    scalar_jvp_roots(2) = scalar_f_p_tp(1)
    scalar_vjp_roots(1) = scalar_f_x(1)
    scalar_vjp_roots(2) = scalar_f_p_t_u(1)
    scalar_vjp_roots(3) = scalar_f_p_t_u(2)
    scalar_jvp_args(1) = str("x")
    scalar_jvp_args(2) = str("p1")
    scalar_jvp_args(3) = str("tp1")
    scalar_jvp_args(4) = str("tp2")
    scalar_jvp_outputs(1) = str("f_x")
    scalar_jvp_outputs(2) = str("f_p_tp")
    scalar_vjp_args(1) = str("x")
    scalar_vjp_args(2) = str("p1")
    scalar_vjp_args(3) = str("u")
    scalar_vjp_outputs(1) = str("f_x")
    scalar_vjp_outputs(2) = str("f_p1_t_u")
    scalar_vjp_outputs(3) = str("f_p2_t_u")
    call write_product( &
        "fortnum_scalar_root_residual_jvp_kernel.f90", &
        "fortnum_scalar_root_residual_jvp_kernel", &
        "fortnum_generated_scalar_root_residual_jvp", &
        scalar_jvp_args, scalar_jvp_outputs, scalar_jvp_outputs, &
        scalar_jvp_roots)
    call write_product( &
        "fortnum_scalar_root_residual_vjp_kernel.f90", &
        "fortnum_scalar_root_residual_vjp_kernel", &
        "fortnum_generated_scalar_root_residual_vjp", &
        scalar_vjp_args, scalar_vjp_outputs, scalar_vjp_outputs, &
        scalar_vjp_roots)

    x1 = sym(arena, "x(1)")
    x2 = sym(arena, "x(2)")
    u1 = sym(arena, "u(1)")
    u2 = sym(arena, "u(2)")
    vector_residuals(1) = x1**2 + x2 - p1
    vector_residuals(2) = x1 + x2**2 - p2
    vector_x_variables(1) = x1
    vector_x_variables(2) = x2
    vector_p_variables(1) = p1
    vector_p_variables(2) = p2
    vector_p_tangents(1) = sym(arena, "tp(1)")
    vector_p_tangents(2) = sym(arena, "tp(2)")
    vector_state_direction_one(1) = num(arena, 1)
    vector_state_direction_one(2) = num(arena, 0)
    vector_state_direction_two(1) = num(arena, 0)
    vector_state_direction_two(2) = num(arena, 1)
    vector_cotangents(1) = u1
    vector_cotangents(2) = u2
    vector_jacobian_column_one = jvp( &
        vector_residuals, vector_x_variables, vector_state_direction_one)
    vector_jacobian_column_two = jvp( &
        vector_residuals, vector_x_variables, vector_state_direction_two)
    vector_parameter_jvp = jvp( &
        vector_residuals, vector_p_variables, vector_p_tangents)
    vector_parameter_vjp = vjp( &
        vector_residuals, vector_p_variables, vector_cotangents)
    vector_jvp_roots(1) = vector_jacobian_column_one(1)
    vector_jvp_roots(2) = vector_jacobian_column_two(1)
    vector_jvp_roots(3) = vector_jacobian_column_one(2)
    vector_jvp_roots(4) = vector_jacobian_column_two(2)
    vector_jvp_roots(5) = vector_parameter_jvp(1)
    vector_jvp_roots(6) = vector_parameter_jvp(2)
    vector_jacobian_roots = vector_jvp_roots(1:4)
    vector_vjp_roots = vector_parameter_vjp
    vector_jvp_args(1) = str("x")
    vector_jvp_args(2) = str("tp")
    vector_jvp_arg_shapes(1) = str("(2)")
    vector_jvp_arg_shapes(2) = str("(2)")
    vector_jvp_outputs(1) = str("jac_x")
    vector_jvp_outputs(2) = str("f_p_tp")
    vector_jvp_output_shapes(1) = str("(2,2)")
    vector_jvp_output_shapes(2) = str("(2)")
    vector_jvp_references(1) = str("jac_x(1,1)")
    vector_jvp_references(2) = str("jac_x(1,2)")
    vector_jvp_references(3) = str("jac_x(2,1)")
    vector_jvp_references(4) = str("jac_x(2,2)")
    vector_jvp_references(5) = str("f_p_tp(1)")
    vector_jvp_references(6) = str("f_p_tp(2)")
    vector_jacobian_args(1) = str("x")
    vector_jacobian_arg_shapes(1) = str("(2)")
    vector_jacobian_outputs(1) = str("jac_x")
    vector_jacobian_output_shapes(1) = str("(2,2)")
    vector_jacobian_references = vector_jvp_references(1:4)
    vector_vjp_args(1) = str("u")
    vector_vjp_arg_shapes(1) = str("(2)")
    vector_vjp_outputs(1) = str("f_p_t_u")
    vector_vjp_output_shapes(1) = str("(2)")
    vector_vjp_references(1) = str("f_p_t_u(1)")
    vector_vjp_references(2) = str("f_p_t_u(2)")
    call write_product( &
        "fortnum_vector_root_residual_jvp_kernel.f90", &
        "fortnum_vector_root_residual_jvp_kernel", &
        "fortnum_generated_vector_root_residual_jvp", &
        vector_jvp_args, vector_jvp_outputs, vector_jvp_references, &
        vector_jvp_roots, vector_jvp_arg_shapes, vector_jvp_output_shapes)
    call write_product( &
        "fortnum_vector_root_residual_jacobian_kernel.f90", &
        "fortnum_vector_root_residual_jacobian_kernel", &
        "fortnum_generated_vector_root_residual_jacobian", &
        vector_jacobian_args, vector_jacobian_outputs, &
        vector_jacobian_references, vector_jacobian_roots, &
        vector_jacobian_arg_shapes, vector_jacobian_output_shapes)
    call write_product( &
        "fortnum_vector_root_residual_vjp_kernel.f90", &
        "fortnum_vector_root_residual_vjp_kernel", &
        "fortnum_generated_vector_root_residual_vjp", &
        vector_vjp_args, vector_vjp_outputs, vector_vjp_references, &
        vector_vjp_roots, vector_vjp_arg_shapes, vector_vjp_output_shapes)

contains

    subroutine write_product(filename, name, module_name, arguments, &
            outputs, output_references, expressions, arg_shapes, output_shapes)
        character(*), intent(in) :: filename, name, module_name
        type(str_t), intent(in) :: arguments(:), outputs(:), output_references(:)
        type(expr_t), intent(in) :: expressions(:)
        type(str_t), intent(in), optional :: arg_shapes(:)
        type(str_t), intent(in), optional :: output_shapes(:)
        type(kernel_spec_t) :: product_spec
        type(expr_t) :: product_roots(size(expressions))
        type(engine_result_t) :: result
        type(operation_count_t) :: product_operations
        character(:), allocatable :: product_code
        integer :: product_unit, product_ios, k

        product_roots = expressions
        do k = 1, size(product_roots)
            result = engine%simplify(product_roots(k))
            if (result%ok) product_roots(k) = result%value
        end do
        product_spec%name = str(name)
        product_spec%module_name = str(module_name)
        product_spec%mode = KERNEL_SUBROUTINE
        product_spec%temp_prefix = str("t")
        product_spec%generator = str("gen_implicit_root_residual")
        product_spec%generator_revision = str(fortsym_revision())
        product_spec%regenerate_command = str( &
            "cd tools/codegen && fo exec gen_implicit_root_residual")
        product_spec%pure_procedure = .true.
        product_spec%openmp_declare_target = .true.
        product_spec%openacc_routine_seq = .true.
        product_spec%args = arguments
        product_spec%outputs = outputs
        product_spec%output_references = output_references
        if (present(arg_shapes)) product_spec%arg_shapes = arg_shapes
        if (present(output_shapes)) product_spec%output_shapes = output_shapes

        open (newunit=product_unit, file=generated_path(filename), &
            status="replace", action="write", iostat=product_ios)
        if (product_ios /= 0) error stop "cannot write "//filename
        product_code = chars(emit_kernel(product_roots, product_spec))
        write (product_unit, "(a)") product_code(:len(product_code) - 1)
        close (product_unit)
        product_operations = count_operations(product_roots)
        call codegen_log("wrote "//generated_path(filename))
        call codegen_log_count( &
            "post-CSE structural operations: ", product_operations%total)
    end subroutine write_product

end program gen_implicit_root_residual
