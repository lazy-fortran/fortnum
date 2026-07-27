program gen_implicit_root_residual
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(-)
    use fortsym_products, only: jvp, implicit_tangent_rhs
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, count_operations, &
        operation_count_t, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, &
        make_symengine_engine
    use fortnum_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(expr_t) :: x(1), p(1), dp(1), residual(1), residual_x(1)
    type(expr_t) :: tangent_rhs(1), roots(3)
    type(symengine_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(kernel_spec_t) :: spec
    type(operation_count_t) :: operations
    character(:), allocatable :: output, code
    integer :: unit, ios, i

    call arena%init()
    engine = make_symengine_engine(arena)
    x(1) = sym(arena, "x")
    p(1) = sym(arena, "p")
    dp(1) = sym(arena, "tp")
    residual(1) = x(1)**2 - p(1)
    residual_x = jvp(residual, x, [num(arena, 1)])
    tangent_rhs = implicit_tangent_rhs(residual, p, dp)
    roots = [residual(1), residual_x(1), -tangent_rhs(1)]
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
    spec%args = [str("x"), str("p"), str("tp")]
    spec%outputs = [str("residual"), str("f_x"), str("f_p_tp")]

    output = generated_path("fortnum_implicit_root_residual_kernel.f90")
    open (newunit=unit, file=output, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write "//output
    code = chars(emit_kernel(roots, spec))
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)

    operations = count_operations(roots)
    write (output_unit, "(a)") "wrote "//output
    write (output_unit, "(a,i0)") "post-CSE structural operations: ", &
        operations%total
end program gen_implicit_root_residual
