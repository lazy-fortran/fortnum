program gen_dawson_outer
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(**), sin
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, count_operations, &
        operation_count_t, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortnum_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    character(:), allocatable :: output
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: x, f, v, value, roots(2), candidate(2)
    type(kernel_spec_t) :: spec
    type(operation_count_t) :: operations, candidate_operations
    type(engine_result_t) :: simplified
    character(:), allocatable :: code
    integer :: unit, ios, k

    call arena%init()
    engine = make_symengine_engine(arena)
    output = generated_path("fortnum_dawson_outer_kernel.f90")
    x = sym(arena, "x")
    f = sym(arena, "f")
    v = sym(arena, "v")

    ! f is the already-computed Dawson value. This preserves the numerical
    ! approximation as an operator boundary while fortsym applies the outer
    ! product and chain rules and fuses value plus JVP.
    value = sin(f) + f**2
    roots(1) = value
    roots(2) = diff(value, f)*(1 - 2*x*f)*v
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

    spec%name = str("fortnum_dawson_outer_kernel")
    spec%module_name = str("fortnum_generated_dawson_outer")
    spec%mode = KERNEL_SUBROUTINE
    spec%temp_prefix = str("t")
    spec%generator = str("gen_dawson_outer")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str( &
        "cd tools/codegen && fo exec gen_dawson_outer")
    spec%pure_procedure = .true.
    allocate (spec%args(3), spec%outputs(2))
    spec%args = [str("x"), str("f"), str("v")]
    spec%outputs = [str("value"), str("jvp")]

    open (newunit=unit, file=output, status="replace", action="write", iostat=ios)
    if (ios /= 0) error stop "cannot write "//output
    code = chars(emit_kernel(roots, spec))
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)

    write (output_unit, "(a)") "wrote "//output
    write (output_unit, "(a,i0)") "post-CSE structural operations: ", &
        operations%total
end program gen_dawson_outer
