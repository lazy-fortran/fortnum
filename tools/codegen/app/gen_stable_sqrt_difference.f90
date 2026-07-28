program gen_stable_sqrt_difference
    use fortsym_string, only: str_t, str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, sqrt, operator(+), operator(-)
    use fortsym_stability, only: rationalize_sqrt_difference
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortnum_codegen_provenance, only: codegen_log, fortsym_revision, &
        generated_path
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(engine_result_t) :: simplified
    type(expr_t) :: x, a, b
    type(expr_t) :: raw(1), stable(1)

    call arena%init()
    engine = make_symengine_engine(arena)
    x = sym(arena, "x")
    a = num(arena, 1) + x
    b = num(arena, 1)
    raw(1) = sqrt(a) - sqrt(b)
    stable(1) = rationalize_sqrt_difference(a, b)
    simplified = engine%simplify(raw(1))
    if (.not. simplified%ok) error stop "cannot simplify raw sqrt difference"
    raw(1) = simplified%value
    simplified = engine%simplify(stable(1))
    if (.not. simplified%ok) error stop "cannot simplify stable sqrt difference"
    stable(1) = simplified%value
    call write_candidate("raw", raw)
    call write_candidate("stable", stable)

contains

    subroutine write_candidate(stem, expression)
        character(*), intent(in) :: stem
        type(expr_t), intent(in) :: expression(:)
        type(kernel_spec_t) :: spec
        type(str_t) :: arguments(1), outputs(1)
        character(:), allocatable :: filename, code
        integer :: unit, ios

        filename = "fortnum_sqrt1pm1_"//stem//"_kernel.f90"
        spec%name = str("fortnum_sqrt1pm1_"//stem//"_kernel")
        spec%module_name = str("fortnum_generated_sqrt1pm1_"//stem)
        spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str("t")
        spec%generator = str("gen_stable_sqrt_difference")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && fo exec gen_stable_sqrt_difference")
        spec%pure_procedure = .true.
        spec%elemental_procedure = .true.
        spec%openmp_declare_target = .true.
        spec%openacc_routine_seq = .true.
        arguments(1) = str("x")
        outputs(1) = str("value")
        spec%args = arguments
        spec%outputs = outputs

        open (newunit=unit, file=generated_path(filename), status="replace", &
            action="write", iostat=ios)
        if (ios /= 0) error stop "cannot write "//filename
        code = chars(emit_kernel(expression, spec))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)
        call codegen_log("wrote "//generated_path(filename))
    end subroutine write_candidate

end program gen_stable_sqrt_difference
