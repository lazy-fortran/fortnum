program gen_jacobi_recurrence
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, operator(*), operator(+), &
        operator(-), operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    use fortnum_codegen_provenance, only: codegen_log, fortsym_revision, &
        generated_path
    implicit none

    type(arena_t), target :: arena
    type(expr_t) :: alpha, beta, current, degree, denominator
    type(expr_t) :: numerator, previous, roots(1), x
    type(kernel_spec_t) :: spec
    type(str_t) :: arguments(6), outputs(1)
    character(:), allocatable :: code, filename
    integer :: ios, unit

    call arena%init()
    degree = sym(arena, "degree")
    alpha = sym(arena, "alpha")
    beta = sym(arena, "beta")
    x = sym(arena, "x")
    previous = sym(arena, "previous")
    current = sym(arena, "current")

    denominator = 2*degree*(degree + alpha + beta)* &
        (2*degree + alpha + beta - 2)
    numerator = (2*degree + alpha + beta - 1)* &
        ((2*degree + alpha + beta)*(2*degree + alpha + beta - 2)*x + &
        alpha*alpha - beta*beta)*current - &
        2*(degree + alpha - 1)*(degree + beta - 1)* &
        (2*degree + alpha + beta)*previous
    roots(1) = numerator/denominator

    spec%name = str("jacobi_recurrence_kernel")
    spec%module_name = str("fortnum_generated_jacobi_recurrence")
    spec%mode = KERNEL_SUBROUTINE
    spec%temp_prefix = str("t")
    spec%generator = str("gen_jacobi_recurrence")
    spec%generator_revision = str(fortsym_revision())
    spec%regenerate_command = str( &
        "cd tools/codegen && fo exec gen_jacobi_recurrence")
    spec%pure_procedure = .true.
    spec%elemental_procedure = .true.
    spec%openmp_declare_target = .true.
    spec%openacc_routine_seq = .true.
    arguments = [ &
        str("degree"), str("alpha"), str("beta"), str("x"), &
        str("previous"), str("current")]
    outputs(1) = str("next")
    spec%args = arguments
    spec%outputs = outputs

    filename = generated_path("fortnum_jacobi_recurrence_kernel.f90")
    open (newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write "//filename
    code = chars(emit_kernel(roots, spec))
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)
    call codegen_log("wrote "//filename)
end program gen_jacobi_recurrence
