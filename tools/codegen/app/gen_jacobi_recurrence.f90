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
    type(expr_t) :: alpha, beta, current, current_scale, current_x
    type(expr_t) :: current_scale_scale, current_x_scale, current_x_x
    type(expr_t) :: degree, denominator, factor, linear_scale, linear_x
    type(expr_t) :: numerator, previous, previous_factor
    type(expr_t) :: previous_scale, previous_x, roots(6), scale, x
    type(expr_t) :: previous_scale_scale, previous_x_scale, previous_x_x
    type(kernel_spec_t) :: spec
    type(str_t) :: arguments(17), outputs(6)
    character(:), allocatable :: code, filename
    integer :: ios, unit

    call arena%init()
    degree = sym(arena, "degree")
    alpha = sym(arena, "alpha")
    beta = sym(arena, "beta")
    x = sym(arena, "x")
    scale = sym(arena, "scale")
    previous = sym(arena, "previous")
    current = sym(arena, "current")
    previous_x = sym(arena, "previous_x")
    current_x = sym(arena, "current_x")
    previous_scale = sym(arena, "previous_scale")
    current_scale = sym(arena, "current_scale")
    previous_x_x = sym(arena, "previous_x_x")
    current_x_x = sym(arena, "current_x_x")
    previous_x_scale = sym(arena, "previous_x_scale")
    current_x_scale = sym(arena, "current_x_scale")
    previous_scale_scale = sym(arena, "previous_scale_scale")
    current_scale_scale = sym(arena, "current_scale_scale")

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
    arguments(:6) = [ &
        str("degree"), str("alpha"), str("beta"), str("x"), &
        str("previous"), str("current")]
    outputs(1) = str("next")
    spec%args = arguments(:6)
    spec%outputs = outputs(:1)

    filename = generated_path("fortnum_jacobi_recurrence_kernel.f90")
    open (newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write "//filename
    code = chars(emit_kernel(roots(:1), spec))
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)
    call codegen_log("wrote "//filename)

    numerator = (2*degree + alpha + beta - 1)* &
        ((2*degree + alpha + beta)*(2*degree + alpha + beta - 2)*x + &
        (alpha*alpha - beta*beta)*scale)*current - &
        2*(degree + alpha - 1)*(degree + beta - 1)* &
        (2*degree + alpha + beta)*scale*scale*previous
    roots(1) = numerator/denominator
    spec%name = str("scaled_jacobi_recurrence_kernel")
    spec%module_name = str("fortnum_generated_scaled_jacobi_recurrence")
    arguments(:7) = [ &
        str("degree"), str("alpha"), str("beta"), str("x"), &
        str("scale"), str("previous"), str("current")]
    spec%args = arguments(:7)
    filename = generated_path( &
        "fortnum_scaled_jacobi_recurrence_kernel.f90")
    open (newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write "//filename
    code = chars(emit_kernel(roots(:1), spec))
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)
    call codegen_log("wrote "//filename)

    factor = 2*degree + alpha + beta - 1
    linear_x = (2*degree + alpha + beta)* &
        (2*degree + alpha + beta - 2)
    linear_scale = alpha*alpha - beta*beta
    previous_factor = 2*(degree + alpha - 1)*(degree + beta - 1)* &
        (2*degree + alpha + beta)
    roots(2) = (factor*(linear_x*current + &
        (linear_x*x + linear_scale*scale)*current_x) - &
        previous_factor*scale*scale*previous_x)/denominator
    roots(3) = (factor*(linear_scale*current + &
        (linear_x*x + linear_scale*scale)*current_scale) - &
        previous_factor*(2*scale*previous + &
        scale*scale*previous_scale))/denominator
    spec%name = str("scaled_jacobi_gradient_recurrence_kernel")
    spec%module_name = str( &
        "fortnum_generated_scaled_jacobi_gradient_recurrence")
    arguments(:11) = [ &
        str("degree"), str("alpha"), str("beta"), str("x"), &
        str("scale"), str("previous"), str("current"), &
        str("previous_x"), str("current_x"), str("previous_scale"), &
        str("current_scale")]
    outputs(:3) = [str("next"), str("next_x"), str("next_scale")]
    spec%args = arguments(:11)
    spec%outputs = outputs(:3)
    filename = generated_path( &
        "fortnum_scaled_jacobi_gradient_recurrence_kernel.f90")
    open (newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write "//filename
    code = chars(emit_kernel(roots(:3), spec))
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)
    call codegen_log("wrote "//filename)

    roots(4) = (factor*(2*linear_x*current_x + &
        (linear_x*x + linear_scale*scale)*current_x_x) - &
        previous_factor*scale*scale*previous_x_x)/denominator
    roots(5) = (factor*(linear_x*current_scale + &
        linear_scale*current_x + &
        (linear_x*x + linear_scale*scale)*current_x_scale) - &
        previous_factor*(2*scale*previous_x + &
        scale*scale*previous_x_scale))/denominator
    roots(6) = (factor*(2*linear_scale*current_scale + &
        (linear_x*x + linear_scale*scale)*current_scale_scale) - &
        previous_factor*(2*previous + 4*scale*previous_scale + &
        scale*scale*previous_scale_scale))/denominator
    spec%name = str("scaled_jacobi_hessian_recurrence_kernel")
    spec%module_name = str( &
        "fortnum_generated_scaled_jacobi_hessian_recurrence")
    arguments = [ &
        str("degree"), str("alpha"), str("beta"), str("x"), &
        str("scale"), str("previous"), str("current"), &
        str("previous_x"), str("current_x"), str("previous_scale"), &
        str("current_scale"), str("previous_x_x"), str("current_x_x"), &
        str("previous_x_scale"), str("current_x_scale"), &
        str("previous_scale_scale"), str("current_scale_scale")]
    outputs = [str("next"), str("next_x"), str("next_scale"), &
        str("next_x_x"), str("next_x_scale"), str("next_scale_scale")]
    spec%args = arguments
    spec%outputs = outputs
    filename = generated_path( &
        "fortnum_scaled_jacobi_hessian_recurrence_kernel.f90")
    open (newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write "//filename
    code = chars(emit_kernel(roots, spec))
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)
    call codegen_log("wrote "//filename)
end program gen_jacobi_recurrence
