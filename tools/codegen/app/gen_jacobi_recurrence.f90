program gen_jacobi_recurrence
    use fortsym_expr, only: expr_t, operator(*), operator(+), &
        operator(-), operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t, words
    use fortnum_codegen_provenance, only: codegen_log, fortsym_revision, &
        generated_path, cost_block_text, insert_cost_block
    implicit none

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

    degree = sym("degree")
    alpha = sym("alpha")
    beta = sym("beta")
    x = sym("x")
    scale = sym("scale")
    previous = sym("previous")
    current = sym("current")
    previous_x = sym("previous_x")
    current_x = sym("current_x")
    previous_scale = sym("previous_scale")
    current_scale = sym("current_scale")
    previous_x_x = sym("previous_x_x")
    current_x_x = sym("current_x_x")
    previous_x_scale = sym("previous_x_scale")
    current_x_scale = sym("current_x_scale")
    previous_scale_scale = sym("previous_scale_scale")
    current_scale_scale = sym("current_scale_scale")

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
    arguments(:6) = words("degree alpha beta x previous current")
    outputs(1) = str("next")
    spec%args = arguments(:6)
    spec%outputs = outputs(:1)

    filename = generated_path("fortnum_jacobi_recurrence_kernel.f90")
    open (newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write "//filename
    code = chars(emit_kernel(roots(:1), spec))
    code = insert_cost_block(code, cost_block_text(roots(:1), roots(:1)))
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
    arguments(:7) = words("degree alpha beta x scale previous current")
    spec%args = arguments(:7)
    filename = generated_path( &
        "fortnum_scaled_jacobi_recurrence_kernel.f90")
    open (newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write "//filename
    code = chars(emit_kernel(roots(:1), spec))
    code = insert_cost_block(code, cost_block_text(roots(:1), roots(:1)))
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
    arguments(:11) = words("degree alpha beta x scale previous current " &
        //"previous_x current_x previous_scale current_scale")
    outputs(:3) = words("next next_x next_scale")
    spec%args = arguments(:11)
    spec%outputs = outputs(:3)
    filename = generated_path( &
        "fortnum_scaled_jacobi_gradient_recurrence_kernel.f90")
    open (newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write "//filename
    code = chars(emit_kernel(roots(:3), spec))
    code = insert_cost_block(code, cost_block_text(roots(:3), roots(:3)))
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
    arguments = words("degree alpha beta x scale previous current " &
        //"previous_x current_x previous_scale current_scale " &
        //"previous_x_x current_x_x previous_x_scale current_x_scale " &
        //"previous_scale_scale current_scale_scale")
    outputs = words("next next_x next_scale next_x_x next_x_scale " &
        //"next_scale_scale")
    spec%args = arguments
    spec%outputs = outputs
    filename = generated_path( &
        "fortnum_scaled_jacobi_hessian_recurrence_kernel.f90")
    open (newunit=unit, file=filename, status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) error stop "cannot write "//filename
    code = chars(emit_kernel(roots, spec))
    code = insert_cost_block(code, cost_block_text(roots, roots))
    write (unit, "(a)") code(:len(code) - 1)
    close (unit)
    call codegen_log("wrote "//filename)
end program gen_jacobi_recurrence
