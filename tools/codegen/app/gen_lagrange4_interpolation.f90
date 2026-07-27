program gen_lagrange4_interpolation
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_string, only: str_t, str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/)
    use fortsym_products, only: jvp, vjp
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, count_operations, &
        operation_count_t, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, &
        make_symengine_engine
    use fortnum_codegen_provenance, only: fortsym_revision, generated_path
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: x, value, basis, u
    type(expr_t) :: samples(4), nodes(4), variables(5), tangents(5)
    type(expr_t) :: values(1), cotangents(1)
    type(expr_t) :: jvp_product(1), vjp_products(5)
    type(expr_t) :: jvp_roots(2), vjp_roots(6)
    integer :: i, j

    call arena%init()
    engine = make_symengine_engine(arena)
    x = sym(arena, "x")
    variables(1) = x
    tangents(1) = sym(arena, "tx")
    nodes = [num(arena, -1), num(arena, 0), num(arena, 1), num(arena, 2)]
    value = num(arena, 0)
    do i = 1, 4
        samples(i) = sym(arena, indexed_name("y", i))
        variables(i + 1) = samples(i)
        tangents(i + 1) = sym(arena, indexed_name("ty", i))
        basis = num(arena, 1)
        do j = 1, 4
            if (j == i) cycle
            basis = basis*(x - nodes(j))/(nodes(i) - nodes(j))
        end do
        value = value + samples(i)*basis
    end do
    values(1) = value
    u = sym(arena, "u")
    cotangents(1) = u
    jvp_product = jvp(values, variables, tangents)
    vjp_products = vjp(values, variables, cotangents)
    jvp_roots = [value, jvp_product(1)]
    vjp_roots = [value, vjp_products]
    call simplify_roots(engine, jvp_roots)
    call simplify_roots(engine, vjp_roots)

    call write_kernel( &
        "fortnum_lagrange4_jvp_kernel.f90", &
        "fortnum_lagrange4_jvp_kernel", &
        "fortnum_generated_lagrange4_jvp", &
        [str("x"), str("y1"), str("y2"), str("y3"), str("y4"), &
        str("tx"), str("ty1"), str("ty2"), str("ty3"), str("ty4")], &
        [str("value"), str("jvp")], jvp_roots)
    call write_kernel( &
        "fortnum_lagrange4_vjp_kernel.f90", &
        "fortnum_lagrange4_vjp_kernel", &
        "fortnum_generated_lagrange4_vjp", &
        [str("x"), str("y1"), str("y2"), str("y3"), str("y4"), str("u")], &
        [str("value"), str("adjoint_x"), str("adjoint_y1"), &
        str("adjoint_y2"), str("adjoint_y3"), str("adjoint_y4")], &
        vjp_roots)

contains

    function indexed_name(prefix, index) result(name)
        character(*), intent(in) :: prefix
        integer, intent(in) :: index
        character(:), allocatable :: name
        character(1) :: index_text

        write (index_text, "(i0)") index
        name = prefix//index_text
    end function indexed_name

    subroutine simplify_roots(symbolic_engine, roots)
        type(symengine_engine_t), intent(inout) :: symbolic_engine
        type(expr_t), intent(inout) :: roots(:)
        type(engine_result_t) :: simplified
        integer :: root_index

        do root_index = 1, size(roots)
            simplified = symbolic_engine%simplify(roots(root_index))
            if (simplified%ok) roots(root_index) = simplified%value
        end do
    end subroutine simplify_roots

    subroutine write_kernel( &
            filename, name, module_name, arguments, outputs, roots)
        character(*), intent(in) :: filename, name, module_name
        type(str_t), intent(in) :: arguments(:), outputs(:)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t) :: spec
        type(operation_count_t) :: operations
        character(:), allocatable :: output, code
        integer :: unit, ios

        spec%name = str(name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str("t")
        spec%generator = str("gen_lagrange4_interpolation")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && fo exec gen_lagrange4_interpolation")
        spec%pure_procedure = .true.
        spec%openmp_declare_target = .true.
        spec%openacc_routine_seq = .true.
        spec%args = arguments
        spec%outputs = outputs

        output = generated_path(filename)
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
    end subroutine write_kernel

end program gen_lagrange4_interpolation
