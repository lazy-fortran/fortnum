program gen_multi_input_scalar
    use fortsym_string, only: str_t, str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(*), &
        operator(/), operator(**), sin
    use fortsym_products, only: jvp, vjp
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, count_operations, &
        operation_count_t, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortnum_codegen_provenance, only: codegen_log, codegen_log_count, &
        fortsym_revision, generated_path, cost_block_text, insert_cost_block
    implicit none

    integer, parameter :: active_sizes(4) = [2, 4, 8, 16]
    integer :: k

    do k = 1, size(active_sizes)
        call generate_size(active_sizes(k))
    end do

contains

    subroutine generate_size(nactive)
        integer, intent(in) :: nactive
        type(arena_t), target :: arena
        type(symengine_engine_t) :: engine
        type(expr_t), allocatable :: x(:), direction(:), variables(:)
        type(expr_t), allocatable :: tangent(:), cotangent(:)
        type(expr_t), allocatable :: jvp_product(:), vjp_products(:)
        type(expr_t), allocatable :: jvp_roots(:), vjp_roots(:)
        type(expr_t) :: value, sum_x, u
        type(str_t), allocatable :: jvp_args(:), vjp_args(:)
        type(str_t), allocatable :: jvp_outputs(:), vjp_outputs(:)
        character(2) :: size_text
        integer :: i

        call arena%init()
        engine = make_symengine_engine(arena)
        allocate (x(nactive), direction(nactive), variables(nactive))
        allocate (tangent(nactive), cotangent(1))
        allocate (jvp_args(2*nactive), vjp_args(nactive + 1))
        allocate (jvp_outputs(2), vjp_outputs(nactive + 1))
        write (size_text, "(i0)") nactive

        sum_x = sym(arena, "zero")
        value = sym(arena, "zero")
        do i = 1, nactive
            x(i) = sym(arena, indexed_name("x", i))
            direction(i) = sym(arena, indexed_name("v", i))
            if (i == 1) then
                sum_x = x(i)
                value = sin(x(i))
            else
                sum_x = sum_x + x(i)
                value = value + sin(x(i))
            end if
            variables(i) = x(i)
            tangent(i) = direction(i)
            jvp_args(i) = str(indexed_name("x", i))
            jvp_args(nactive + i) = str(indexed_name("v", i))
            vjp_args(i) = str(indexed_name("x", i))
            vjp_outputs(i + 1) = str(indexed_name("adjoint", i))
        end do
        value = value + sum_x**2/2
        u = sym(arena, "u")
        cotangent(1) = u
        vjp_args(nactive + 1) = str("u")
        jvp_outputs = [str("value"), str("jvp")]
        vjp_outputs(1) = str("value")

        jvp_product = jvp([value], variables, tangent)
        vjp_products = vjp([value], variables, cotangent)
        allocate (jvp_roots(2), vjp_roots(nactive + 1))
        jvp_roots = [value, jvp_product(1)]
        vjp_roots(1) = value
        vjp_roots(2:) = vjp_products

        call write_variant( &
            "fortnum_multi_input_p"//trim(size_text)//"_jvp_kernel.f90", &
            "fortnum_multi_input_p"//trim(size_text)//"_jvp_kernel", &
            "fortnum_generated_multi_input_p"//trim(size_text)//"_jvp", &
            jvp_args, jvp_outputs, jvp_roots, engine)
        call write_variant( &
            "fortnum_multi_input_p"//trim(size_text)//"_vjp_kernel.f90", &
            "fortnum_multi_input_p"//trim(size_text)//"_vjp_kernel", &
            "fortnum_generated_multi_input_p"//trim(size_text)//"_vjp", &
            vjp_args, vjp_outputs, vjp_roots, engine)
    end subroutine generate_size

    function indexed_name(prefix, index) result(name)
        character(*), intent(in) :: prefix
        integer, intent(in) :: index
        character(:), allocatable :: name
        character(2) :: index_text

        write (index_text, "(i0)") index
        name = prefix//trim(index_text)
    end function indexed_name

    subroutine write_variant(filename, name, module_name, arguments, outputs, &
            expressions, engine)
        character(*), intent(in) :: filename, name, module_name
        type(str_t), intent(in) :: arguments(:), outputs(:)
        type(expr_t), intent(in) :: expressions(:)
        type(symengine_engine_t), intent(inout) :: engine
        character(:), allocatable :: output, code
        type(expr_t) :: roots(size(expressions)), candidate(size(expressions))
        type(expr_t) :: symbolic(size(expressions))
        type(kernel_spec_t) :: spec
        type(operation_count_t) :: operations, candidate_operations
        type(engine_result_t) :: simplified
        integer :: unit, ios, i

        output = generated_path(filename)
        roots = expressions
        symbolic = roots
        operations = count_operations(roots)
        candidate = roots
        do i = 1, size(candidate)
            simplified = engine%simplify(candidate(i))
            if (simplified%ok) candidate(i) = simplified%value
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
        spec%generator = str("gen_multi_input_scalar")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && fo exec gen_multi_input_scalar")
        spec%pure_procedure = .true.
        spec%openmp_declare_target = .true.
        spec%openacc_routine_seq = .true.
        spec%args = arguments
        spec%outputs = outputs

        open (newunit=unit, file=output, status="replace", action="write", &
            iostat=ios)
        if (ios /= 0) error stop "cannot write "//output
        code = chars(emit_kernel(roots, spec))
        code = insert_cost_block(code, cost_block_text(symbolic, roots))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)

        call codegen_log("wrote "//output)
        call codegen_log_count("post-CSE structural operations: ", &
            operations%total)
    end subroutine write_variant

end program gen_multi_input_scalar
