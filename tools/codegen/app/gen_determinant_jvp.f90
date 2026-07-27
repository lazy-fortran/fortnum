program gen_determinant_jvp
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym
    use fortsym_parse, only: parse_expr
    use fortsym_products, only: directional_derivative
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, count_operations, &
        operation_count_t, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine

    call arena%init()
    engine = make_symengine_engine(arena)
    call generate_det2()
    call generate_det3()

contains

    subroutine generate_det2()
        character(*), parameter :: names(4) = ["a", "b", "c", "d"]
        character(*), parameter :: tangent_names(4) = &
            ["va", "vb", "vc", "vd"]
        type(expr_t) :: variables(4), tangents(4), determinant

        call make_symbols(names, variables)
        call make_symbols(tangent_names, tangents)
        determinant = parsed("a*d - c*b")
        call write_kernel( &
            "../../src/generated/fortnum_det2_jvp_kernel.f90", &
            "fortnum_det2_jvp_kernel", "fortnum_generated_det2_jvp", &
            names, tangent_names, determinant, variables, tangents)
    end subroutine generate_det2

    subroutine generate_det3()
        character(*), parameter :: names(9) = [ &
            "a", "b", "c", "d", "f", "g", "h", "j", "k"]
        character(*), parameter :: tangent_names(9) = [ &
            "ta", "tb", "tc", "td", "tf", "tg", "th", "tj", "tk"]
        type(expr_t) :: variables(9), tangents(9), determinant

        call make_symbols(names, variables)
        call make_symbols(tangent_names, tangents)
        determinant = parsed( &
            "a*(f*k-j*g)-d*(b*k-j*c)+h*(b*g-f*c)")
        call write_kernel( &
            "../../src/generated/fortnum_det3_jvp_kernel.f90", &
            "fortnum_det3_jvp_kernel", "fortnum_generated_det3_jvp", &
            names, tangent_names, determinant, variables, tangents)
    end subroutine generate_det3

    subroutine make_symbols(names, expressions)
        character(*), intent(in) :: names(:)
        type(expr_t), intent(out) :: expressions(:)
        integer :: k

        do k = 1, size(names)
            expressions(k) = sym(arena, names(k))
        end do
    end subroutine make_symbols

    function parsed(text) result(expression)
        character(*), intent(in) :: text
        type(expr_t) :: expression
        character(:), allocatable :: message
        logical :: ok

        expression = parse_expr(arena, text, ok, message)
        if (.not. ok) error stop message
    end function parsed

    subroutine write_kernel(path, name, module_name, names, tangent_names, &
            determinant, variables, tangents)
        character(*), intent(in) :: path, name, module_name
        character(*), intent(in) :: names(:), tangent_names(:)
        type(expr_t), intent(in) :: determinant, variables(:), tangents(:)
        type(expr_t) :: roots(1), candidate(1)
        type(kernel_spec_t) :: spec
        type(operation_count_t) :: operations, candidate_operations
        type(engine_result_t) :: simplified
        character(:), allocatable :: code
        integer :: unit, ios, k

        roots(1) = directional_derivative(determinant, variables, tangents)
        operations = count_operations(roots)
        simplified = engine%simplify(roots(1))
        if (simplified%ok) then
            candidate(1) = simplified%value
            candidate_operations = count_operations(candidate)
            if (candidate_operations%total < operations%total) then
                roots(1) = simplified%value
                operations = candidate_operations
            end if
        end if

        spec%name = str(name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str("t")
        spec%generator = str("gen_determinant_jvp")
        spec%regenerate_command = str( &
            "fpm run -C tools/codegen --target gen_determinant_jvp")
        spec%openacc_routine_seq = .true.
        spec%pure_procedure = .true.
        allocate (spec%args(size(names) + size(tangent_names)), spec%outputs(1))
        do k = 1, size(names)
            spec%args(k) = str(names(k))
        end do
        do k = 1, size(tangent_names)
            spec%args(size(names) + k) = str(tangent_names(k))
        end do
        spec%outputs(1) = str("jvp")

        open (newunit=unit, file=path, status="replace", action="write", iostat=ios)
        if (ios /= 0) error stop "cannot write "//path
        code = chars(emit_kernel(roots, spec))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)

        write (output_unit, "(a)") "wrote "//path
        write (output_unit, "(a,i0)") "post-CSE structural operations: ", &
            operations%total
    end subroutine write_kernel

end program gen_determinant_jvp
