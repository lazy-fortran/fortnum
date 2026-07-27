program gen_inverse_jvp
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(-), operator(*)
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, count_operations, &
        operation_count_t, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine

    call arena%init()
    engine = make_symengine_engine(arena)
    call generate_kernel(2, &
        "../../src/generated/fortnum_inv2_jvp_kernel.f90", &
        "fortnum_inv2_jvp_kernel", "fortnum_generated_inv2_jvp")
    call generate_kernel(3, &
        "../../src/generated/fortnum_inv3_jvp_kernel.f90", &
        "fortnum_inv3_jvp_kernel", "fortnum_generated_inv3_jvp")

contains

    subroutine generate_kernel(n, path, name, module_name)
        integer, intent(in) :: n
        character(*), intent(in) :: path, name, module_name
        type(expr_t) :: inverse_entries(n*n), directions(n*n)
        type(expr_t) :: roots(n*n), candidate(n*n)
        type(kernel_spec_t) :: spec
        type(operation_count_t) :: operations, candidate_operations
        type(engine_result_t) :: simplified
        character(:), allocatable :: code
        character(8) :: symbol_name
        integer :: i, j, k, l, q, unit, ios

        do q = 1, n*n
            write (symbol_name, "(a,i0)") "r", q
            inverse_entries(q) = sym(arena, trim(symbol_name))
            write (symbol_name, "(a,i0)") "v", q
            directions(q) = sym(arena, trim(symbol_name))
        end do

        do j = 1, n
            do i = 1, n
                q = index_of(i, j, n)
                roots(q) = -inverse_entries(index_of(i, 1, n))* &
                    directions(index_of(1, 1, n))* &
                    inverse_entries(index_of(1, j, n))
                do l = 1, n
                    do k = 1, n
                        if (k == 1 .and. l == 1) cycle
                        roots(q) = roots(q) - &
                            inverse_entries(index_of(i, k, n))* &
                            directions(index_of(k, l, n))* &
                            inverse_entries(index_of(l, j, n))
                    end do
                end do
            end do
        end do

        operations = count_operations(roots)
        candidate = roots
        do q = 1, size(candidate)
            simplified = engine%simplify(candidate(q))
            if (simplified%ok) candidate(q) = simplified%value
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
        spec%generator = str("gen_inverse_jvp")
        spec%regenerate_command = str( &
            "fpm run -C tools/codegen --target gen_inverse_jvp")
        spec%openacc_routine_seq = .true.
        spec%pure_procedure = .true.
        allocate (spec%args(2*n*n), spec%outputs(n*n))
        do q = 1, n*n
            write (symbol_name, "(a,i0)") "r", q
            spec%args(q) = str(trim(symbol_name))
            write (symbol_name, "(a,i0)") "v", q
            spec%args(n*n + q) = str(trim(symbol_name))
            write (symbol_name, "(a,i0)") "w", q
            spec%outputs(q) = str(trim(symbol_name))
        end do

        open (newunit=unit, file=path, status="replace", action="write", iostat=ios)
        if (ios /= 0) error stop "cannot write "//path
        code = chars(emit_kernel(roots, spec))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)

        write (output_unit, "(a)") "wrote "//path
        write (output_unit, "(a,i0)") "post-CSE structural operations: ", &
            operations%total
    end subroutine generate_kernel

    pure function index_of(i, j, n) result(q)
        integer, intent(in) :: i, j, n
        integer :: q

        q = i + (j - 1)*n
    end function index_of

end program gen_inverse_jvp
