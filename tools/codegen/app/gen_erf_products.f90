program gen_erf_products
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, erf, erfc, &
        operator(*), operator(/)
    use fortsym_parse, only: parse_expr
    use fortsym_products, only: directional_derivative, vjp
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortnum_codegen_provenance, only: codegen_log, fortsym_revision, &
        generated_path, cost_block_text, insert_cost_block
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: x

    call arena%init()
    engine = make_symengine_engine(arena)
    x = sym(arena, "x(1)")
    call generate_products("erf", erf(x))
    call generate_products("erfc", erfc(x))

contains

    subroutine generate_products(stem, value)
        character(*), intent(in) :: stem
        type(expr_t), intent(in) :: value
        type(expr_t) :: variables(1), directions(1), values(1), cotangents(1)
        type(expr_t) :: product(1)

        variables(1) = x
        directions(1) = sym(arena, "v(1)")
        product(1) = directional_derivative(value, variables, directions)
        call simplify_product(product)
        call write_product(stem, "jvp", "v", product)

        values(1) = value
        cotangents(1) = sym(arena, "u(1)")
        product = vjp(values, variables, cotangents)
        call simplify_product(product)
        call write_product(stem, "vjp", "u", product)
    end subroutine generate_products

    subroutine simplify_product(product)
        type(expr_t), intent(inout) :: product(:)
        type(engine_result_t) :: result
        type(expr_t) :: symbolic_prefactor, folded_constant
        character(:), allocatable :: message
        logical :: ok
        integer :: k

        symbolic_prefactor = parse_expr(arena, "2/sqrt(pi)", ok, message)
        if (.not. ok) error stop message
        folded_constant = real_expr(arena, 2.0_dp/sqrt(acos(-1.0_dp)))
        do k = 1, size(product)
            product(k) = product(k)/symbolic_prefactor*folded_constant
            result = engine%simplify(product(k))
            if (result%ok) product(k) = result%value
        end do
    end subroutine simplify_product

    subroutine write_product(stem, product_name, seed_name, product)
        character(*), intent(in) :: stem, product_name, seed_name
        type(expr_t), intent(in) :: product(:)
        type(kernel_spec_t) :: spec
        character(:), allocatable :: filename, procedure_name, module_name
        character(:), allocatable :: code
        integer :: unit, ios

        procedure_name = "fortnum_"//stem//"_"//product_name
        module_name = "fortnum_generated_"//stem//"_"//product_name
        filename = procedure_name//"_kernel.f90"
        spec%name = str(procedure_name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str("t")
        spec%generator = str("gen_erf_products")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && fo exec gen_erf_products")
        spec%openacc_routine_seq = .true.
        spec%openmp_declare_target = .true.
        spec%pure_procedure = .true.
        allocate (spec%args(2), spec%outputs(1))
        allocate (spec%arg_shapes(2), spec%output_shapes(1), &
            spec%output_references(1))
        spec%args(1) = str("x")
        spec%args(2) = str(seed_name)
        spec%arg_shapes = [str("(:)"), str("(:)")]
        spec%outputs(1) = str(product_name)
        spec%output_shapes(1) = str("(:)")
        spec%output_references(1) = str(product_name//"(1)")

        open (newunit=unit, file=generated_path(filename), status="replace", &
            action="write", iostat=ios)
        if (ios /= 0) error stop "cannot write "//filename
        code = chars(emit_kernel(product, spec))
        code = insert_cost_block(code, cost_block_text(product, product))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)
        call codegen_log("wrote "//generated_path(filename))
    end subroutine write_product

end program gen_erf_products
