program gen_special_region_candidates
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, exp, log, &
        operator(+), operator(*), operator(**)
    use fortsym_products, only: directional_derivative, vjp
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortnum_codegen_provenance, only: codegen_log, fortsym_revision, &
        generated_path
    implicit none

    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: x, value, m
    type(expr_t) :: variables(1), directions(1), cotangents(1), values(1)
    type(expr_t) :: product(1)
    real(dp) :: gamma_25

    call arena%init()
    engine = make_symengine_engine(arena)
    x = sym(arena, "x")
    gamma_25 = gamma(2.5_dp)
    m = real_expr(arena, gamma_25)*exp(x)*x**(-1.5_dp)
    value = log(real_expr(arena, 1.0_dp) + m*m)
    variables(1) = x
    values(1) = value

    directions(1) = sym(arena, "v")
    product(1) = directional_derivative(value, variables, directions)
    call simplify(product)
    call write_product("jvp", "v", product)

    cotangents(1) = sym(arena, "u")
    product = vjp(values, variables, cotangents)
    call simplify(product)
    call write_product("vjp", "u", product)

contains

    subroutine simplify(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        type(engine_result_t) :: result
        integer :: i

        do i = 1, size(expressions)
            result = engine%simplify(expressions(i))
            if (result%ok) expressions(i) = result%value
        end do
    end subroutine simplify

    subroutine write_product(product_name, seed_name, expressions)
        character(*), intent(in) :: product_name, seed_name
        type(expr_t), intent(in) :: expressions(:)
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code, filename, module_name, procedure_name
        integer :: unit, ios

        procedure_name = "fortnum_hyperg_asymptotic_outer_"//product_name
        module_name = "fortnum_generated_hyperg_asymptotic_outer_"//product_name
        filename = procedure_name//"_kernel.f90"
        spec%name = str(procedure_name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%generator = str("gen_special_region_candidates")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && fo exec gen_special_region_candidates")
        spec%temp_prefix = str("t")
        spec%pure_procedure = .true.
        spec%elemental_procedure = .true.
        spec%openmp_declare_target = .true.
        spec%openacc_routine_seq = .true.
        spec%args = [str("x"), str(seed_name)]
        spec%outputs = [str(product_name)]

        open (newunit=unit, file=generated_path(filename), status="replace", &
            action="write", iostat=ios)
        if (ios /= 0) error stop "cannot write "//filename
        code = chars(emit_kernel(expressions, spec))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)
        call codegen_log("wrote "//generated_path(filename))
    end subroutine write_product

end program gen_special_region_candidates
