program gen_lagrange4_interpolation
    use fortsym_string, only: str_t, str, chars, itoa, words
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/)
    use fortsym_products, only: jvp, vjp
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, count_operations, &
        operation_count_t, KERNEL_SUBROUTINE
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_kernel_emit, only: kernel_emit_spec_t, &
        emit_fortran_kernel_ir, emit_cuda_device_ir, &
        TARGET_FORTRAN_CPU, TARGET_FORTRAN_OPENMP_TARGET, &
        TARGET_FORTRAN_OPENACC, TARGET_CUDA
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, &
        make_symengine_engine
    use fortnum_codegen_provenance, only: codegen_log, codegen_log_count, &
        fortsym_revision, generated_path, cost_block_text, insert_cost_block
    implicit none

    type(symengine_engine_t) :: engine
    type(expr_t) :: x, value, basis, u
    type(expr_t) :: samples(4), nodes(4), variables(5), tangents(5)
    type(expr_t) :: values(1), cotangents(1)
    type(expr_t) :: jvp_product(1), vjp_products(5)
    type(expr_t) :: jvp_roots(2), vjp_roots(6)
    integer :: i, j

    engine = make_symengine_engine()
    x = sym("x")
    variables(1) = x
    tangents(1) = sym("tx")
    nodes = [num(-1), num(0), num(1), num(2)]
    value = num(0)
    do i = 1, 4
        samples(i) = sym("y"//itoa(i))
        variables(i + 1) = samples(i)
        tangents(i + 1) = sym("ty"//itoa(i))
        basis = num(1)
        do j = 1, 4
            if (j == i) cycle
            basis = basis*(x - nodes(j))/(nodes(i) - nodes(j))
        end do
        value = value + samples(i)*basis
    end do
    values(1) = value
    u = sym("u")
    cotangents(1) = u
    jvp_product = jvp(values, variables, tangents)
    vjp_products = vjp(values, variables, cotangents)
    jvp_roots = [value, jvp_product(1)]
    vjp_roots = [value, vjp_products]
    call simplify_roots(engine, jvp_roots)
    call simplify_roots(engine, vjp_roots)

    ! The VJP leaf keeps the historical dual-annotated spelling so existing
    ! consumers and tests are unchanged.
    call write_kernel( &
        "fortnum_lagrange4_vjp_kernel.f90", &
        "fortnum_lagrange4_vjp_kernel", &
        "fortnum_generated_lagrange4_vjp", &
        words("x y1 y2 y3 y4 u"), &
        words("value adjoint_x adjoint_y1 adjoint_y2 adjoint_y3 adjoint_y4"), &
        vjp_roots)

    ! The historical dual-annotated JVP leaf is retained so existing consumers
    ! and tests keep working while the four per-target artifacts are added.
    call write_kernel( &
        "fortnum_lagrange4_jvp_kernel.f90", &
        "fortnum_lagrange4_jvp_kernel", &
        "fortnum_generated_lagrange4_jvp", &
        words("x y1 y2 y3 y4 tx ty1 ty2 ty3 ty4"), &
        words("value jvp"), jvp_roots)

    ! The JVP kernel is the demonstration kernel for the multi-target claim:
    ! one expression lowered once, then emitted to all four targets.
    call write_kernel_targets( &
        "lagrange4_jvp_kernel", &
        words("x y1 y2 y3 y4 tx ty1 ty2 ty3 ty4"), &
        words("value jvp"), jvp_roots)

contains

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
        code = insert_cost_block(code, cost_block_text(roots, roots))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)

        operations = count_operations(roots)
        call codegen_log("wrote "//output)
        call codegen_log_count("post-CSE structural operations: ", &
            operations%total)
    end subroutine write_kernel

    subroutine write_kernel_targets(stem, arguments, outputs, roots)
        !! Emit one kernel expression to all four targets from a single
        !! lowered IR: fortran_cpu, fortran_openmp_target, fortran_openacc,
        !! and cuda. Four artifacts, one expression, one answer.
        character(*), intent(in) :: stem
        type(str_t), intent(in) :: arguments(:), outputs(:)
        type(expr_t), intent(in) :: roots(:)

        type(kernel_ir_t) :: ir
        logical :: ok
        character(:), allocatable :: message

        call lower_kernel_ir(roots, ir, ok, message)
        if (.not. ok) error stop "cannot lower kernel IR: "//message

        call write_fortran_target("_cpu", "fortran_cpu", &
            TARGET_FORTRAN_CPU, stem, arguments, outputs, ir, roots)
        call write_fortran_target("_openmp_target", "fortran_openmp_target", &
            TARGET_FORTRAN_OPENMP_TARGET, stem, arguments, outputs, ir, roots)
        call write_fortran_target("_openacc", "fortran_openacc", &
            TARGET_FORTRAN_OPENACC, stem, arguments, outputs, ir, roots)
        call write_cuda_target(stem, arguments, outputs, ir, roots)
    end subroutine write_kernel_targets

    subroutine write_fortran_target(suffix, target_tag, target, stem, &
            arguments, outputs, ir, roots)
        character(*), intent(in) :: suffix, target_tag
        integer, intent(in) :: target
        character(*), intent(in) :: stem
        type(str_t), intent(in) :: arguments(:), outputs(:)
        type(kernel_ir_t), intent(in) :: ir
        type(expr_t), intent(in) :: roots(:)

        type(kernel_emit_spec_t) :: spec
        type(str_t) :: source
        logical :: ok
        character(:), allocatable :: message, code
        character(:), allocatable :: module_name, subroutine_name, filename
        integer :: banner_end

        subroutine_name = "fortnum_"//stem//suffix
        module_name = "fortnum_generated_lagrange4_jvp"//suffix
        filename = "fortnum_"//stem//suffix//".f90"

        spec%name = str(subroutine_name)
        spec%args = arguments
        spec%outputs = outputs
        spec%temp_prefix = str("t")
        spec%target = target
        spec%pure_procedure = .true.
        spec%generator = str("gen_lagrange4_interpolation")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && fo exec gen_lagrange4_interpolation")

        source = emit_fortran_kernel_ir(ir, spec, ok, message)
        if (.not. ok) error stop "cannot emit Fortran "//target_tag//": "//message

        ! The IR emitter writes a bare subroutine; wrap it in a module so the
        ! artifact is `use`able by Fortran consumers and tests. The banner
        ! stays on top, where fo looks for the generated-file marker.
        code = insert_cost_block(chars(source), cost_block_text(roots, roots))
        banner_end = banner_length(code)
        code = code(:banner_end)// &
            "module "//module_name//new_line("a")// &
            "    implicit none"//new_line("a")// &
            "    private"//new_line("a")// &
            "    public :: "//subroutine_name//new_line("a")// &
            "contains"//new_line("a")// &
            code(banner_end + 1:)// &
            "end module "//module_name//new_line("a")
        call write_artifact(filename, code)
    end subroutine write_fortran_target

    subroutine write_cuda_target(stem, arguments, outputs, ir, roots)
        character(*), intent(in) :: stem
        type(str_t), intent(in) :: arguments(:), outputs(:)
        type(kernel_ir_t), intent(in) :: ir
        type(expr_t), intent(in) :: roots(:)

        type(kernel_emit_spec_t) :: spec
        type(str_t) :: source
        logical :: ok
        character(:), allocatable :: message, code

        spec%name = str("fortnum_"//stem//"_cuda")
        spec%args = arguments
        spec%outputs = outputs
        spec%temp_prefix = str("t")
        spec%target = TARGET_CUDA
        spec%generator = str("gen_lagrange4_interpolation")
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && fo exec gen_lagrange4_interpolation")

        source = emit_cuda_device_ir(ir, spec, ok, message)
        if (.not. ok) error stop "cannot emit CUDA leaf: "//message

        code = insert_cost_block(chars(source), &
            cost_block_text(roots, roots, "//"), "//")
        call write_artifact("fortnum_"//stem//"_cuda.cu", code)
    end subroutine write_cuda_target

    pure function banner_length(code) result(length)
        character(*), intent(in) :: code
        integer :: length, line_end

        length = 0
        do while (length < len(code))
            if (code(length + 1:length + 1) /= "!") exit
            line_end = index(code(length + 1:), new_line("a"))
            if (line_end == 0) exit
            length = length + line_end
        end do
    end function banner_length

    subroutine write_artifact(filename, code)
        character(*), intent(in) :: filename, code
        character(:), allocatable :: output
        integer :: unit, ios

        output = generated_path(filename)
        open (newunit=unit, file=output, status="replace", action="write", &
            iostat=ios)
        if (ios /= 0) error stop "cannot write "//output
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)
        call codegen_log("wrote "//output)
    end subroutine write_artifact

end program gen_lagrange4_interpolation
