program gen_bicubic
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
                            operator(*), operator(**)
    use fortsym_diff, only: diff
    use fortsym_products, only: jvp, vjp
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_string, only: str_t, str, chars
    use fortnum_codegen_provenance, only: generated_path, cost_block_text, &
                                          insert_cost_block
    implicit none

    ! Independently validated generator stack; does not change older generators.
    character(40) :: revision
    type(arena_t), target :: arena
    type(native_engine_t) :: engine
    type(expr_t) :: variables(18), directions(18), weights(3)
    type(expr_t) :: c(4, 4), column(4), polynomial, horner, jet(3)
    type(expr_t) :: tangent(3), adjoint(18), identity
    type(str_t) :: names(18), direction_names(18), weight_names(3)
    type(str_t) :: result_names(3), tangent_names(3), adjoint_names(18)
    type(engine_result_t) :: checked
    integer :: i, j, k, unit, lock_unit, ios
    character(2) :: ij

    open (newunit=lock_unit, file='fortsym-bicubic.lock', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'cannot read fortsym-bicubic.lock'
    read (lock_unit, '(a)', iostat=ios) revision
    close (lock_unit)
    if (ios /= 0) error stop 'cannot parse fortsym-bicubic.lock'
    call arena%init()
    engine = make_native_engine(arena)
    names(1) = str('x')
    names(2) = str('y')
    do j = 1, 4
        do i = 1, 4
            write (ij, '(2i1)') i, j
            names(2 + i + 4*(j - 1)) = str('c'//ij)
        end do
    end do
    do k = 1, 18
        variables(k) = sym(arena, chars(names(k)))
        direction_names(k) = str('d'//chars(names(k)))
        directions(k) = sym(arena, chars(direction_names(k)))
        adjoint_names(k) = str('bar_'//chars(names(k)))
    end do
    c = reshape(variables(3:), [4, 4])

    ! Authoritative polynomial in increasing x/first-index, y/second-index order.
    polynomial = num(arena, 0)
    do j = 1, 4
        do i = 1, 4
            polynomial = polynomial + c(i, j)*variables(1)**(i - 1) &
                         *variables(2)**(j - 1)
        end do
    end do
    ! Hardware-independent fixed-degree evaluation schedule.
    do j = 1, 4
        column(j) = c(4, j)
        do i = 3, 1, -1
            column(j) = c(i, j) + variables(1)*column(j)
        end do
    end do
    horner = column(4)
    do j = 3, 1, -1
        horner = column(j) + variables(2)*horner
    end do
    jet = [horner, diff(horner, variables(1)), diff(horner, variables(2))]
    call prove(horner - polynomial, 'Horner identity')
    call prove(jet(2) - diff(polynomial, variables(1)), 'x partial identity')
    call prove(jet(3) - diff(polynomial, variables(2)), 'y partial identity')

    do k = 1, 3
        result_names(k) = str('jet'//digit(k))
        tangent_names(k) = str('djet'//digit(k))
        weight_names(k) = str('w'//digit(k))
        weights(k) = sym(arena, chars(weight_names(k)))
    end do
    tangent = jvp(jet, variables, directions)
    adjoint = vjp(jet, variables, weights)
    identity = num(arena, 0)
    do k = 1, 3
        identity = identity + weights(k)*tangent(k)
    end do
    do k = 1, 18
        identity = identity - directions(k)*adjoint(k)
    end do
    call prove(identity, 'contracted JVP/VJP identity')

    open (newunit=unit, file=generated_path('fortnum_bicubic_kernel.f90'), &
          status='replace')
    write (unit, '(a)') '! Generated from tools/codegen/app/gen_bicubic.f90.'
    write (unit, '(a)') '! Tested FortGen: 42a4c6cce9272a6bc92e0539b21f3ffc0c97a068.'
    write (unit, '(a)') 'module fortnum_bicubic_kernel'
    write (unit, '(a)') '    implicit none'
    write (unit, '(a)') '    private'
    write (unit, '(a)') &
        '    public :: bicubic_kernel, bicubic_jvp_kernel, bicubic_vjp_kernel'
    write (unit, '(a)') 'contains'
    call emit('bicubic_kernel', names, result_names, jet)
    call emit('bicubic_jvp_kernel', [names, direction_names], &
              [result_names, tangent_names], [jet, tangent])
    call emit('bicubic_vjp_kernel', [names, weight_names], &
              [result_names, adjoint_names], [jet, adjoint])
    write (unit, '(a)') 'end module fortnum_bicubic_kernel'
    close (unit)
    print '(a)', 'PASS bicubic polynomial/partials and contracted JVP/VJP identities'
contains
    function digit(number) result(text)
        integer, intent(in) :: number
        character(1) :: text
        write (text, '(i1)') number
    end function digit

    subroutine prove(expression, label)
        type(expr_t), intent(in) :: expression
        character(*), intent(in) :: label
        checked = engine%zero_test(expression)
        if (checked%verdict /= VERDICT_TRUE) then
            print '(a)', label
            error stop 'bicubic symbolic identity unproved'
        end if
    end subroutine prove

    subroutine emit(name, args, results, roots)
        character(*), intent(in) :: name
        type(str_t), intent(in) :: args(:), results(:)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_ir_t) :: ir
        type(kernel_emit_spec_t) :: spec
        type(str_t) :: source
        character(:), allocatable :: message
        logical :: ok

        call lower_kernel_ir(roots, ir, ok, message)
        if (.not. ok) error stop 'bicubic IR lowering failed'
        spec%name = str(name)
        spec%args = args
        spec%outputs = results
        spec%temp_prefix = str('t')
        spec%pure_procedure = .true.
        spec%generator = str('gen_bicubic')
        spec%generator_revision = str('fortsym@'//revision)
        spec%regenerate_command = str('cd tools/codegen; fo exec gen_bicubic')
        source = emit_fortran_kernel_ir(ir, spec, ok, message)
        if (.not. ok) error stop 'bicubic Fortran emission failed'
        write (unit, '(a)') insert_cost_block(chars(source), &
                                              cost_block_text(roots, roots))
    end subroutine emit
end program gen_bicubic
