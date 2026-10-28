program gen_multifidelity
    use fortsym_string, only: str_t, str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, exp, operator(+), operator(-), &
                            operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortnum_codegen_provenance, only: fortsym_revision, cost_block_text, &
                                          insert_cost_block, generated_path
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: x(3), f
    call arena%init()
    x = [sym(arena, 'x'), sym(arena, 'center'), sym(arena, 'length')]
    f = ((x(1) - x(2))/x(3))**2/num(arena, 2)
    call emit('mf_distance', x, f)
    x = [sym(arena, 'distance'), sym(arena, 'amplitude'), sym(arena, 'unused')]
    f = x(2)**2*exp(-x(1))
    call emit('mf_covariance', x(:2), f)
    x = [sym(arena, 'rho'), sym(arena, 'low'), sym(arena, 'correction')]
    f = x(1)*x(2) + x(3)
    call emit('mf_autoregressive', x, f)
contains
    subroutine emit(name, args, value)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: args(:), value
        type(expr_t) :: roots(size(args) + 1)
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code
        integer :: i, unit
        roots(1) = value
        allocate (spec%args(size(args)), spec%outputs(size(roots)))
        spec%outputs(1) = str('value')
        do i = 1, size(args)
            select case (name)
            case ('mf_distance')
                block
                    type(str_t) :: names(3)
                    names = [str('x'), str('center'), str('length')]
                    spec%args(i) = names(i)
                end block
            case ('mf_covariance')
                block
                    type(str_t) :: names(2)
                    names = [str('distance'), str('amplitude')]
                    spec%args(i) = names(i)
                end block
            case default
                block
                    type(str_t) :: names(3)
                    names = [str('rho'), str('low'), str('correction')]
                    spec%args(i) = names(i)
                end block
            end select
            spec%outputs(i + 1) = str('d_'//chars(spec%args(i)))
            roots(i + 1) = diff(value, args(i))
        end do
        spec%name = str(name)
        spec%module_name = str('fortnum_generated_'//name)
        spec%generator = str('gen_multifidelity')
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str('cd tools/codegen && fo exec gen_multifidelity')
        spec%pure_procedure = .true.
        spec%elemental_procedure = .true.
        code = chars(emit_kernel(roots, spec))
        code = insert_cost_block(code, cost_block_text(roots, roots))
        open (newunit=unit, file=generated_path('fortnum_'//name//'.f90'), &
              status='replace')
        write (unit, '(A)') code(:len(code) - 1)
        close (unit)
    end subroutine
end program
