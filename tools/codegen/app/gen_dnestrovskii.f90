program gen_dnestrovskii
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, real_expr, exp, erfc, sqrt, operator(+), &
                            operator(-), operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortnum_codegen_provenance, only: generated_path, cost_block_text, &
                                          insert_cost_block, fortsym_revision
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: x, f, g, term
    real(dp) :: c(33), q
    real(dp), parameter :: loss(2) = [2.36327180120735470306422331112_dp, &
                                      0.945308720482941881225689324449_dp]
    integer :: i, j, k, input
    character(40) :: name

    call arena%init()
    x = sym(arena, 'x')
    open (newunit=input, file='data/dnestrovskii.dat', status='old')
    do i = 1, 2
        q = real(i, dp) + 1.5_dp
        do j = 1, 8
            read (input, *) c
            f = real_expr(arena, c(1))
            do k = 2, size(c)
                f = f*x + real_expr(arena, c(k))
            end do
            write (name, '(A,I0,A,I0)') 'fortnum_dnestro_re_', 2*i + 3, '_', j
            call emit(trim(name), f)
        end do
        ! Inverse-argument expansion avoids cancellation on the far negative axis.
        term = -x
        f = term
        do k = 1, 31
            term = term*x*real_expr(arena, q + k - 1)
            f = f + term
        end do
        write (name, '(A,I0)') 'fortnum_dnestro_tail_', 2*i + 3
        call emit(trim(name), f)
        g = -real_expr(arena, loss(i))* &
            x**real_expr(arena, q - 1)*exp(-x)
        write (name, '(A,I0)') 'fortnum_dnestro_im_', 2*i + 3
        call emit(trim(name), g)
    end do
    do i = 1, 2
        q = real(i, dp) + 1.5_dp
        do j = 2, 8
            read (input, *) c
            f = real_expr(arena, c(1))
            do k = 2, size(c)
                f = f*x + real_expr(arena, c(k))
            end do
            write (name, '(A,I0,A,I0)') 'fortnum_dnestro_pos_', 2*i + 3, '_', j
            call emit(trim(name), f)
        end do
        term = x
        f = term
        do k = 1, 31
            term = -term*x*real_expr(arena, q + k - 1)
            f = f + term
        end do
        write (name, '(A,I0)') 'fortnum_dnestro_pos_tail_', 2*i + 3
        call emit(trim(name), f)
    end do
    f = real_expr(arena, 2.0_dp/3.0_dp)*(1 - 2*x + &
                 real_expr(arena, 3.54490770181103205459633496668_dp)*x**real_expr(arena, 1.5_dp)* &
                                         exp(x)*erfc(sqrt(x)))
    call emit('fortnum_dnestro_pos_5_1', f)
    f = real_expr(arena, 0.4_dp)*(1 - x*f)
    call emit('fortnum_dnestro_pos_7_1', f)
    close (input)
contains
    subroutine emit(name, value)
        character(*), intent(in) :: name
        type(expr_t), intent(in) :: value
        type(kernel_spec_t) :: spec
        type(expr_t) :: roots(2)
        character(:), allocatable :: code
        integer :: unit, last
        roots = [value, diff(value, x)]
        spec%name = str(name)
        spec%module_name = str(name//'_m')
        spec%args = [str('x')]
        spec%outputs = [str('value'), str('derivative')]
        spec%generator = str('gen_dnestrovskii')
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str('cd tools/codegen && fo exec gen_dnestrovskii')
        spec%pure_procedure = .true.
        spec%elemental_procedure = .true.
        code = chars(emit_kernel(roots, spec))
        code = insert_cost_block(code, cost_block_text(roots, roots))
        last = len_trim(code)
        do while (last > 0)
            if (code(last:last) /= new_line('a')) exit
            last = last - 1
        end do
        open (newunit=unit, file=generated_path(name//'.f90'), status='replace')
        write (unit, '(A)') code(:last)
        close (unit)
    end subroutine
end program
