program gen_coulomb_series
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
                            operator(*), operator(/)
    use fortsym_kernel, only: kernel_spec_t, emit_kernel
    use fortsym_string, only: str, chars
    implicit none
    type(arena_t), target :: a
    type(expr_t) :: n, eta, rho, f1, f2, v1, v2, roots(2)
    type(kernel_spec_t) :: spec
    integer :: unit
    call a%init()
    n = sym(a, 'n'); eta = sym(a, 'eta'); rho = sym(a, 'rho')
    f1 = sym(a, 'f1'); f2 = sym(a, 'f2')
    v1 = sym(a, 'v1'); v2 = sym(a, 'v2')
    roots = [(2*eta*rho*f1 - rho*rho*f2)/(n*(n + 1)), &
             (2*eta*rho*v1 - rho*rho*v2 - 2*eta*rho*(2*n - 1)*f1)/(n*(n - 1))]
    spec%name = str('coulomb_series_step')
    spec%generator = str('gen_coulomb_series')
    spec%generator_revision = str('7e8ef7abdd90c49bb5f75e2c119e7d181656bcc6')
    spec%regenerate_command = str('see docs/coulomb.md')
    spec%scalar_type = str('real(real128)')
    spec%pure_procedure = .true.
    spec%args = [str('n'), str('eta'), str('rho'), str('f1'), str('f2'), &
                 str('v1'), str('v2')]
    spec%outputs = [str('fn'), str('vn')]
    open (newunit=unit, file='src/generated/fortnum_coulomb_series_kernel.f90', &
          status='replace')
    write (unit, '(a)') 'module fortnum_coulomb_series_kernel'
    write (unit, '(a)') 'use, intrinsic :: iso_fortran_env, only: real128'
    write (unit, '(a)') 'implicit none'
    write (unit, '(a)') 'contains'
    write (unit, '(a)') chars(emit_kernel(roots, spec))
    write (unit, '(a)') 'end module'
    close (unit)
end program
