program gen_eigen_bounds
    use fortsym_string, only: str, str_t, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, num, sqrt, sin, operator(+), operator(-), &
                            operator(*), operator(/), operator(**)
    use fortsym_diff, only: diff
    use fortsym_products, only: jvp, vjp
    use fortsym_engine, only: engine_result_t, VERDICT_TRUE
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE, &
                              PRECISION_REAL64, PRECISION_REAL32
    use fortnum_codegen_provenance, only: fortsym_revision, generated_path, &
                                          cost_block_text, insert_cost_block
    implicit none
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: r, ml, xn, rho, t, a, b, z, d, c, lam, x, k, f, identity, term
    type(expr_t) :: roots(2), vars(4), tangents(4), cot(2)
    type(expr_t), allocatable :: forward(:), reverse(:), products(:)
    integer :: precision, j
    call arena%init()
    engine = make_symengine_engine(arena)
    r = sym(arena, 'r')
    ml = sym(arena, 'ml')
    xn = sym(arena, 'xn')
    rho = r/(sqrt(ml)*xn)
    t = sym(arena, 'theta')
    a = sym(arena, 'left')
    b = sym(arena, 'right')
    ! Verify exact algebraic identities, never accepting a numerical probe.
    call zero('metric residual scaling', rho*rho*ml*xn*xn - r*r)
    roots = [t - r*r/(b - t), t + r*r/(t - a)]
    call zero('Temple lower', (t - roots(1))*(b - t) - r*r)
    call zero('Temple upper', (roots(2) - t)*(t - a) - r*r)
    ! Toeplitz modes v_j=z^j-z^-j: lambda=d+c(z+z^-1).
    ! Dirichlet endpoints give z=exp(i*pi*k/(n+1)), hence lambda=d+2c cos(...).
    z = sym(arena, 'z')
    d = sym(arena, 'd')
    c = sym(arena, 'c')
    lam = d + c*(z + 1/z)
    do j = 1, 6
        identity = c*(z**(j - 1) - z**(1 - j)) + &
                   (d - lam)*(z**j - z**(-j)) + c*(z**(j + 1) - z**(-j - 1))
        call zero('Toeplitz recurrence', identity)
    end do
    ! Continuous interval Laplacian: -d2 sin(k*x)=k^2 sin(k*x).
    x = sym(arena, 'x')
    k = sym(arena, 'k')
    f = sin(k*x)
    call zero('interval Laplacian', diff(diff(f, x), x) + k*k*f)
    ! Axisymmetric disc Laplacian: J0(sqrt(lambda)*r) power series.
    ! The degree-16 partial sum satisfies the ODE up to its final term.
    term = num(arena, 1)
    f = term
    do j = 1, 8
        term = -term*lam*x*x/num(arena, 4*j*j)
        f = f + term
    end do
    call zero('disc Bessel recurrence',diff(diff(f,x),x)+diff(f,x)/x+lam*f-lam*term)
    ! Congruence K=S T S, M=S^2 gives identical characteristic roots.
    ! In the diagonal manufactured pencil (K_ii=M_ii lambda_i), the exact
    ! shift-invert value is 1/(lambda_i-sigma), independently of M_ii.
    call zero('generalized shifted pencil', ml/(ml*lam - t*ml) - 1/(lam - t))
    do precision = PRECISION_REAL64, PRECISION_REAL32
        call emit('radius', [rho], [str('r'), str('ml'), str('xn')], &
                  [str('rho')], precision)
       call emit('temple', roots, [str('r'), str('theta'), str('left'), str('right')], &
                  [str('lower'), str('upper')], precision)
        vars = [r, t, a, b]
        tangents = [sym(arena, 'tr'), sym(arena, 'ttheta'), &
                    sym(arena, 'tleft'), sym(arena, 'tright')]
        cot = [sym(arena, 'ulower'), sym(arena, 'uupper')]
        forward = jvp(roots, vars, tangents)
        reverse = vjp(roots, vars, cot)
        products = [roots, forward, reverse]
        call emit('temple_products', products, &
                  [str('r'), str('theta'), str('left'), str('right'), &
                   str('tr'), str('ttheta'), str('tleft'), str('tright'), &
                   str('ulower'), str('uupper')], &
                  [str('lower'), str('upper'), str('tlower'), str('tupper'), &
                   str('ar'), str('atheta'), str('aleft'), str('aright')], precision)
    end do
    print *, 'eigen bound identities and analytic spectrum identities: symbolic PASS'
contains
    subroutine zero(label, e)
        character(*), intent(in) :: label
        type(expr_t), intent(in) :: e
        type(engine_result_t) :: result
        result = engine%zero_test(e)
        if (.not. result%ok) error stop 'symbolic backend failed'
        if (result%verdict /= VERDICT_TRUE) then
            print *, label
            error stop 'identity not proved symbolically'
        end if
    end subroutine
    subroutine emit(stem, expressions, args, outputs, precision)
        character(*), intent(in) :: stem
        type(expr_t), intent(in) :: expressions(:)
        type(str_t), intent(in) :: args(:), outputs(:)
        integer, intent(in) :: precision
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code, name, suffix
        integer :: unit
        suffix = '_dp'
        if (precision == PRECISION_REAL32) suffix = '_sp'
        name = 'fortnum_eigen_'//stem//suffix
        spec%name = str(name)
        spec%module_name = str(name//'_kernel')
        spec%generator = str('gen_eigen_bounds')
        spec%generator_revision = str(fortsym_revision())
        spec%regenerate_command = str('cd tools/codegen && fo exec gen_eigen_bounds')
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%elemental_procedure = .true.
        spec%precision = precision
        spec%args = args
        spec%outputs = outputs
        spec%temp_prefix = str('t')
        code = chars(emit_kernel(expressions, spec))
        code = insert_cost_block(code, cost_block_text(expressions, expressions))
        open (newunit=unit, file=generated_path(name//'.f90'), status='replace')
        write (unit, '(a)') code(:len(code) - 1)
        close (unit)
    end subroutine
end program
