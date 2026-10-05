program gen_verified_midpoint
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(-), &
        operator(*), operator(/), operator(**)
    use fortsym_string, only: str, chars
    use fortsym_rigorous_emit, only: rigorous_kernel_spec_t, interval_runtime, &
        emit_rigorous_kernel
    use fortnum_codegen_provenance, only: generated_path
    implicit none

    type(arena_t), target :: arena
    type(expr_t) :: left, right, width, value, curvature, acc, center, remainder
    character(40) :: revision
    integer :: unit, lock_unit, ios

    open (newunit=lock_unit, file='fortsym-verified-midpoint.lock', &
          status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'cannot read fortsym-verified-midpoint.lock'
    read (lock_unit, '(a)', iostat=ios) revision
    close (lock_unit)
    if (ios /= 0) error stop 'cannot parse fortsym-verified-midpoint.lock'
    call arena%init()
    left = sym(arena, 'left'); right = sym(arena, 'right')
    width = sym(arena, 'width'); value = sym(arena, 'value')
    curvature = sym(arena, 'curvature'); acc = sym(arena, 'acc')
    center = acc + width*value
    remainder = curvature*width**3/24

    open (newunit=unit, &
          file=generated_path('fortnum_verified_midpoint_kernel.f90'), &
          status='replace', action='write')
    write (unit, '(a)') '! Generator: gen_verified_midpoint; do not edit.'
    write (unit, '(a)') '! Generator revision: fortsym@'//revision
    write (unit, '(a)') '! Ordered finite point endpoints; curvature bounds |f''''|.'
    write (unit, '(a)') '! Midpoint theorem: integral error <= width**3*M2/24.'
    write (unit, '(a)') 'module fortnum_generated_verified_midpoint'
    write (unit, '(a)') 'implicit none'
    write (unit, '(a)') 'private'
    write (unit, '(a)') 'public :: midpoint_geometry, midpoint_accumulate'
    write (unit, '(a)') 'contains'
    call emit('midpoint_geometry', ['left ', 'right'], &
              ['midpoint', 'width   '], [left + (right - left)/2, right - left])
    call emit('midpoint_accumulate', &
              ['acc      ', 'width    ', 'value    ', 'curvature'], &
              ['lower', 'upper'], [center - remainder, center + remainder])
    write (unit, '(a)') 'end module fortnum_generated_verified_midpoint'
    close (unit)

contains

    subroutine emit(name, arguments, outputs, roots)
        character(*), intent(in) :: name, arguments(:), outputs(:)
        type(expr_t), intent(in) :: roots(:)
        type(rigorous_kernel_spec_t) :: spec
        character(:), allocatable :: message
        integer :: j
        logical :: ok

        spec%name = str(name)
        spec%generator = str('gen_verified_midpoint')
        spec%runtime = interval_runtime('fortnum_interval', 'interval_t')
        spec%elemental_procedure = .true.
        allocate (spec%args(size(arguments)), spec%outputs(size(outputs)))
        do j = 1, size(arguments)
            spec%args(j) = str(trim(arguments(j)))
        end do
        do j = 1, size(outputs)
            spec%outputs(j) = str(trim(outputs(j)))
        end do
        block
            character(:), allocatable :: code
            code = chars(emit_rigorous_kernel(roots, spec, ok, message))
            if (.not. ok) error stop message
            write (unit, '(a)') code
        end block
    end subroutine emit

end program gen_verified_midpoint
