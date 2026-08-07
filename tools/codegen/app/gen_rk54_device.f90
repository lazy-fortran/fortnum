program gen_rk54_device
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_parse, only: parse_expr
    use fortsym_string, only: str, chars
    use fortnum_codegen_provenance, only: codegen_log, generated_path
    implicit none

    type(arena_t), target :: arena

    call arena%init()
    call generate_cash_karp()
    call generate_dormand_prince()

contains

    subroutine generate_cash_karp()
        call write_stage("ck", 2, [character(len=2) :: "y", "h", "k1"], &
            "y+h*(1/5*k1)")
        call write_stage("ck", 3, [character(len=2) :: "y", "h", "k1", "k2"], &
            "y+h*(3/40*k1+9/40*k2)")
        call write_stage("ck", 4, &
            [character(len=2) :: "y", "h", "k1", "k2", "k3"], &
            "y+h*(3/10*k1-9/10*k2+6/5*k3)")
        call write_stage("ck", 5, &
            [character(len=2) :: "y", "h", "k1", "k2", "k3", "k4"], &
            "y+h*(-11/54*k1+5/2*k2-70/27*k3+35/27*k4)")
        call write_stage("ck", 6, &
            [character(len=2) :: "y", "h", "k1", "k2", "k3", "k4", "k5"], &
            "y+h*(1631/55296*k1+175/512*k2+575/13824*k3"// &
            "+44275/110592*k4+253/4096*k5)")
        call write_finish("ck", &
            [character(len=2) :: "y", "h", "k1", "k3", "k4", "k5", "k6"], &
            "y+h*(37/378*k1+250/621*k3+125/594*k4+512/1771*k6)", &
            "h*((37/378-2825/27648)*k1+(250/621-18575/48384)*k3"// &
            "+(125/594-13525/55296)*k4-277/14336*k5"// &
            "+(512/1771-1/4)*k6)")
    end subroutine generate_cash_karp

    subroutine generate_dormand_prince()
        call write_stage("dp", 2, [character(len=2) :: "y", "h", "k1"], &
            "y+h*(1/5*k1)")
        call write_stage("dp", 3, [character(len=2) :: "y", "h", "k1", "k2"], &
            "y+h*(3/40*k1+9/40*k2)")
        call write_stage("dp", 4, &
            [character(len=2) :: "y", "h", "k1", "k2", "k3"], &
            "y+h*(44/45*k1-56/15*k2+32/9*k3)")
        call write_stage("dp", 5, &
            [character(len=2) :: "y", "h", "k1", "k2", "k3", "k4"], &
            "y+h*(19372/6561*k1-25360/2187*k2+64448/6561*k3"// &
            "-212/729*k4)")
        call write_stage("dp", 6, &
            [character(len=2) :: "y", "h", "k1", "k2", "k3", "k4", "k5"], &
            "y+h*(9017/3168*k1-355/33*k2+46732/5247*k3"// &
            "+49/176*k4-5103/18656*k5)")
        call write_stage("dp", 7, &
            [character(len=2) :: "y", "h", "k1", "k2", "k3", "k4", "k5", "k6"], &
            "y+h*(35/384*k1+500/1113*k3+125/192*k4"// &
            "-2187/6784*k5+11/84*k6)")
        call write_finish("dp", &
            [character(len=2) :: "y", "h", "k1", "k3", "k4", "k5", "k6", "k7"], &
            "y+h*(35/384*k1+500/1113*k3+125/192*k4"// &
            "-2187/6784*k5+11/84*k6)", &
            "h*(71/57600*k1-71/16695*k3+71/1920*k4"// &
            "-17253/339200*k5+22/525*k6-1/40*k7)")
    end subroutine generate_dormand_prince

    subroutine write_stage(prefix, stage, names, expression_text)
        character(*), intent(in) :: prefix
        integer, intent(in) :: stage
        character(*), intent(in) :: names(:)
        character(*), intent(in) :: expression_text
        type(expr_t) :: roots(1)
        character(len=2) :: stage_text

        write (stage_text, "(i0)") stage
        roots(1) = parsed(expression_text)
        call write_kernel(prefix//"_stage"//trim(stage_text), names, &
            [character(len=6) :: "ystage"], roots)
    end subroutine write_stage

    subroutine write_finish(prefix, names, solution_text, error_text)
        character(*), intent(in) :: prefix
        character(*), intent(in) :: names(:)
        character(*), intent(in) :: solution_text, error_text
        type(expr_t) :: roots(2)

        roots(1) = parsed(solution_text)
        roots(2) = parsed(error_text)
        call write_kernel(prefix//"_finish", names, &
            [character(len=6) :: "ynew", "yerror"], roots)
    end subroutine write_finish

    function parsed(text) result(expression)
        character(*), intent(in) :: text
        type(expr_t) :: expression
        character(:), allocatable :: message
        logical :: ok

        expression = parse_expr(arena, text, ok, message)
        if (.not. ok) error stop message
    end function parsed

    subroutine write_kernel(short_name, names, outputs, roots)
        character(*), intent(in) :: short_name
        character(*), intent(in) :: names(:), outputs(:)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t) :: spec
        character(:), allocatable :: code, path
        integer :: unit, ios, i

        spec%name = str("fortnum_rk54_"//short_name)
        spec%module_name = str("fortnum_generated_rk54_"//short_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%temp_prefix = str("t")
        spec%generator = str("gen_rk54_device")
        spec%generator_revision = str(rk54_fortsym_revision())
        spec%regenerate_command = str( &
            "cd tools/codegen && fo exec gen_rk54_device")
        spec%openacc_routine_seq = .true.
        spec%openmp_declare_target = .true.
        spec%nvfortran_inline = .true.
        spec%pure_procedure = .true.
        spec%elemental_procedure = .true.
        allocate (spec%args(size(names)), spec%outputs(size(outputs)))
        do i = 1, size(names)
            spec%args(i) = str(trim(names(i)))
        end do
        do i = 1, size(outputs)
            spec%outputs(i) = str(trim(outputs(i)))
        end do

        path = generated_path("fortnum_rk54_"//short_name//"_kernel.f90")
        open (newunit=unit, file=path, status="replace", action="write", &
            iostat=ios)
        if (ios /= 0) error stop "cannot write "//path
        code = chars(emit_kernel(roots, spec))
        write (unit, "(a)") code(:len(code) - 1)
        close (unit)
        call codegen_log("wrote "//path)
    end subroutine write_kernel

    function rk54_fortsym_revision() result(revision)
        character(:), allocatable :: revision
        character(128) :: line
        integer :: unit, ios

        open (newunit=unit, file="fortsym-rk54.lock", status="old", &
            action="read", iostat=ios)
        if (ios /= 0) error stop "cannot read fortsym-rk54.lock"
        read (unit, "(a)", iostat=ios) line
        close (unit)
        if (ios /= 0) error stop "cannot parse fortsym-rk54.lock"
        if (len_trim(line) /= 40) then
            error stop "fortsym-rk54.lock must contain a full SHA"
        end if
        revision = "fortsym@"//trim(line)
    end function rk54_fortsym_revision

end program gen_rk54_device
