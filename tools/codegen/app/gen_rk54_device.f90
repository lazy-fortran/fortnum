program gen_rk54_device
    !> Device kernels for the Cash-Karp and Dormand-Prince 5(4) pairs.
    !>
    !> Each method is stated once, as its published Butcher tableau in exact
    !> rationals. Everything else is derived: the stage expressions from the
    !> coupling matrix, the embedded error weights from b - bhat, and the FSAL
    !> property from the last stage row. No coefficient is written twice and no
    !> arithmetic is done by hand.
    !>
    !> Before anything is emitted the tableaux are put through fortsym's order
    !> conditions. A method that does not attain the order it claims, or whose
    !> rows do not sum to their nodes, stops the build. That check is exact:
    !> the residuals are rationals through FLINT, so it is a proof and not a
    !> tolerance.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortsym_arena, only: arena_t
    use fortsym_exact, only: exact_to_real
    use fortsym_expr, only: expr_t, sym, real_expr, operator(+), operator(*)
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_rk, only: butcher_t, rk_from_rows, rk_is_consistent, &
        rk_attains_order, rk_error_weights, rk_is_fsal, rk_is_zero
    use fortsym_string, only: str_t, str, chars
    use fortnum_codegen_provenance, only: codegen_log, generated_path
    implicit none

    type(arena_t), target :: arena
    type(butcher_t) :: cash_karp, dormand_prince

    call arena%init()

    call build_cash_karp(cash_karp)
    call verify(cash_karp, "Cash-Karp", 5, 4, expect_fsal=.false.)
    call emit_method(cash_karp, "ck")

    call build_dormand_prince(dormand_prince)
    call verify(dormand_prince, "Dormand-Prince", 5, 4, expect_fsal=.true.)
    call emit_method(dormand_prince, "dp")

contains

    !> Cash, Karp, ACM TOMS 16(3), 201-222 (1990).
    subroutine build_cash_karp(tableau)
        type(butcher_t), intent(out) :: tableau
        logical :: ok

        call rk_from_rows(6, &
            [character(len=16) :: &
             "1/5", &
             "3/40", "9/40", &
             "3/10", "-9/10", "6/5", &
             "-11/54", "5/2", "-70/27", "35/27", &
             "1631/55296", "175/512", "575/13824", "44275/110592", &
             "253/4096"], &
            [character(len=16) :: "37/378", "0", "250/621", "125/594", "0", &
             "512/1771"], &
            [character(len=16) :: "0", "1/5", "3/10", "3/5", "1", "7/8"], &
            tableau, ok, &
            bhat=[character(len=16) :: "2825/27648", "0", "18575/48384", &
                  "13525/55296", "277/14336", "1/4"])
        if (.not. ok) error stop "Cash-Karp tableau is not exact rational data"
    end subroutine build_cash_karp

    !> Dormand, Prince, J. Comp. Appl. Math. 6(1), 19-26 (1980).
    subroutine build_dormand_prince(tableau)
        type(butcher_t), intent(out) :: tableau
        logical :: ok

        call rk_from_rows(7, &
            [character(len=16) :: &
             "1/5", &
             "3/40", "9/40", &
             "44/45", "-56/15", "32/9", &
             "19372/6561", "-25360/2187", "64448/6561", "-212/729", &
             "9017/3168", "-355/33", "46732/5247", "49/176", "-5103/18656", &
             "35/384", "0", "500/1113", "125/192", "-2187/6784", "11/84"], &
            [character(len=16) :: "35/384", "0", "500/1113", "125/192", &
             "-2187/6784", "11/84", "0"], &
            [character(len=16) :: "0", "1/5", "3/10", "4/5", "8/9", "1", "1"], &
            tableau, ok, &
            bhat=[character(len=16) :: "5179/57600", "0", "7571/16695", &
                  "393/640", "-92097/339200", "187/2100", "1/40"])
        if (.not. ok) error stop "Dormand-Prince tableau is not exact rational data"
    end subroutine build_dormand_prince

    !> Refuse to emit a method that is not what it claims to be.
    subroutine verify(tableau, name, order, embedded_order, expect_fsal)
        type(butcher_t), intent(in) :: tableau
        character(*), intent(in) :: name
        integer, intent(in) :: order, embedded_order
        logical, intent(in) :: expect_fsal

        if (.not. rk_is_consistent(tableau)) then
            error stop name//": stage rows do not sum to their nodes"
        end if
        if (.not. rk_attains_order(tableau, tableau%b, order)) then
            error stop name//": weights do not attain the claimed order"
        end if
        if (.not. rk_attains_order(tableau, tableau%bhat, embedded_order)) then
            error stop name//": embedded weights do not attain their order"
        end if
        ! An embedded pair whose two solutions agreed to the full order would
        ! estimate nothing, so this has to fail.
        if (rk_attains_order(tableau, tableau%bhat, order)) then
            error stop name//": embedded weights match the main order"
        end if
        if (rk_is_fsal(tableau) .neqv. expect_fsal) then
            error stop name//": FSAL does not match the method"
        end if
        call codegen_log(name//": consistent, order verified, FSAL "// &
                         merge("yes", "no ", expect_fsal))
    end subroutine verify

    subroutine emit_method(tableau, prefix)
        type(butcher_t), intent(in) :: tableau
        character(*), intent(in) :: prefix
        integer :: stage

        do stage = 2, tableau%stages
            call emit_stage(tableau, prefix, stage)
        end do
        call emit_finish(tableau, prefix)
    end subroutine emit_method

    !> ystage = y + h * sum_j a(stage, j) k_j.
    subroutine emit_stage(tableau, prefix, stage)
        type(butcher_t), intent(in) :: tableau
        character(*), intent(in) :: prefix
        integer, intent(in) :: stage
        type(expr_t) :: roots(1)
        type(str_t), allocatable :: names(:)
        character(len=2) :: stage_text

        roots(1) = sym(arena, "y") + scaled_sum(tableau%a(stage, :))
        names = argument_names(reshape(tableau%a(stage, :), &
                                       [tableau%stages, 1]))
        write (stage_text, "(i0)") stage
        call write_kernel(prefix//"_stage"//trim(stage_text), names, &
                          [str("ystage")], roots)
    end subroutine emit_stage

    !> ynew from b, yerror from the derived b - bhat.
    subroutine emit_finish(tableau, prefix)
        type(butcher_t), intent(in) :: tableau
        character(*), intent(in) :: prefix
        type(expr_t) :: roots(2)
        type(str_t), allocatable :: error_weights(:), names(:)
        logical :: ok

        error_weights = rk_error_weights(tableau, ok)
        if (.not. ok) error stop prefix//": embedded weights are missing"

        roots(1) = sym(arena, "y") + scaled_sum(tableau%b)
        roots(2) = scaled_sum(error_weights)
        names = argument_names(reshape([tableau%b, error_weights], &
                                       [tableau%stages, 2]))
        call write_kernel(prefix//"_finish", names, &
                          [str("ynew"), str("yerror")], roots)
    end subroutine emit_finish

    !> y, h, then every stage any of these weight vectors actually uses, in
    !> ascending stage order. Consumers bind these positionally, so the order
    !> has to depend only on the tableau and not on which expression named a
    !> stage first.
    function argument_names(weight_sets) result(names)
        type(str_t), intent(in) :: weight_sets(:, :)
        type(str_t), allocatable :: names(:)
        integer :: j, w
        logical :: used

        names = [str("y"), str("h")]
        do j = 1, size(weight_sets, 1)
            used = .false.
            do w = 1, size(weight_sets, 2)
                if (.not. is_zero_weight(weight_sets(j, w))) used = .true.
            end do
            if (used) names = [names, str("k"//integer_text(j))]
        end do
    end function argument_names

    !> Exact zeros must not change spelling with the compiler or FLINT build.
    !> Convert only for this boolean classification: zero is represented
    !> exactly as 0.0_dp, while every RK54 coefficient is far from underflow.
    logical function is_zero_weight(weight)
        type(str_t), intent(in) :: weight
        real(dp) :: value
        logical :: ok

        value = exact_to_real(chars(weight), ok)
        if (.not. ok) error stop "weight is not an exact rational"
        is_zero_weight = value == 0.0_dp
    end function is_zero_weight

    !> h * sum_j w_j k_j over the non-zero weights, in ascending stage order.
    function scaled_sum(weights) result(expression)
        type(str_t), intent(in) :: weights(:)
        type(expr_t) :: expression
        type(expr_t) :: total, coefficient
        type(str_t) :: stage_name
        logical :: started, ok
        integer :: j

        started = .false.
        do j = 1, size(weights)
            if (is_zero_weight(weights(j))) cycle
            stage_name = str("k"//integer_text(j))
            ! Round the exact weight to a double here, once, rather than
            ! emitting the rational and leaving a division in the inner loop.
            ! Under strict IEEE a compiler may not fold 37.0_dp/378.0_dp,
            ! because (k*37)/378 and k*(37/378) differ in the last place, so
            ! the emitted rational costs a division per stage per component on
            ! every step. The rounding error is one ulp on a coefficient, far
            ! under the truncation error of the method itself.
            coefficient = real_expr(arena, exact_to_real(chars(weights(j)), ok))
            if (.not. ok) error stop "weight is not an exact rational"
            if (started) then
                total = total + coefficient*sym(arena, chars(stage_name))
            else
                total = coefficient*sym(arena, chars(stage_name))
                started = .true.
            end if
        end do
        if (.not. started) error stop "every weight vanished"
        expression = sym(arena, "h")*total
    end function scaled_sum

    pure function integer_text(value) result(text)
        integer, intent(in) :: value
        character(:), allocatable :: text
        character(len=12) :: buffer

        write (buffer, "(i0)") value
        text = trim(buffer)
    end function integer_text

    subroutine write_kernel(short_name, names, outputs, roots)
        character(*), intent(in) :: short_name
        type(str_t), intent(in) :: names(:), outputs(:)
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
            spec%args(i) = names(i)
        end do
        do i = 1, size(outputs)
            spec%outputs(i) = outputs(i)
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
