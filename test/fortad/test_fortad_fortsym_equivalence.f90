program test_fortad_fortsym_equivalence
    !! fortad's derivative kernels against fortsym's, entry by entry.
    !!
    !! fortnum's production kernels are generated symbolically by fortsym. This
    !! checks that the fortad path is an alternative and not a variant: for the
    !! same primal, at the same points, both products must agree to rounding.
    !!
    !! Agreement is not the only gate. A wrong primal in `tools/fortad/kernels`
    !! would make both sides of a fortad-only comparison agree with each other
    !! and with nothing else, so fortsym's kernel - written from an independent
    !! symbolic expression - is the reference. The adjoint identity is checked
    !! as well, which catches the case where both products are consistently
    !! wrong in the same way.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_generated_det2_jvp, only: fortnum_det2_jvp_kernel
    use fortnum_generated_det2_vjp, only: fortnum_det2_vjp_kernel
    use fortnum_fortad_det2_jvp, only: fortnum_det2_jvp_fortad
    use fortnum_fortad_det2_vjp, only: fortnum_det2_vjp_fortad
    use fortnum_generated_det3_jvp, only: fortnum_det3_jvp_kernel
    use fortnum_generated_det3_vjp, only: fortnum_det3_vjp_kernel
    use fortnum_fortad_det3_jvp, only: fortnum_det3_jvp_fortad
    use fortnum_fortad_det3_vjp, only: fortnum_det3_vjp_fortad
    use fortnum_generated_multi_input_p2_jvp, only: fortnum_multi_input_p2_jvp_kernel
    use fortnum_generated_multi_input_p2_vjp, only: fortnum_multi_input_p2_vjp_kernel
    use fortnum_fortad_multi_input_p2_jvp, only: fortnum_multi_input_p2_jvp_fortad
    use fortnum_fortad_multi_input_p2_vjp, only: fortnum_multi_input_p2_vjp_fortad
    use fortnum_generated_multi_input_p4_jvp, only: fortnum_multi_input_p4_jvp_kernel
    use fortnum_generated_multi_input_p4_vjp, only: fortnum_multi_input_p4_vjp_kernel
    use fortnum_fortad_multi_input_p4_jvp, only: fortnum_multi_input_p4_jvp_fortad
    use fortnum_fortad_multi_input_p4_vjp, only: fortnum_multi_input_p4_vjp_fortad
    use fortnum_generated_multi_input_p8_jvp, only: fortnum_multi_input_p8_jvp_kernel
    use fortnum_generated_multi_input_p8_vjp, only: fortnum_multi_input_p8_vjp_kernel
    use fortnum_fortad_multi_input_p8_jvp, only: fortnum_multi_input_p8_jvp_fortad
    use fortnum_fortad_multi_input_p8_vjp, only: fortnum_multi_input_p8_vjp_fortad
    use fortnum_generated_multi_input_p16_jvp, only: fortnum_multi_input_p16_jvp_kernel
    use fortnum_generated_multi_input_p16_vjp, only: fortnum_multi_input_p16_vjp_kernel
    use fortnum_fortad_multi_input_p16_jvp, only: fortnum_multi_input_p16_jvp_fortad
    use fortnum_fortad_multi_input_p16_vjp, only: fortnum_multi_input_p16_vjp_fortad
    implicit none

    integer :: failures

    failures = 0
    call check(failures)

    if (failures == 0) then
        print *, "test_fortad_fortsym_equivalence: all cases passed"
    else
        print *, "test_fortad_fortsym_equivalence: ", failures, " case(s) FAILED"
        error stop 1
    end if

contains

    subroutine check(failures)
        integer, intent(inout) :: failures
        real(dp) :: x(16), v(16), g_ref(16), g_new(16)
        real(dp) :: value_ref, value_new, jvp_ref, jvp_new, u

        u = 1.7_dp

        ! det2
        call fill(x, v, 4)
        call fortnum_det2_jvp_kernel(x(1), x(2), x(3), x(4), v(1), v(2), v(3),  &
            v(4), jvp_ref)
        call fortnum_det2_jvp_fortad(x(1), v(1), x(2), v(2), x(3), v(3), x(4),  &
            v(4), jvp_new)
        call same("det2 jvp", jvp_ref, jvp_new, failures)

        call fortnum_det2_vjp_kernel(x(1), x(2), x(3), x(4), u, g_ref(1), g_ref(2),  &
            g_ref(3), g_ref(4))
        call fortnum_det2_vjp_fortad(x(1), x(2), x(3), x(4), u, g_new(1), g_new(2),  &
            g_new(3), g_new(4))
        call same_vector("det2 vjp", g_ref(1:4), g_new(1:4), failures)
        ! The two products must also agree with each other.
        call same("det2 adjoint identity", &
                  dot_product(g_new(1:4), v(1:4)), jvp_new*u, failures)

        ! det3
        call fill(x, v, 9)
        call fortnum_det3_jvp_kernel(x(1), x(2), x(3), x(4), x(5), x(6), x(7),  &
            x(8), x(9), v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9),  &
            jvp_ref)
        call fortnum_det3_jvp_fortad(x(1), v(1), x(2), v(2), x(3), v(3), x(4),  &
            v(4), x(5), v(5), x(6), v(6), x(7), v(7), x(8), v(8), x(9), v(9),  &
            jvp_new)
        call same("det3 jvp", jvp_ref, jvp_new, failures)

        call fortnum_det3_vjp_kernel(x(1), x(2), x(3), x(4), x(5), x(6), x(7),  &
            x(8), x(9), u, g_ref(1), g_ref(2), g_ref(3), g_ref(4), g_ref(5),  &
            g_ref(6), g_ref(7), g_ref(8), g_ref(9))
        call fortnum_det3_vjp_fortad(x(1), x(2), x(3), x(4), x(5), x(6), x(7),  &
            x(8), x(9), u, g_new(1), g_new(2), g_new(3), g_new(4), g_new(5),  &
            g_new(6), g_new(7), g_new(8), g_new(9))
        call same_vector("det3 vjp", g_ref(1:9), g_new(1:9), failures)
        ! The two products must also agree with each other.
        call same("det3 adjoint identity", &
                  dot_product(g_new(1:9), v(1:9)), jvp_new*u, failures)

        ! multi_input_p2
        call fill(x, v, 2)
        call fortnum_multi_input_p2_jvp_kernel(x(1), x(2), v(1), v(2), value_ref,  &
            jvp_ref)
        call fortnum_multi_input_p2_jvp_fortad(x(1), v(1), x(2), v(2), value_new,  &
            jvp_new)
        call same("multi_input_p2 jvp", jvp_ref, jvp_new, failures)
        call same("multi_input_p2 jvp value", value_ref, value_new, failures)
        call fortnum_multi_input_p2_vjp_kernel(x(1), x(2), u, value_ref, g_ref(1),  &
            g_ref(2))
        call fortnum_multi_input_p2_vjp_fortad(x(1), x(2), value_new, u, g_new(1),  &
            g_new(2))
        call same_vector("multi_input_p2 vjp", g_ref(1:2), g_new(1:2), failures)
        ! The two products must also agree with each other.
        call same("multi_input_p2 adjoint identity", &
                  dot_product(g_new(1:2), v(1:2)), jvp_new*u, failures)

        ! multi_input_p4
        call fill(x, v, 4)
        call fortnum_multi_input_p4_jvp_kernel(x(1), x(2), x(3), x(4), v(1), v(2),  &
            v(3), v(4), value_ref, jvp_ref)
        call fortnum_multi_input_p4_jvp_fortad(x(1), v(1), x(2), v(2), x(3), v(3),  &
            x(4), v(4), value_new, jvp_new)
        call same("multi_input_p4 jvp", jvp_ref, jvp_new, failures)
        call same("multi_input_p4 jvp value", value_ref, value_new, failures)
        call fortnum_multi_input_p4_vjp_kernel(x(1), x(2), x(3), x(4), u,  &
            value_ref, g_ref(1), g_ref(2), g_ref(3), g_ref(4))
        call fortnum_multi_input_p4_vjp_fortad(x(1), x(2), x(3), x(4), value_new,  &
            u, g_new(1), g_new(2), g_new(3), g_new(4))
        call same_vector("multi_input_p4 vjp", g_ref(1:4), g_new(1:4), failures)
        ! The two products must also agree with each other.
        call same("multi_input_p4 adjoint identity", &
                  dot_product(g_new(1:4), v(1:4)), jvp_new*u, failures)

        ! multi_input_p8
        call fill(x, v, 8)
        call fortnum_multi_input_p8_jvp_kernel(x(1), x(2), x(3), x(4), x(5), x(6),  &
            x(7), x(8), v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), value_ref,  &
            jvp_ref)
        call fortnum_multi_input_p8_jvp_fortad(x(1), v(1), x(2), v(2), x(3), v(3),  &
            x(4), v(4), x(5), v(5), x(6), v(6), x(7), v(7), x(8), v(8), value_new,  &
            jvp_new)
        call same("multi_input_p8 jvp", jvp_ref, jvp_new, failures)
        call same("multi_input_p8 jvp value", value_ref, value_new, failures)
        call fortnum_multi_input_p8_vjp_kernel(x(1), x(2), x(3), x(4), x(5), x(6),  &
            x(7), x(8), u, value_ref, g_ref(1), g_ref(2), g_ref(3), g_ref(4),  &
            g_ref(5), g_ref(6), g_ref(7), g_ref(8))
        call fortnum_multi_input_p8_vjp_fortad(x(1), x(2), x(3), x(4), x(5), x(6),  &
            x(7), x(8), value_new, u, g_new(1), g_new(2), g_new(3), g_new(4),  &
            g_new(5), g_new(6), g_new(7), g_new(8))
        call same_vector("multi_input_p8 vjp", g_ref(1:8), g_new(1:8), failures)
        ! The two products must also agree with each other.
        call same("multi_input_p8 adjoint identity", &
                  dot_product(g_new(1:8), v(1:8)), jvp_new*u, failures)

        ! multi_input_p16
        call fill(x, v, 16)
        call fortnum_multi_input_p16_jvp_kernel(x(1), x(2), x(3), x(4), x(5), x(6),  &
            x(7), x(8), x(9), x(10), x(11), x(12), x(13), x(14), x(15), x(16),  &
            v(1), v(2), v(3), v(4), v(5), v(6), v(7), v(8), v(9), v(10), v(11),  &
            v(12), v(13), v(14), v(15), v(16), value_ref, jvp_ref)
        call fortnum_multi_input_p16_jvp_fortad(x(1), v(1), x(2), v(2), x(3), v(3),  &
            x(4), v(4), x(5), v(5), x(6), v(6), x(7), v(7), x(8), v(8), x(9), v(9),  &
            x(10), v(10), x(11), v(11), x(12), v(12), x(13), v(13), x(14), v(14),  &
            x(15), v(15), x(16), v(16), value_new, jvp_new)
        call same("multi_input_p16 jvp", jvp_ref, jvp_new, failures)
        call same("multi_input_p16 jvp value", value_ref, value_new, failures)
        call fortnum_multi_input_p16_vjp_kernel(x(1), x(2), x(3), x(4), x(5), x(6),  &
            x(7), x(8), x(9), x(10), x(11), x(12), x(13), x(14), x(15), x(16), u,  &
            value_ref, g_ref(1), g_ref(2), g_ref(3), g_ref(4), g_ref(5), g_ref(6),  &
            g_ref(7), g_ref(8), g_ref(9), g_ref(10), g_ref(11), g_ref(12),  &
            g_ref(13), g_ref(14), g_ref(15), g_ref(16))
        call fortnum_multi_input_p16_vjp_fortad(x(1), x(2), x(3), x(4), x(5), x(6),  &
            x(7), x(8), x(9), x(10), x(11), x(12), x(13), x(14), x(15), x(16),  &
            value_new, u, g_new(1), g_new(2), g_new(3), g_new(4), g_new(5),  &
            g_new(6), g_new(7), g_new(8), g_new(9), g_new(10), g_new(11),  &
            g_new(12), g_new(13), g_new(14), g_new(15), g_new(16))
        call same_vector("multi_input_p16 vjp", g_ref(1:16), g_new(1:16), failures)
        ! The two products must also agree with each other.
        call same("multi_input_p16 adjoint identity", &
                  dot_product(g_new(1:16), v(1:16)), jvp_new*u, failures)
    end subroutine check

    subroutine fill(x, v, n)
        !! A fixed, spread-out point. Deterministic on purpose: a failure has to
        !! be reproducible without recording a seed.
        real(dp), intent(out) :: x(:), v(:)
        integer, intent(in) :: n
        integer :: i

        do i = 1, n
            x(i) = 0.4_dp + 0.7_dp*sin(1.3_dp*i)
            v(i) = cos(0.9_dp*i) - 0.2_dp
        end do
    end subroutine fill

    subroutine same(label, reference, produced, failures)
        !! Agreement to rounding. The two kernels evaluate the same expression
        !! in a different order, so they are not bit-identical and should not be
        !! required to be.
        character(*), intent(in) :: label
        real(dp), intent(in) :: reference, produced
        integer, intent(inout) :: failures

        if (abs(produced - reference) > 1.0e-13_dp*max(1.0_dp, abs(reference))) then
            print *, "FAIL ", label, ": ", reference, " vs ", produced
            failures = failures + 1
        else
            print *, "pass ", label
        end if
    end subroutine same

    subroutine same_vector(label, reference, produced, failures)
        character(*), intent(in) :: label
        real(dp), intent(in) :: reference(:), produced(:)
        integer, intent(inout) :: failures
        integer :: i
        logical :: bad

        bad = .false.
        do i = 1, size(reference)
            if (abs(produced(i) - reference(i)) > &
                1.0e-13_dp*max(1.0_dp, abs(reference(i)))) then
                print *, "FAIL ", label, " entry ", i, ": ", reference(i), &
                    " vs ", produced(i)
                bad = .true.
            end if
        end do
        if (bad) then
            failures = failures + 1
        else
            print *, "pass ", label
        end if
    end subroutine same_vector

end program test_fortad_fortsym_equivalence
