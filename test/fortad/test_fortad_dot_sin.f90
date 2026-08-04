program test_fortad_dot_sin
    !! fortnum's acceptance test for a fortad-generated derivative kernel.
    !!
    !! fortnum selects derivative candidates on measured evidence, never on the
    !! name of the mechanism that produced them. A generated `autodiff`
    !! candidate therefore has to clear the same bar as any other: an
    !! independent oracle first, and only then a timing.
    !!
    !! The oracle here is central finite differences with a step-size
    !! convergence check, plus agreement between the scalar and vector forms of
    !! the same tangent. Neither compares fortad against another AD tool.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    ! The derivative kernels come as modules, so the compiler checks these
    ! calls. Only the primal keeps an explicit interface: it is a bare
    ! subroutine by design, being the source the generator reads.
    use fortnum_fortad_dot_sin_jvp, only: fortnum_dot_sin_jvp
    use fortnum_fortad_dot_sin_jvp_v, only: fortnum_dot_sin_jvp_v
    use fortnum_fortad_dot_sin_vjp, only: fortnum_dot_sin_vjp
    implicit none

    interface
        pure subroutine fortnum_dot_sin(n, a, b, s)
            import :: dp
            integer, intent(in) :: n
            real(dp), intent(in) :: a(n), b(n)
            real(dp), intent(out) :: s
        end subroutine fortnum_dot_sin
    end interface

    integer, parameter :: n = 64, n_dir = 4
    real(dp) :: a(n), b(n), a_d(n), b_d(n)
    real(dp) :: av(n_dir, n), bv(n_dir, n), sv(n_dir)
    real(dp) :: s, s_d, sp, sm, fd1, fd2, e1, e2, h
    integer :: i, j
    logical :: bad

    bad = .false.
    do i = 1, n
        a(i) = 0.3_dp + 0.11_dp*i
        b(i) = 0.7_dp + 0.07_dp*i
        a_d(i) = sin(0.9_dp*i)
        b_d(i) = cos(1.3_dp*i)
        do j = 1, n_dir
            av(j, i) = sin(0.9_dp*i + 0.4_dp*j)
            bv(j, i) = cos(1.3_dp*i - 0.2_dp*j)
        end do
    end do

    call fortnum_dot_sin_jvp(n, a, a_d, b, b_d, s, s_d)
    call fortnum_dot_sin(n, a, b, sp)
    if (abs(s - sp) > 1.0e-12_dp*max(1.0_dp, abs(s))) then
        print *, "primal mismatch:", s, sp
        bad = .true.
    end if

    h = 1.0e-5_dp
    call fortnum_dot_sin(n, a + h*a_d, b + h*b_d, sp)
    call fortnum_dot_sin(n, a - h*a_d, b - h*b_d, sm)
    fd1 = (sp - sm)/(2.0_dp*h)
    h = 0.5e-5_dp
    call fortnum_dot_sin(n, a + h*a_d, b + h*b_d, sp)
    call fortnum_dot_sin(n, a - h*a_d, b - h*b_d, sm)
    fd2 = (sp - sm)/(2.0_dp*h)
    e1 = abs(fd1 - s_d)
    e2 = abs(fd2 - s_d)

    if (e2 > 1.0e-6_dp*max(1.0_dp, abs(s_d)) + 1.0e-9_dp) then
        print *, "tangent mismatch: ad =", s_d, " fd =", fd2
        bad = .true.
    end if
    ! The convergence check only means something while truncation error still
    ! dominates. Once the central difference reaches its roundoff floor -
    ! about 1e-11 relative for double precision - halving the step makes the
    ! error worse, and testing that would be testing noise.
    if (e1 > 1.0e-8_dp*max(1.0_dp, abs(s_d)) .and. e2 > 0.40_dp*e1) then
        print *, "no second-order convergence:", e1, e2
        bad = .true.
    end if

    ! Vector mode must reproduce the scalar tangent direction by direction, or
    ! carrying directions together has changed an answer.
    call fortnum_dot_sin_jvp_v(n_dir, n, a, av, b, bv, s, sv)
    do j = 1, n_dir
        call fortnum_dot_sin_jvp(n, a, av(j, :), b, bv(j, :), sp, s_d)
        if (abs(sv(j) - s_d) > 1.0e-13_dp*max(1.0_dp, abs(s_d))) then
            print *, "vector vs scalar mismatch, direction", j, sv(j), s_d
            bad = .true.
        end if
    end do

    ! Reverse mode: one sweep gives the gradient with respect to all 2n
    ! inputs. Checked against the analytical gradient, which for this kernel is
    ! d/da_i = sin(b_i) and d/db_i = a_i cos(b_i), and against the adjoint
    ! identity with the tangent already verified above.
    block
        real(dp) :: a_b(n), b_b(n), s_b, lhs, rhs
        s_b = 1.0_dp
        call fortnum_dot_sin_vjp(n, a, b, s, s_b, a_b, b_b)
        if (maxval(abs(a_b - sin(b))) > 1.0e-13_dp) then
            print *, "reverse d/da mismatch:", maxval(abs(a_b - sin(b)))
            bad = .true.
        end if
        if (maxval(abs(b_b - a*cos(b))) > 1.0e-13_dp) then
            print *, "reverse d/db mismatch:", maxval(abs(b_b - a*cos(b)))
            bad = .true.
        end if
        call fortnum_dot_sin_jvp(n, a, a_d, b, b_d, s, s_d)
        lhs = 0.83_dp*s_d
        s_b = 0.83_dp
        call fortnum_dot_sin_vjp(n, a, b, s, s_b, a_b, b_b)
        rhs = sum(a_b*a_d) + sum(b_b*b_d)
        if (abs(lhs - rhs) > 1.0e-12_dp*max(1.0_dp, abs(lhs))) then
            print *, "adjoint identity violated:", lhs, rhs
            bad = .true.
        end if
    end block

    if (bad) then
        print *, "test_fortad_dot_sin: FAILED"
        error stop 1
    end if
    print *, "test_fortad_dot_sin: passed"

end program test_fortad_dot_sin
