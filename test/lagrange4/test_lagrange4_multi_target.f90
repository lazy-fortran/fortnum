program test_lagrange4_multi_target
    ! Numerical agreement across the four emission targets of the Lagrange-4
    ! interpolation JVP kernel. One expression is emitted to fortran_cpu,
    ! fortran_openmp_target, fortran_openacc, and cuda. On this runner only
    ! the CPU artifact is compiled and checked; the other targets are emitted
    ! (committed) but recorded as not run, never as passing, because no
    ! offload-capable compiler or CUDA toolkit is available here.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_generated_lagrange4_jvp_cpu, only: &
        fortnum_lagrange4_jvp_kernel_cpu
    implicit none

    real(dp) :: x, y(4), tx, ty(4)
    real(dp) :: value, jvp
    real(dp) :: oracle_value, oracle_jvp
    integer :: i

    ! Same fixed, spread-out inputs as the kernel's own equivalence tests.
    x = 0.4_dp + 0.7_dp*sin(1.3_dp*1)
    do i = 1, 4
        y(i) = cos(0.9_dp*i) - 0.2_dp
    end do
    tx = 0.3_dp
    do i = 1, 4
        ty(i) = sin(0.5_dp*i)
    end do

    call fortnum_lagrange4_jvp_kernel_cpu(x, y(1), y(2), y(3), y(4), &
        tx, ty(1), ty(2), ty(3), ty(4), value, jvp)

    call oracle(x, y, tx, ty, oracle_value, oracle_jvp)

    print "(a)", "fortran_cpu:            run"
    if (abs(value - oracle_value) <= 1.0e-13_dp*max(1.0_dp, abs(oracle_value)) &
        .and. abs(jvp - oracle_jvp) <= 1.0e-13_dp*max(1.0_dp, abs(oracle_jvp))) then
        print "(a,es23.16,a,es23.16)", "  value=", value, "  jvp=", jvp
        print "(a)", "  PASS (agrees with independent oracle)"
    else
        print "(a,es23.16,a,es23.16)", "  FAIL value=", value, " oracle=", oracle_value
        print "(a,es23.16,a,es23.16)", "  FAIL jvp=", jvp, " oracle=", oracle_jvp
        error stop 1
    end if

    ! The OpenMP offload, OpenACC, and CUDA artifacts are emitted and
    ! committed, but the runner lacks an offload-capable compiler and CUDA
    ! toolkit, so their comparisons are recorded as not run.
    print "(a)", "fortran_openmp_target: not run (no offload-capable compiler)"
    print "(a)", "fortran_openacc:       not run (no OpenACC device compiler)"
    print "(a)", "cuda:                  not run (no CUDA toolkit)"

contains

    subroutine oracle(x, y, tx, ty, value, jvp)
        ! Independent closed-form Lagrange-4 interpolation value and
        ! directional derivative. Nodes are -1, 0, 1, 2.
        real(dp), intent(in) :: x, y(4), tx, ty(4)
        real(dp), intent(out) :: value, jvp
        integer :: i, j
        real(dp) :: basis, dbasis

        value = 0.0_dp
        jvp = 0.0_dp
        do i = 1, 4
            basis = 1.0_dp
            do j = 1, 4
                if (j == i) cycle
                basis = basis*(x - node(j))/(node(i) - node(j))
            end do
            value = value + y(i)*basis
            dbasis = 0.0_dp
            do j = 1, 4
                if (j == i) cycle
                dbasis = dbasis + basis/(x - node(j))
            end do
            jvp = jvp + ty(i)*basis + y(i)*dbasis*tx
        end do
    end subroutine oracle

    real(dp) function node(k)
        integer, intent(in) :: k
        real(dp) :: n(4)
        n = [-1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
        node = n(k)
    end function node

end program test_lagrange4_multi_target
