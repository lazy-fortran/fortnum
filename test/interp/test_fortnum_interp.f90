program test_fortnum_interp
    ! Behavioral tests for fortnum_interp / grid_search.
    !
    ! Invariants:
    !   1. Interior: p(i-1) < xi <= p(i) for the returned i.
    !   2. Clamp low:  xi <= p(nmin)  =>  i = nmin + 1.
    !   3. Clamp high: xi >= p(nmax)  =>  i = nmax.
    !   4. Exact grid point: xi == p(k) for some interior k must return an i
    !      such that p(i-1) < xi <= p(i), i.e. i == k  (right-closed convention).
    !   5. Works for a two-element grid (smallest non-trivial case).
    !   6. Non-uniform spacing: correct index found.
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_interp, only: grid_search

    implicit none

    integer :: nfail
    nfail = 0

    call test_clamp_low(nfail)
    call test_clamp_high(nfail)
    call test_right_endpoint_clamp(nfail)
    call test_interior_uniform(nfail)
    call test_interior_nonuniform(nfail)
    call test_exact_node(nfail)
    call test_two_element(nfail)

    if (nfail > 0) then
        write (error_unit, "(i0,a)") nfail, " test(s) failed"
        stop 1
    end if
    write (*, "(a)") "PASS"
    stop 0

contains

    subroutine check_idx(label, got, expected, nfail)
        character(*), intent(in)    :: label
        integer,      intent(in)    :: got, expected
        integer,      intent(inout) :: nfail

        if (got /= expected) then
            nfail = nfail + 1
            write (error_unit, "(a,a,a,i0,a,i0)") &
                "FAIL [", label, "] got=", got, " expected=", expected
        end if
    end subroutine check_idx

    ! xi well below grid minimum -> clamp to nmin+1.
    subroutine test_clamp_low(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: p(1:5)
        integer  :: i

        p = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
        call grid_search(p, 1, 5, 0.0_dp, i)
        call check_idx("clamp_low", i, 2, nfail)
    end subroutine test_clamp_low

    ! xi at exactly the minimum -> clamp to nmin+1.
    subroutine test_right_endpoint_clamp(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: p(0:4)
        integer  :: i

        p = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
        call grid_search(p, 0, 4, 0.0_dp, i)
        call check_idx("clamp_at_min", i, 1, nfail)
    end subroutine test_right_endpoint_clamp

    ! xi >= grid maximum -> clamp to nmax; also xi == p(nmax).
    subroutine test_clamp_high(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: p(1:5)
        integer  :: i

        p = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
        call grid_search(p, 1, 5, 6.0_dp, i)
        call check_idx("clamp_high", i, 5, nfail)

        call grid_search(p, 1, 5, 5.0_dp, i)
        call check_idx("clamp_at_max", i, 5, nfail)
    end subroutine test_clamp_high

    ! Sweep all interior cells on a uniform 10-point grid.
    subroutine test_interior_uniform(nfail)
        integer, intent(inout) :: nfail

        integer,  parameter :: n = 10
        real(dp) :: p(1:n), xi
        integer  :: k, i
        character(len=32) :: label

        do k = 1, n
            p(k) = real(k, dp)
        end do

        ! Query midpoint of each cell [k, k+1] for k = 1..n-1.
        do k = 1, n - 1
            xi = p(k) + 0.5_dp*(p(k + 1) - p(k))
            call grid_search(p, 1, n, xi, i)
            write (label, "(a,i0)") "interior_cell_", k
            ! xi in (p(k), p(k+1)] => i must be k+1
            call check_idx(trim(label), i, k + 1, nfail)
        end do
    end subroutine test_interior_uniform

    ! Non-uniform grid; verify a selection of cells.
    subroutine test_interior_nonuniform(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: p(0:5)
        integer  :: i

        p = [0.0_dp, 0.1_dp, 1.0_dp, 3.0_dp, 5.0_dp, 10.0_dp]

        ! Between p(0)=0 and p(1)=0.1: expect i=1.
        call grid_search(p, 0, 5, 0.05_dp, i)
        call check_idx("nonunif_cell0", i, 1, nfail)

        ! Between p(2)=1 and p(3)=3: expect i=3.
        call grid_search(p, 0, 5, 2.0_dp, i)
        call check_idx("nonunif_cell2", i, 3, nfail)

        ! Between p(4)=5 and p(5)=10: expect i=5.
        call grid_search(p, 0, 5, 7.5_dp, i)
        call check_idx("nonunif_cell4", i, 5, nfail)
    end subroutine test_interior_nonuniform

    ! xi == p(k) for interior k uses the right-closed convention: i == k.
    subroutine test_exact_node(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: p(1:6)
        integer  :: k, i
        character(len=32) :: label

        p = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp, 7.0_dp, 10.0_dp]

        ! For interior nodes k=2..5, xi==p(k) => i should be k
        ! (p(k-1) < p(k) = xi <= p(k), right-closed).
        do k = 2, 5
            call grid_search(p, 1, 6, p(k), i)
            write (label, "(a,i0)") "exact_node_k=", k
            call check_idx(trim(label), i, k, nfail)
        end do
    end subroutine test_exact_node

    ! Minimal two-element grid: one cell, both clamp cases.
    subroutine test_two_element(nfail)
        integer, intent(inout) :: nfail

        real(dp) :: p(1:2)
        integer  :: i

        p = [0.0_dp, 1.0_dp]

        call grid_search(p, 1, 2, -0.5_dp, i)
        call check_idx("two_el_clamp_low", i, 2, nfail)

        call grid_search(p, 1, 2, 0.5_dp, i)
        call check_idx("two_el_interior", i, 2, nfail)

        call grid_search(p, 1, 2, 1.5_dp, i)
        call check_idx("two_el_clamp_high", i, 2, nfail)
    end subroutine test_two_element

end program test_fortnum_interp
