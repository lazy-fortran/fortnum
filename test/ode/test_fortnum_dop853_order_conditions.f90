program test_fortnum_dop853_order_conditions
    ! Spec-derived invariant check on the generated DOP853 tableau.
    !
    ! The tableau in src/generated/fortnum_dop853_tableau.f90 is derived from
    ! four free nodes (see tools/codegen/src/dop853_construction.f90), not
    ! transcribed. That construction claims an eighth-order method. This test
    ! is an independent, algebraic check that the *generated* coefficients
    ! satisfy the mathematical specification of an order-8 Runge-Kutta method:
    !
    !   * stage consistency:  c_i = sum_j a(i,j)   for every stage i
    !   * weight consistency: sum_i b_i = 1
    !   * the Butcher order conditions: for every rooted tree tau with
    !     |tau| <= 8, the elementary weight b(tau) equals 1/gamma(tau).
    !
    ! The elementary weights are computed from the generated coefficients by
    ! the standard recursive definition (Hairer, Norsett & Wanner, Solving
    ! Ordinary Differential Equations I, section II.1), and the exact values
    ! 1/gamma(tau) come from the tree structure alone -- not from the
    ! construction. If a generated tableau entry had a transcription or
    ! rounding error that changed any order condition, this test fails before
    ! any integration is run. This is the "generated diagnostic invariant /
    ! convergence check derived from the specification" half of the synthesis
    ! pilot.
    !
    ! Trees are represented as canonical strings: "()" for a leaf and
    ! "(" ++ nondecreasing(child keys) ++ ")" for an internal node. The tree
    ! generator is validated against the number of unlabeled rooted trees by
    ! node count (OEIS A000081): 1, 1, 2, 4, 9, 20, 48, 115 for n = 1..8.
    ! That count is an independent oracle for the generator, not for the
    ! tableau.

    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_dop853_tableau, only: dop853_c, dop853_b, dop853_a, &
        DOP853_STAGES
    implicit none

    integer, parameter :: MAX_ORDER = 8
    integer, parameter :: KNOWN_TREE_COUNTS(MAX_ORDER) = [1, 1, 2, 4, 9, 20, &
        48, 115]
    integer, parameter :: MAXKEY = 64

    type :: str_list_t
        character(len=MAXKEY) :: items(1024)
        integer :: n = 0
    end type str_list_t

    type(str_list_t) :: trees(MAX_ORDER)
    real(dp) :: tol
    integer :: n, i, nfail, ntrees

    nfail = 0
    ! The generated coefficients are correctly-rounded doubles of exact
    ! elements of Q(sqrt 6), so each carries up to about an ulp of error.
    ! An elementary weight is a sum of products of up to 7 coefficients; the
    ! accumulated floating-point error is a few ulps. A genuine transcription
    ! error changes an order condition by O(1), many orders above this bound.
    tol = 1.0e-12_dp

    call check_stage_consistency(nfail, tol)
    call check_weight_consistency(nfail, tol)

    ! Generate all rooted trees of each size, deduplicated by canonical form.
    do n = 1, MAX_ORDER
        call gen_trees(n, trees(n))
        call dedupe(trees(n))
        if (trees(n)%n /= KNOWN_TREE_COUNTS(n)) then
            write (error_unit, "(a,i0,a,i0,a,i0)") &
                "FAIL tree enumeration: size ", n, " count ", trees(n)%n, &
                " expected ", KNOWN_TREE_COUNTS(n)
            nfail = nfail + 1
        end if
        do i = 1, trees(n)%n
            if (tree_node_count(trees(n)%items(i)) /= n) then
                write (error_unit, "(a,i0,a,i0,a,a)") &
                    "FAIL node-count mismatch: order ", n, " got ", &
                    tree_node_count(trees(n)%items(i)), " tree=", &
                    trim(trees(n)%items(i))
                nfail = nfail + 1
            end if
            call check_tree_weight(trees(n)%items(i), tol, nfail)
        end do
        write (*, "(a,i0,a,i0,a,i0)") "order ", n, &
            " (trees=", trees(n)%n, ") conditions checked"
    end do

    ntrees = 0
    do n = 1, MAX_ORDER
        ntrees = ntrees + trees(n)%n
    end do
    write (*, "(a,i0)") "total rooted trees checked: ", ntrees

    if (nfail /= 0) then
        write (error_unit, "(i0,a)") nfail, &
            " DOP853 order-condition check(s) failed"
        error stop 1
    end if
    write (*, "(a)") "fortnum_dop853_order_conditions: all passed"

contains

    !> Generate every rooted tree with exactly n nodes as a canonical string.
    !> A tree of n nodes is "(" ++ multiset of child canonical keys ++ ")".
    recursive subroutine gen_trees(n, out)
        integer, intent(in) :: n
        type(str_list_t), intent(out) :: out
        integer :: parts(MAX_ORDER)
        type(str_list_t) :: cand

        call init_list(out)
        if (n <= 1) then
            call push_str(out, "()")
            return
        end if

        parts = 0
        call init_list(cand)
        call enumerate_partitions(n - 1, 1, parts, 0, cand)
        out = cand
    end subroutine gen_trees

    !> Enumerate partitions of total into nondecreasing parts >= minpart,
    !> recorded in parts(1:nparts), and for each emit every cartesian
    !> assignment of concrete child trees to the slots.
    recursive subroutine enumerate_partitions(total, minpart, parts, nparts, out)
        integer, intent(in) :: total, minpart
        integer, intent(inout) :: parts(:)
        integer, intent(in) :: nparts
        type(str_list_t), intent(inout) :: out
        integer :: part

        if (total == 0) then
            call emit_assignments(parts, nparts, out)
            return
        end if

        do part = minpart, total
            parts(nparts + 1) = part
            call enumerate_partitions(total - part, part, parts, nparts + 1, out)
        end do
    end subroutine enumerate_partitions

    !> For each slot in parts(1:nparts), pick a concrete child tree of that
    !> size (all orderings), and assemble the parent canonical key. The caller
    !> deduplicates.
    recursive subroutine emit_assignments(parts, nparts, out)
        integer, intent(in) :: parts(:)
        integer, intent(in) :: nparts
        type(str_list_t), intent(inout) :: out
        character(len=MAXKEY) :: chosen(MAX_ORDER)

        if (nparts == 0) then
            call push_str(out, "()")
            return
        end if

        call emit_slot(parts, nparts, 1, chosen, out)
    end subroutine emit_assignments

    !> Recursively choose a concrete tree for slot `slot` of parts(1:nparts),
    !> accumulating the chosen child keys in chosen(:). When all slots are
    !> assigned, assemble the parent canonical key and push it.
    recursive subroutine emit_slot(parts, nparts, slot, chosen, out)
        integer, intent(in) :: parts(:)
        integer, intent(in) :: nparts, slot
        character(len=MAXKEY), intent(inout) :: chosen(:)
        type(str_list_t), intent(inout) :: out
        type(str_list_t) :: sized
        character(len=MAXKEY) :: key
        integer :: i

        if (slot > nparts) then
            ! Canonicalize into a LOCAL copy so the shared `chosen` array
            ! (reused across sibling iterations and across recursion levels
            ! with different nparts) is never mutated.
            call build_key(chosen, nparts, key)
            call push_str(out, key)
            return
        end if

        call trees_of_size(parts(slot), sized)
        do i = 1, sized%n
            chosen(slot) = sized%items(i)
            call emit_slot(parts, nparts, slot + 1, chosen, out)
        end do
    end subroutine emit_slot

    !> Cached access to trees(n): regenerate and deduplicate on demand.
    subroutine trees_of_size(n, out)
        integer, intent(in) :: n
        type(str_list_t), intent(out) :: out
        call gen_trees(n, out)
        call dedupe(out)
    end subroutine trees_of_size

    !> Remove duplicate canonical keys, keeping the first occurrence of each.
    subroutine dedupe(list)
        type(str_list_t), intent(inout) :: list
        type(str_list_t) :: unique
        integer :: i, j
        logical :: seen
        call init_list(unique)
        do i = 1, list%n
            seen = .false.
            do j = 1, unique%n
                if (unique%items(j) == list%items(i)) then
                    seen = .true.
                    exit
                end if
            end do
            if (.not. seen) call push_str(unique, list%items(i))
        end do
        list = unique
    end subroutine dedupe

    subroutine check_stage_consistency(nfail, tol)
        integer, intent(inout) :: nfail
        real(dp), intent(in) :: tol
        integer :: i
        real(dp) :: rowsum
        do i = 1, DOP853_STAGES
            rowsum = sum(dop853_a(i, :))
            if (abs(rowsum - dop853_c(i)) > tol) then
                write (error_unit, "(a,i0,a,es12.4,a,es12.4)") &
                    "FAIL stage consistency c_", i, " sum=", rowsum, &
                    " c=", dop853_c(i)
                nfail = nfail + 1
            end if
        end do
    end subroutine check_stage_consistency

    subroutine check_weight_consistency(nfail, tol)
        integer, intent(inout) :: nfail
        real(dp), intent(in) :: tol
        if (abs(sum(dop853_b) - 1.0_dp) > tol) then
            write (error_unit, "(a,es12.4)") "FAIL weight sum ", sum(dop853_b)
            nfail = nfail + 1
        end if
    end subroutine check_weight_consistency

    !> Check b(tau) = 1/gamma(tau) for one tree given by its canonical key.
    subroutine check_tree_weight(key, tol, nfail)
        character(*), intent(in) :: key
        real(dp), intent(in) :: tol
        integer, intent(inout) :: nfail
        real(dp) :: b, gamma, exact, err
        b = elementary_weight(key)
        gamma = tree_density(key)
        exact = 1.0_dp / gamma
        err = abs(b - exact)
        if (err > tol) then
            write (error_unit, "(a,es12.4,a,es12.4,a,es12.4,a,i0,a,a)") &
                "FAIL order condition: b=", b, " exact=", exact, &
                " err=", err, " nodes=", tree_node_count(key), &
                " tree=", trim(key)
            nfail = nfail + 1
        end if
    end subroutine check_tree_weight

    !> Number of nodes in a tree given by its canonical key (parens / 2).
    function tree_node_count(key) result(s)
        character(*), intent(in) :: key
        integer :: s, i
        s = 0
        do i = 1, len_trim(key)
            if (key(i:i) == '(') s = s + 1
        end do
    end function tree_node_count

    !> b(tau) = sum_i b_i * phi_i(tau), with phi_i from the canonical key.
    function elementary_weight(key) result(b)
        character(*), intent(in) :: key
        real(dp) :: b
        integer :: pos
        real(dp) :: phi(DOP853_STAGES)
        pos = 1
        call parse_phi(key, pos, phi)
        b = sum(dop853_b * phi)
    end function elementary_weight

    !> Parse a canonical key starting at pos (which must point to '(') and
    !> return phi_i(tau) = prod over children (sum_j a(i,j) phi_j(child)).
    recursive subroutine parse_phi(key, pos, phi)
        character(*), intent(in) :: key
        integer, intent(inout) :: pos
        real(dp), intent(out) :: phi(DOP853_STAGES)
        real(dp) :: child_phi(DOP853_STAGES), acc
        integer :: i, j

        pos = pos + 1                      ! consume '('
        phi = 1.0_dp
        do while (key(pos:pos) == '(')
            call parse_phi(key, pos, child_phi)
            do i = 1, DOP853_STAGES
                acc = 0.0_dp
                do j = 1, DOP853_STAGES
                    acc = acc + dop853_a(i, j) * child_phi(j)
                end do
                phi(i) = phi(i) * acc
            end do
        end do
        pos = pos + 1                      ! consume ')'
    end subroutine parse_phi

    !> gamma(tau) = |tau| * prod_k gamma(child_k); gamma(leaf) = 1.
    function tree_density(key) result(g)
        character(*), intent(in) :: key
        real(dp) :: g
        integer :: pos, nodes
        pos = 1
        nodes = 0
        call parse_density(key, pos, g, nodes)
    end function tree_density

    recursive subroutine parse_density(key, pos, g, nodes)
        character(*), intent(in) :: key
        integer, intent(inout) :: pos
        real(dp), intent(out) :: g
        integer, intent(out) :: nodes
        real(dp) :: cg
        integer :: cnodes

        pos = pos + 1
        g = 1.0_dp
        nodes = 1
        do while (key(pos:pos) == '(')
            call parse_density(key, pos, cg, cnodes)
            g = g * cg
            nodes = nodes + cnodes
        end do
        pos = pos + 1
        g = g * real(nodes, dp)
    end subroutine parse_density

    !> Build the canonical key "(" ++ sorted(chosen(1:nparts)) ++ ")" into a
    !> local array so the caller's shared `chosen` buffer is left untouched.
    subroutine build_key(chosen, nparts, key)
        character(len=MAXKEY), intent(in) :: chosen(:)
        integer, intent(in) :: nparts
        character(len=MAXKEY), intent(out) :: key
        character(len=MAXKEY) :: sorted(MAX_ORDER), tmp
        integer :: k, j
        do k = 1, nparts
            sorted(k) = chosen(k)
        end do
        do k = 1, nparts - 1
            do j = k + 1, nparts
                if (sorted(j) < sorted(k)) then
                    tmp = sorted(k)
                    sorted(k) = sorted(j)
                    sorted(j) = tmp
                end if
            end do
        end do
        key = "("
        do k = 1, nparts
            key = trim(key) // trim(sorted(k))
        end do
        key = trim(key) // ")"
    end subroutine build_key

    subroutine init_list(list)
        type(str_list_t), intent(out) :: list
        list%n = 0
        list%items = ""
    end subroutine init_list

    subroutine push_str(list, s)
        type(str_list_t), intent(inout) :: list
        character(*), intent(in) :: s
        if (list%n >= size(list%items)) then
            write (error_unit, "(a)") "push_str: list overflow"
            error stop 1
        end if
        list%n = list%n + 1
        list%items(list%n) = s
    end subroutine push_str

end program test_fortnum_dop853_order_conditions
