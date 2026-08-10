module dop853_construction
    !> Dormand and Prince's reduced system for an eighth order Runge-Kutta
    !> method, solved exactly in Q(sqrt(d)).
    !>
    !> This is fortnum's construction of one specific method family, written in
    !> fortsym's framework. The generic algebra it stands on -- the quadratic
    !> field, exact linear algebra over it, univariate and multivariate
    !> polynomials, and Sylvester elimination -- lives in fortsym, where it
    !> serves any consumer. The reduced system below is particular to this
    !> family, so it belongs here rather than there.
    !>
    !> Hairer, Norsett and Wanner, Solving Ordinary Differential Equations I,
    !> 2nd edition, section II.5: equations (5.20a-m), (5.25a-c), (5.26),
    !> (5.27). Twelve stages. Four free nodes go in and all 78 remaining
    !> coefficients come out. Nothing is transcribed, fitted, or recovered from
    !> a published decimal.
    !>
    !> Why the field. (5.25a) sets c4 = (6 - sqrt 6)c6/10, so the nodes are
    !> irrational and so is everything built on them. In floating point such a
    !> tableau can only be compared against a tolerance; over Q(sqrt 6) the
    !> order conditions either vanish exactly or the method is wrong.
    !>
    !> How the last step closes. After the weights and the rows up to eight,
    !> two families remain: (5.20i-k) give three equations for the four
    !> unknowns of column four, and three more for column five, so each column
    !> carries one free parameter. Three conditions close the system. Row nine
    !> has four coefficients against the five equations (5.20b-f), so its last
    !> equation constrains the parameters rather than defining a coefficient --
    !> that is the condition (5.25c) is chosen to make solvable -- and (5.20l)
    !> supplies one more for each of k = 4, 5.
    !>
    !> The parameters are carried as polynomial variables throughout, so those
    !> three conditions come out as exact polynomials rather than samples. The
    !> row-nine condition is linear and eliminates one parameter by
    !> substitution; the two remaining conditions are then univariate and their
    !> gcd is the factor they agree on. A linear gcd gives the parameter by one
    !> division, and the tableau never leaves the field.
    use fortsym_quadratic, only: quadratic_t, quad_rational, quad_root, &
        quad_add, quad_sub, quad_mul, quad_div, quad_pow, quad_neg, &
        quad_is_zero
    use fortsym_quadratic_linalg, only: quad_solve, quad_nullspace, &
        SOLVE_OK, SOLVE_UNDERDETERMINED
    use fortsym_quadratic_poly, only: qpoly_t, qpoly_gcd, qpoly_degree, &
        qpoly_linear_root
    use fortsym_quadratic_multipoly, only: qmpoly_t, qmpoly_zero, qmpoly_const, &
        qmpoly_var, qmpoly_add, qmpoly_sub, qmpoly_mul, qmpoly_neg, &
        qmpoly_is_zero, qmpoly_degree, qmpoly_coefficient_of, &
        qmpoly_eval_var, qmpoly_substitute, qmpoly_to_univariate, &
        qmpoly_eval_all
    implicit none
    private

    public :: rk_dp8_t, rk_dormand_prince_8, rk_dp8_error_weights
    public :: RK_DP8_OK, RK_DP8_NODES, RK_DP8_WEIGHTS, RK_DP8_ROWS
    public :: RK_DP8_FAMILY, RK_DP8_CLOSURE, RK_DP8_ARITHMETIC

    integer, parameter :: RK_DP8_OK = 0
    integer, parameter :: RK_DP8_NODES = 1
    integer, parameter :: RK_DP8_WEIGHTS = 2
    integer, parameter :: RK_DP8_ROWS = 3
    integer, parameter :: RK_DP8_FAMILY = 4
    integer, parameter :: RK_DP8_CLOSURE = 5
    integer, parameter :: RK_DP8_ARITHMETIC = 6

    !> Stages of the propagating method.
    integer, parameter, public :: RK_DP8_STAGES = 12
    integer, parameter :: S = RK_DP8_STAGES
    !> The two free parameters of the deferred columns.
    integer, parameter :: NVARS = 2

    type :: rk_dp8_t
        type(quadratic_t) :: c(S)
        type(quadratic_t) :: b(S)
        type(quadratic_t) :: a(S, S)
    end type rk_dp8_t

    !> The construction while the two parameters are still symbolic.
    type :: pending_t
        type(qmpoly_t) :: a(S, S)
    end type pending_t

contains

    !> The two embedded estimators DOP853 uses, from section II.10.
    !>
    !> The section II.5 pair is 8(6), and its sixth order embedding shares the
    !> last weight with the main method: bbar(12) = b(12). Hairer and Wanner
    !> record that this makes the estimate blind to the twelfth stage -- the
    !> only one evaluated at the end of the step -- so a step whose internal
    !> stages have run away is accepted, and a discontinuity just left of a grid
    !> point is missed. DOP853 therefore drops that embedding.
    !>
    !> In its place, a fifth order estimator with bbar(6), bbar(7) and bbar(12)
    !> taken as free parameters at the values the book gives,
    !>
    !>     bbar(6) = b(6)/2 + 1, bbar(7) = b(7)/2 + 9/20, bbar(12) = b(12)/2
    !>
    !> with bbar(2..5) = 0 and the rest fixed by the order five conditions; and
    !> a third order estimator carried by the three nodes c1 = 0, c9 and
    !> c12 = 1 alone. The two together let the code damp an over-estimate.
    !>
    !> Both are pure linear solves against the derived tableau, so they stay in
    !> the field: 9/20 is exact, and no decimal enters.
    subroutine rk_dp8_error_weights(tableau, bbar, btilde, status)
        type(rk_dp8_t), intent(in) :: tableau
        type(quadratic_t), intent(out) :: bbar(S), btilde(S)
        integer, intent(out) :: status

        integer, parameter :: five_idx(5) = [1, 8, 9, 10, 11]
        integer, parameter :: three_idx(3) = [1, 9, 12]
        type(quadratic_t) :: matrix(5, 5), rhs(5), fixed(3), contribution
        type(quadratic_t) :: m3(3, 3), r3(3)
        type(quadratic_t), allocatable :: solution(:)
        integer, parameter :: fixed_idx(3) = [6, 7, 12]
        integer :: q, j, local, d
        logical :: ok
        character(len=32) :: text

        d = tableau%c(1)%radicand
        status = RK_DP8_CLOSURE
        do j = 1, S
            bbar(j) = quad_rational("0", d, ok)
            btilde(j) = quad_rational("0", d, ok)
            if (.not. ok) return
        end do

        fixed(1) = quad_add(quad_div(tableau%b(6), quad_rational("2", d, ok), ok), &
                            quad_rational("1", d, ok), ok)
        fixed(2) = quad_add(quad_div(tableau%b(7), quad_rational("2", d, ok), ok), &
                            quad_rational("9/20", d, ok), ok)
        fixed(3) = quad_div(tableau%b(12), quad_rational("2", d, ok), ok)
        if (.not. ok) return
        do j = 1, 3
            bbar(fixed_idx(j)) = fixed(j)
        end do

        ! Order five for the remaining five weights.
        do q = 1, 5
            do j = 1, 5
                matrix(q, j) = quad_pow(tableau%c(five_idx(j)), q - 1, ok)
                if (.not. ok) return
            end do
            write (text, "(a,i0)") "1/", q
            rhs(q) = quad_rational(trim(text), d, ok)
            if (.not. ok) return
            do j = 1, 3
                contribution = quad_mul(fixed(j), &
                                        quad_pow(tableau%c(fixed_idx(j)), q - 1, ok), ok)
                if (.not. ok) return
                rhs(q) = quad_sub(rhs(q), contribution, ok)
                if (.not. ok) return
            end do
        end do
        call quad_solve(matrix, rhs, solution, local)
        if (local /= SOLVE_OK) return
        do j = 1, 5
            bbar(five_idx(j)) = solution(j)
        end do

        ! Order three on the three retained nodes.
        do q = 1, 3
            do j = 1, 3
                m3(q, j) = quad_pow(tableau%c(three_idx(j)), q - 1, ok)
                if (.not. ok) return
            end do
            write (text, "(a,i0)") "1/", q
            r3(q) = quad_rational(trim(text), d, ok)
            if (.not. ok) return
        end do
        call quad_solve(m3, r3, solution, local)
        if (local /= SOLVE_OK) return
        do j = 1, 3
            btilde(three_idx(j)) = solution(j)
        end do
        status = RK_DP8_OK
    end subroutine rk_dp8_error_weights

    subroutine rk_dormand_prince_8(c7, c8, c10, c11, tableau, status)
        type(quadratic_t), intent(in) :: c7, c8, c10, c11
        type(rk_dp8_t), intent(out) :: tableau
        integer, intent(out) :: status

        type(quadratic_t) :: p4(4), n4(4), p5(4), n5(4), u, v
        type(quadratic_t) :: zero
        logical :: ok
        integer :: i, j, d

        d = c7%radicand
        zero = quad_rational("0", d, ok)
        status = RK_DP8_ARITHMETIC
        if (.not. ok) return
        do i = 1, S
            tableau%c(i) = zero
            tableau%b(i) = zero
            do j = 1, S
                tableau%a(i, j) = zero
            end do
        end do

        call build_nodes(c7, c8, c10, c11, tableau, status)
        if (status /= RK_DP8_OK) return
        call build_weights(tableau, status)
        if (status /= RK_DP8_OK) return
        call build_early_rows(tableau, status)
        if (status /= RK_DP8_OK) return
        call build_tail_links(tableau, status)
        if (status /= RK_DP8_OK) return
        call build_families(tableau, p4, n4, p5, n5, status)
        if (status /= RK_DP8_OK) return
        call close_parameters(tableau, p4, n4, p5, n5, u, v, status)
        if (status /= RK_DP8_OK) return
        call apply_parameters(tableau, p4, n4, p5, n5, u, v, status)
    end subroutine rk_dormand_prince_8

    !> (5.25a), (5.25b), (5.25c).
    subroutine build_nodes(c7, c8, c10, c11, tableau, status)
        type(quadratic_t), intent(in) :: c7, c8, c10, c11
        type(rk_dp8_t), intent(inout) :: tableau
        integer, intent(out) :: status

        type(quadratic_t) :: root, six, ten, c6, s1, s2, s3
        type(quadratic_t) :: numerator, denominator, term
        logical :: ok
        integer :: d

        d = c7%radicand
        status = RK_DP8_NODES
        root = quad_root(d, ok); if (.not. ok) return
        six = quad_rational("6", d, ok); if (.not. ok) return
        ten = quad_rational("10", d, ok); if (.not. ok) return

        c6 = quad_mul(quad_rational("4/3", d, ok), c7, ok); if (.not. ok) return
        tableau%c(1) = quad_rational("0", d, ok)
        tableau%c(6) = c6
        tableau%c(7) = c7
        tableau%c(8) = c8
        tableau%c(10) = c10
        tableau%c(11) = c11
        tableau%c(12) = quad_rational("1", d, ok)
        tableau%c(4) = quad_mul(quad_div(quad_sub(six, root, ok), ten, ok), c6, ok)
        tableau%c(5) = quad_mul(quad_div(quad_add(six, root, ok), ten, ok), c6, ok)
        tableau%c(3) = quad_mul(quad_rational("2/3", d, ok), tableau%c(4), ok)
        tableau%c(2) = quad_mul(quad_rational("2/3", d, ok), tableau%c(3), ok)
        if (.not. ok) return

        s1 = quad_add(quad_add(c6, c7, ok), c8, ok)
        s2 = quad_add(quad_add(quad_mul(c6, c7, ok), quad_mul(c6, c8, ok), ok), &
                      quad_mul(c7, c8, ok), ok)
        s3 = quad_mul(quad_mul(c6, c7, ok), c8, ok)
        if (.not. ok) return

        numerator = quad_mul(quad_rational("3", d, ok), s1, ok)
        call accumulate(numerator, "-28", s2, d, ok)
        call accumulate(numerator, "189", s3, d, ok)
        call accumulate(numerator, "14", quad_mul(s1, s2, ok), d, ok)
        call accumulate(numerator, "-168", quad_mul(s1, s3, ok), d, ok)
        call accumulate(numerator, "98", quad_mul(s2, s3, ok), d, ok)

        denominator = quad_rational("6", d, ok)
        call accumulate(denominator, "-21", s1, d, ok)
        call accumulate(denominator, "35", s2, d, ok)
        call accumulate(denominator, "-42", s3, d, ok)
        call accumulate(denominator, "21", quad_pow(s1, 2, ok), d, ok)
        call accumulate(denominator, "98", quad_pow(s2, 2, ok), d, ok)
        call accumulate(denominator, "735", quad_pow(s3, 2, ok), d, ok)
        call accumulate(denominator, "-84", quad_mul(s1, s2, ok), d, ok)
        call accumulate(denominator, "168", quad_mul(s1, s3, ok), d, ok)
        call accumulate(denominator, "-490", quad_mul(s2, s3, ok), d, ok)
        if (.not. ok) return
        if (quad_is_zero(denominator)) return

        term = quad_div(numerator, denominator, ok)
        tableau%c(9) = quad_div(term, quad_rational("2", d, ok), ok)
        if (.not. ok) return
        status = RK_DP8_OK
    end subroutine build_nodes

    subroutine accumulate(total, coefficient, value, d, ok)
        type(quadratic_t), intent(inout) :: total
        character(*), intent(in) :: coefficient
        type(quadratic_t), intent(in) :: value
        integer, intent(in) :: d
        logical, intent(inout) :: ok
        type(quadratic_t) :: scaled

        if (.not. ok) return
        scaled = quad_mul(quad_rational(coefficient, d, ok), value, ok)
        if (.not. ok) return
        total = quad_add(total, scaled, ok)
    end subroutine accumulate

    !> (5.20a) with (5.20g).
    subroutine build_weights(tableau, status)
        type(rk_dp8_t), intent(inout) :: tableau
        integer, intent(out) :: status

        integer, parameter :: idx(8) = [1, 6, 7, 8, 9, 10, 11, 12]
        type(quadratic_t) :: matrix(8, 8), rhs(8)
        type(quadratic_t), allocatable :: solution(:)
        integer :: q, j, local, d
        logical :: ok
        character(len=32) :: text

        d = tableau%c(1)%radicand
        status = RK_DP8_WEIGHTS
        do q = 1, 8
            do j = 1, 8
                matrix(q, j) = quad_pow(tableau%c(idx(j)), q - 1, ok)
                if (.not. ok) return
            end do
            write (text, "(a,i0)") "1/", q
            rhs(q) = quad_rational(trim(text), d, ok)
            if (.not. ok) return
        end do
        call quad_solve(matrix, rhs, solution, local)
        if (local /= SOLVE_OK) return
        do j = 1, 8
            tableau%b(idx(j)) = solution(j)
        end do
        status = RK_DP8_OK
    end subroutine build_weights

    !> Columns of row i that (5.20h) leaves free.
    subroutine row_columns(i, columns, n)
        integer, intent(in) :: i
        integer, intent(out) :: columns(:), n
        integer :: j

        n = 0
        do j = 1, i - 1
            if (j == 2 .and. i >= 4) cycle
            if (j == 3 .and. i >= 6) cycle
            n = n + 1
            columns(n) = j
        end do
    end subroutine row_columns

    !> Rows 2..8 from (5.20b-f); the equation count follows the ranges there.
    subroutine build_early_rows(tableau, status)
        type(rk_dp8_t), intent(inout) :: tableau
        integer, intent(out) :: status

        type(quadratic_t) :: matrix(5, 8), rhs(5)
        type(quadratic_t), allocatable :: solution(:)
        integer :: i, n, used, local, columns(8), j, eq, d
        logical :: ok
        character(len=32) :: text

        d = tableau%c(1)%radicand
        status = RK_DP8_ROWS
        do i = 2, 8
            call row_columns(i, columns, n)
            used = 1
            if (i >= 3) used = 3
            if (i >= 6) used = 5
            do eq = 1, used
                do j = 1, n
                    matrix(eq, j) = quad_pow(tableau%c(columns(j)), eq - 1, ok)
                    if (.not. ok) return
                end do
                if (eq == 1) then
                    rhs(eq) = tableau%c(i)
                else
                    write (text, "(a,i0)") "1/", eq
                    rhs(eq) = quad_mul(quad_pow(tableau%c(i), eq, ok), &
                                       quad_rational(trim(text), d, ok), ok)
                    if (.not. ok) return
                end if
            end do
            call quad_solve(matrix(1:used, 1:n), rhs(1:used), solution, local)
            if (local /= SOLVE_OK) return
            do j = 1, n
                tableau%a(i, columns(j)) = solution(j)
            end do
        end do
        status = RK_DP8_OK
    end subroutine build_early_rows

    !> a(12,11) from (5.20i) at j = 11, the e vector from (5.26) and (5.27),
    !> then a(11,10) and a(12,10).
    subroutine build_tail_links(tableau, status)
        type(rk_dp8_t), intent(inout) :: tableau
        integer, intent(out) :: status

        type(quadratic_t) :: e(S), matrix(6, 6), rhs(6), two, one
        type(quadratic_t) :: m2(2, 2), r2(2)
        type(quadratic_t), allocatable :: solution(:)
        integer, parameter :: idx(6) = [1, 6, 7, 8, 9, 10]
        integer :: q, j, local, d
        logical :: ok

        d = tableau%c(1)%radicand
        status = RK_DP8_ROWS
        one = quad_rational("1", d, ok)
        two = quad_rational("2", d, ok)
        if (.not. ok) return

        tableau%a(12, 11) = quad_div( &
            quad_mul(tableau%b(11), quad_sub(one, tableau%c(11), ok), ok), &
            tableau%b(12), ok)
        if (.not. ok) return

        do j = 1, S
            e(j) = quad_rational("0", d, ok)
        end do
        e(11) = quad_sub( &
            quad_mul(quad_mul(tableau%b(12), tableau%c(12), ok), &
                     tableau%a(12, 11), ok), &
            quad_mul(quad_div(tableau%b(11), two, ok), &
                     quad_sub(one, quad_pow(tableau%c(11), 2, ok), ok), ok), ok)
        if (.not. ok) return

        do q = 1, 6
            do j = 1, 6
                matrix(q, j) = quad_pow(tableau%c(idx(j)), q - 1, ok)
                if (.not. ok) return
            end do
            rhs(q) = quad_neg(quad_mul(quad_pow(tableau%c(11), q - 1, ok), &
                                       e(11), ok), ok)
            if (.not. ok) return
        end do
        call quad_solve(matrix, rhs, solution, local)
        if (local /= SOLVE_OK) return
        do j = 1, 6
            e(idx(j)) = solution(j)
        end do

        m2(1, 1) = tableau%b(11); m2(1, 2) = tableau%b(12)
        m2(2, 1) = quad_mul(tableau%b(11), tableau%c(11), ok)
        m2(2, 2) = quad_mul(tableau%b(12), tableau%c(12), ok)
        r2(1) = quad_mul(tableau%b(10), quad_sub(one, tableau%c(10), ok), ok)
        r2(2) = quad_add(e(10), &
                         quad_mul(quad_div(tableau%b(10), two, ok), &
                                  quad_sub(one, quad_pow(tableau%c(10), 2, ok), &
                                           ok), ok), ok)
        if (.not. ok) return
        call quad_solve(m2, r2, solution, local)
        if (local /= SOLVE_OK) return
        tableau%a(11, 10) = solution(1)
        tableau%a(12, 10) = solution(2)
        status = RK_DP8_OK
    end subroutine build_tail_links

    !> Columns four and five for rows 9..12: (5.20i-k) give three equations for
    !> four unknowns, so each is a one-parameter family.
    subroutine build_families(tableau, p4, n4, p5, n5, status)
        type(rk_dp8_t), intent(in) :: tableau
        type(quadratic_t), intent(out) :: p4(4), n4(4), p5(4), n5(4)
        integer, intent(out) :: status

        status = RK_DP8_FAMILY
        call one_family(tableau, 4, p4, n4, status)
        if (status /= RK_DP8_OK) return
        call one_family(tableau, 5, p5, n5, status)
    end subroutine build_families

    subroutine one_family(tableau, column, particular, direction, status)
        type(rk_dp8_t), intent(in) :: tableau
        integer, intent(in) :: column
        type(quadratic_t), intent(out) :: particular(4), direction(4)
        integer, intent(out) :: status

        type(quadratic_t) :: matrix(3, 4), rhs(3), known
        type(quadratic_t), allocatable :: solution(:), basis(:, :)
        integer :: w, kk, i, local, d
        logical :: ok

        d = tableau%c(1)%radicand
        status = RK_DP8_FAMILY
        do w = 1, 3
            do kk = 1, 4
                matrix(w, kk) = quad_mul(tableau%b(8 + kk), &
                                         quad_pow(tableau%c(8 + kk), w - 1, ok), ok)
                if (.not. ok) return
            end do
            ! b_j vanishes for j = 4, 5 by (5.20g), so the whole right-hand
            ! side is minus what rows up to eight already contribute.
            known = quad_rational("0", d, ok)
            do i = column + 1, 8
                known = quad_add(known, &
                                 quad_mul(quad_mul(tableau%b(i), &
                                                   quad_pow(tableau%c(i), w - 1, ok), &
                                                   ok), tableau%a(i, column), ok), ok)
                if (.not. ok) return
            end do
            rhs(w) = quad_neg(known, ok)
            if (.not. ok) return
        end do
        call quad_solve(matrix, rhs, solution, local)
        if (local /= SOLVE_UNDERDETERMINED) return
        call quad_nullspace(matrix, basis, local)
        if (local /= SOLVE_UNDERDETERMINED) return
        if (size(basis, 2) /= 1) return
        particular = solution
        direction = basis(:, 1)
        status = RK_DP8_OK
    end subroutine one_family

    !> Fill rows 9..12 with the parameters still symbolic, and return the three
    !> closing conditions as exact polynomials in them.
    subroutine pending_rows(tableau, p4, n4, p5, n5, pending, r9, l4, l5, ok)
        type(rk_dp8_t), intent(in) :: tableau
        type(quadratic_t), intent(in) :: p4(4), n4(4), p5(4), n5(4)
        type(pending_t), intent(out) :: pending
        type(qmpoly_t), intent(out) :: r9, l4, l5
        logical, intent(out) :: ok

        type(qmpoly_t) :: uvar, vvar, rhs(5), term, acc
        type(quadratic_t) :: matrix(5, 8), qrhs(5)
        type(quadratic_t), allocatable :: solution(:)
        type(qmpoly_t) :: sym_rhs(5)
        integer :: i, j, kk, n, columns(8), use_columns(8), free_n, d, local
        integer :: eq
        character(len=32) :: text

        d = tableau%c(1)%radicand
        ok = .true.
        uvar = qmpoly_var(1, NVARS, d, ok); if (.not. ok) return
        vvar = qmpoly_var(2, NVARS, d, ok); if (.not. ok) return

        do i = 1, S
            do j = 1, S
                pending%a(i, j) = qmpoly_const(tableau%a(i, j), NVARS)
            end do
        end do
        do kk = 1, 4
            pending%a(8 + kk, 4) = qmpoly_add(qmpoly_const(p4(kk), NVARS), &
                                              qmpoly_mul(qmpoly_const(n4(kk), NVARS), &
                                                         uvar, ok), ok)
            pending%a(8 + kk, 5) = qmpoly_add(qmpoly_const(p5(kk), NVARS), &
                                              qmpoly_mul(qmpoly_const(n5(kk), NVARS), &
                                                         vvar, ok), ok)
            if (.not. ok) return
        end do

        r9 = qmpoly_zero(NVARS, d)
        do i = 9, S
            call row_columns(i, columns, n)
            free_n = 0
            do j = 1, n
                if (columns(j) == 4 .or. columns(j) == 5) cycle
                if (i >= 11 .and. columns(j) == 10) cycle
                if (i == 12 .and. columns(j) == 11) cycle
                free_n = free_n + 1
                use_columns(free_n) = columns(j)
            end do
            ! Build the five equations of (5.20b-f) for this row, moving the
            ! already-fixed columns to the right-hand side. Those columns carry
            ! the parameters, so the right-hand side is polynomial.
            do eq = 1, 5
                do j = 1, free_n
                    matrix(eq, j) = quad_pow(tableau%c(use_columns(j)), eq - 1, ok)
                    if (.not. ok) return
                end do
                if (eq == 1) then
                    qrhs(eq) = tableau%c(i)
                else
                    write (text, "(a,i0)") "1/", eq
                    qrhs(eq) = quad_mul(quad_pow(tableau%c(i), eq, ok), &
                                        quad_rational(trim(text), d, ok), ok)
                    if (.not. ok) return
                end if
                sym_rhs(eq) = qmpoly_const(qrhs(eq), NVARS)
                do j = 1, n
                    if (any(use_columns(1:free_n) == columns(j))) cycle
                    term = qmpoly_mul(qmpoly_const(quad_pow(tableau%c(columns(j)), &
                                                            eq - 1, ok), NVARS), &
                                      pending%a(i, columns(j)), ok)
                    if (.not. ok) return
                    sym_rhs(eq) = qmpoly_sub(sym_rhs(eq), term, ok)
                    if (.not. ok) return
                end do
            end do

            if (free_n == 4) then
                ! Row nine. Four unknowns, five equations: solve on the first
                ! four and keep the fifth as a condition on the parameters.
                call solve_symbolic(matrix(1:4, 1:4), sym_rhs(1:4), &
                                    pending, i, use_columns(1:4), ok)
                if (.not. ok) return
                acc = qmpoly_neg(sym_rhs(5), ok)
                if (.not. ok) return
                do j = 1, 4
                    term = qmpoly_mul(qmpoly_const(matrix(5, j), NVARS), &
                                      pending%a(i, use_columns(j)), ok)
                    if (.not. ok) return
                    acc = qmpoly_add(acc, term, ok)
                    if (.not. ok) return
                end do
                r9 = acc
            else
                call solve_symbolic(matrix(1:5, 1:free_n), sym_rhs(1:5), &
                                    pending, i, use_columns(1:free_n), ok)
                if (.not. ok) return
            end if
        end do

        l4 = condition_l(pending, tableau, 4, ok); if (.not. ok) return
        l5 = condition_l(pending, tableau, 5, ok)
    end subroutine pending_rows

    !> Solve a system whose matrix is in the field but whose right-hand side is
    !> polynomial, by solving once per basis direction of the right-hand side.
    !>
    !> The right-hand sides here are linear in the parameters, so writing each
    !> as constant + u*du + v*dv and solving the three field systems gives the
    !> exact polynomial solution.
    subroutine solve_symbolic(matrix, rhs, pending, row, columns, ok)
        type(quadratic_t), intent(in) :: matrix(:, :)
        type(qmpoly_t), intent(in) :: rhs(:)
        type(pending_t), intent(inout) :: pending
        integer, intent(in) :: row, columns(:)
        logical, intent(out) :: ok

        type(quadratic_t), allocatable :: constant(:), du(:), dv(:)
        type(quadratic_t), allocatable :: sol0(:), solu(:), solv(:)
        type(qmpoly_t) :: uvar, vvar, piece
        integer :: m, j, local, d
        type(qmpoly_t) :: coeff

        m = size(rhs)
        d = matrix(1, 1)%radicand
        allocate (constant(m), du(m), dv(m))
        ok = .true.
        do j = 1, m
            coeff = qmpoly_eval_var(qmpoly_eval_var(rhs(j), 1, &
                                                    quad_rational("0", d, ok), ok), &
                                    2, quad_rational("0", d, ok), ok)
            constant(j) = constant_of(coeff, d, ok); if (.not. ok) return
            coeff = qmpoly_coefficient_of(rhs(j), 1, 1, ok); if (.not. ok) return
            du(j) = constant_of(coeff, d, ok); if (.not. ok) return
            coeff = qmpoly_coefficient_of(rhs(j), 2, 1, ok); if (.not. ok) return
            dv(j) = constant_of(coeff, d, ok); if (.not. ok) return
        end do

        call quad_solve(matrix, constant, sol0, local)
        if (local /= SOLVE_OK) then
            ok = .false.
            return
        end if
        call quad_solve(matrix, du, solu, local)
        if (local /= SOLVE_OK) then
            ok = .false.
            return
        end if
        call quad_solve(matrix, dv, solv, local)
        if (local /= SOLVE_OK) then
            ok = .false.
            return
        end if

        uvar = qmpoly_var(1, NVARS, d, ok); if (.not. ok) return
        vvar = qmpoly_var(2, NVARS, d, ok); if (.not. ok) return
        do j = 1, size(columns)
            piece = qmpoly_const(sol0(j), NVARS)
            piece = qmpoly_add(piece, qmpoly_mul(qmpoly_const(solu(j), NVARS), &
                                                 uvar, ok), ok)
            if (.not. ok) return
            piece = qmpoly_add(piece, qmpoly_mul(qmpoly_const(solv(j), NVARS), &
                                                 vvar, ok), ok)
            if (.not. ok) return
            pending%a(row, columns(j)) = piece
        end do
    end subroutine solve_symbolic

    !> The value of a polynomial that must already be constant.
    function constant_of(p, d, ok) result(value)
        type(qmpoly_t), intent(in) :: p
        integer, intent(in) :: d
        logical, intent(out) :: ok
        type(quadratic_t) :: value

        value = quad_rational("0", d, ok)
        if (.not. ok) return
        if (qmpoly_is_zero(p)) return
        if (size(p%term) /= 1) then
            ok = .false.
            return
        end if
        if (any(p%term(1)%exponent /= 0)) then
            ok = .false.
            return
        end if
        value = p%term(1)%coefficient
    end function constant_of

    !> (5.20l) for a given k, as a polynomial in the parameters.
    function condition_l(pending, tableau, k, ok) result(total)
        type(pending_t), intent(in) :: pending
        type(rk_dp8_t), intent(in) :: tableau
        integer, intent(in) :: k
        logical, intent(out) :: ok
        type(qmpoly_t) :: total
        type(qmpoly_t) :: inner, term
        type(quadratic_t) :: weight
        integer :: i, j, d

        d = tableau%c(1)%radicand
        total = qmpoly_zero(NVARS, d)
        ok = .true.
        do i = k + 2, S
            inner = qmpoly_zero(NVARS, d)
            do j = k + 1, i - 1
                term = qmpoly_mul(pending%a(i, j), pending%a(j, k), ok)
                if (.not. ok) return
                inner = qmpoly_add(inner, term, ok)
                if (.not. ok) return
            end do
            weight = quad_mul(tableau%b(i), tableau%c(i), ok)
            if (.not. ok) return
            term = qmpoly_mul(qmpoly_const(weight, NVARS), inner, ok)
            if (.not. ok) return
            total = qmpoly_add(total, term, ok)
            if (.not. ok) return
        end do
    end function condition_l

    !> Eliminate one parameter with the linear row-nine condition, then take
    !> the gcd of the two remaining conditions.
    subroutine close_parameters(tableau, p4, n4, p5, n5, u, v, status)
        type(rk_dp8_t), intent(in) :: tableau
        type(quadratic_t), intent(in) :: p4(4), n4(4), p5(4), n5(4)
        type(quadratic_t), intent(out) :: u, v
        integer, intent(out) :: status

        type(pending_t) :: pending
        type(qmpoly_t) :: r9, l4, l5, vsub, c0, c1, m4, m5
        type(qpoly_t) :: poly4, poly5, common
        type(quadratic_t) :: a0, a1, zero
        logical :: ok
        integer :: d

        d = tableau%c(1)%radicand
        status = RK_DP8_CLOSURE
        zero = quad_rational("0", d, ok); if (.not. ok) return

        call pending_rows(tableau, p4, n4, p5, n5, pending, r9, l4, l5, ok)
        if (.not. ok) return

        ! Row nine's condition must be linear in v to eliminate it.
        if (qmpoly_degree(r9, 2) /= 1) return
        c1 = qmpoly_coefficient_of(r9, 2, 1, ok); if (.not. ok) return
        c0 = qmpoly_coefficient_of(r9, 2, 0, ok); if (.not. ok) return
        a1 = constant_of(c1, d, ok); if (.not. ok) return
        if (quad_is_zero(a1)) return
        ! c0 may still involve u, so v = -c0(u)/a1 is a polynomial in u.
        vsub = qmpoly_neg(c0, ok); if (.not. ok) return
        vsub = qmpoly_mul(vsub, qmpoly_const(quad_div( &
                                             quad_rational("1", d, ok), a1, ok), NVARS), ok)
        if (.not. ok) return

        m4 = qmpoly_substitute(l4, 2, vsub, ok); if (.not. ok) return
        m5 = qmpoly_substitute(l5, 2, vsub, ok); if (.not. ok) return
        poly4 = qmpoly_to_univariate(m4, 1, ok); if (.not. ok) return
        poly5 = qmpoly_to_univariate(m5, 1, ok); if (.not. ok) return

        common = qpoly_gcd(poly4, poly5, ok); if (.not. ok) return
        if (qpoly_degree(common) /= 1) return
        u = qpoly_linear_root(common, ok); if (.not. ok) return

        v = qmpoly_eval_all(vsub, [u, zero], ok)
        if (.not. ok) return
        status = RK_DP8_OK
    end subroutine close_parameters

    !> Substitute the parameters and insist every closing condition vanishes.
    subroutine apply_parameters(tableau, p4, n4, p5, n5, u, v, status)
        type(rk_dp8_t), intent(inout) :: tableau
        type(quadratic_t), intent(in) :: p4(4), n4(4), p5(4), n5(4), u, v
        integer, intent(out) :: status

        type(pending_t) :: pending
        type(qmpoly_t) :: r9, l4, l5
        type(quadratic_t) :: values(NVARS), value
        logical :: ok
        integer :: i, j

        status = RK_DP8_CLOSURE
        call pending_rows(tableau, p4, n4, p5, n5, pending, r9, l4, l5, ok)
        if (.not. ok) return
        values(1) = u
        values(2) = v
        if (.not. quad_is_zero(qmpoly_eval_all(r9, values, ok))) return
        if (.not. ok) return
        if (.not. quad_is_zero(qmpoly_eval_all(l4, values, ok))) return
        if (.not. ok) return
        if (.not. quad_is_zero(qmpoly_eval_all(l5, values, ok))) return
        if (.not. ok) return

        do i = 1, S
            do j = 1, S
                value = qmpoly_eval_all(pending%a(i, j), values, ok)
                if (.not. ok) return
                tableau%a(i, j) = value
            end do
        end do
        status = RK_DP8_OK
    end subroutine apply_parameters

end module dop853_construction
