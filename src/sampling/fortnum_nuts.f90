module fortnum_nuts
    !! No-U-Turn Hamiltonian Monte Carlo with multinomial trajectory sampling,
    !! a dense Euclidean metric and dual-averaging step-size adaptation.
    !!
    !! The algorithm follows Betancourt (2017, arXiv:1701.02434, App. A) and the
    !! reference implementation of Stan (base_nuts): trajectories double in a
    !! random direction until the generalized no-U-turn criterion fails on the
    !! whole trajectory or on any merged subtree, the new state is drawn from
    !! the trajectory with weights exp(-H) (biased progressive sampling between
    !! subtrees, uniform within a subtree), and leapfrog steps whose energy error
    !! exceeds max_energy_error end the trajectory as divergent. The step size is
    !! tuned during warm-up by the dual averaging of Hoffman and Gelman (2014,
    !! JMLR 15, 1593, Alg. 5) towards target_accept.
    !!
    !! The metric is supplied as the lower Cholesky factor L of the inverse
    !! metric (an estimate of the target covariance, e.g. a Laplace covariance).
    !! Sampling runs in z = L^-1 (x - x_ref) with the identity metric, which is
    !! equivalent to the dense metric M = (L L^T)^-1 in x. The metric is not
    !! adapted.
    !!
    !! The log density and its gradient come from a nuts_log_density callback
    !! with a class(*) context threaded unchanged. A callback that returns a non-OK status
    !! or a non-finite value marks the point as having infinite energy, so a
    !! failed model evaluation becomes a divergent transition, not an abort.
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_status, only: fortnum_status_t, status_set, FORTNUM_OK, &
        FORTNUM_DOMAIN_ERROR, FORTNUM_CONVERGENCE_ERROR
    use fortnum_rng, only: rng_t, rng_uniform, rng_normal
    implicit none
    private

    public :: nuts_options_t, nuts_result_t, nuts_sample, nuts_log_density

    abstract interface
        subroutine nuts_log_density(n, x, f, g, context, status)
            !! f = log density (up to a constant) at x, g = its gradient.
            import :: dp, fortnum_status_t
            integer, intent(in) :: n
            real(dp), intent(in) :: x(n)
            real(dp), intent(out) :: f, g(n)
            class(*), intent(inout) :: context
            type(fortnum_status_t), intent(out) :: status
        end subroutine nuts_log_density
    end interface

    type :: nuts_options_t
        integer :: n_warmup = 500
        integer :: n_samples = 500
        integer :: max_depth = 10
        real(dp) :: target_accept = 0.8_dp
        real(dp) :: step_size = 0.0_dp          !! <= 0: heuristic initial value
        real(dp) :: max_energy_error = 1000.0_dp
        logical :: adapt_step = .true.
    end type nuts_options_t

    type :: nuts_result_t
        real(dp), allocatable :: x(:, :)              !! (n, n_samples) draws after warm-up
        real(dp), allocatable :: log_density(:)       !! at each draw
        real(dp), allocatable :: accept_stat(:)       !! mean Metropolis probability
        real(dp), allocatable :: energy(:)            !! Hamiltonian at each draw
        integer, allocatable :: depth(:), n_leapfrog(:)
        logical, allocatable :: divergent(:)
        real(dp) :: step_size = 0.0_dp                !! after adaptation
        integer :: n_divergent = 0                    !! after warm-up
        integer :: n_evaluations = 0                  !! all gradient calls
    end type nuts_result_t

    ! A phase-space point in whitened coordinates.
    type :: point_t
        real(dp), allocatable :: z(:), p(:), g(:)    !! position, momentum, grad log density in z
        real(dp) :: log_density = 0.0_dp
        logical :: valid = .true.
    end type point_t

    ! Everything the recursion needs besides the trajectory state.
    type :: sampler_t
        integer :: n = 0
        real(dp), allocatable :: x_ref(:), chol(:, :)
        real(dp) :: eps = 1.0_dp, max_energy_error = 1000.0_dp
        integer :: n_leapfrog = 0, n_evaluations = 0
        real(dp) :: sum_metro = 0.0_dp
        logical :: divergent = .false.
    end type sampler_t

    real(dp), parameter :: NEG_INF = -huge(1.0_dp)

contains

    subroutine nuts_sample(f, n, x0, chol, opts, rng, context, result, status)
        !! Draw opts%n_samples states after opts%n_warmup warm-up transitions
        !! from x0. chol: lower Cholesky factor of the inverse metric (n, n).
        procedure(nuts_log_density) :: f
        integer, intent(in) :: n
        real(dp), intent(in) :: x0(n), chol(n, n)
        type(nuts_options_t), intent(in) :: opts
        type(rng_t), intent(inout) :: rng
        class(*), intent(inout) :: context
        type(nuts_result_t), intent(out) :: result
        type(fortnum_status_t), intent(out) :: status
        type(sampler_t) :: s
        type(point_t) :: current
        real(dp) :: mu, s_bar, x_bar, log_eps, eta, accept, h
        integer :: it, depth, k

        call status_set(status, FORTNUM_OK, '')
        if (n < 1 .or. opts%n_samples < 0 .or. opts%n_warmup < 0 .or. opts%max_depth < 1) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'nuts_sample: bad sizes or options')
            return
        end if
        if (opts%target_accept <= 0.0_dp .or. opts%target_accept >= 1.0_dp) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'nuts_sample: target_accept not in (0, 1)')
            return
        end if
        do k = 1, n
            if (.not. (chol(k, k) > 0.0_dp)) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, 'nuts_sample: chol needs a positive diagonal')
                return
            end if
        end do
        s%n = n
        s%x_ref = x0
        s%chol = chol
        do k = 1, n
            s%chol(1:k - 1, k) = 0.0_dp
        end do
        s%max_energy_error = opts%max_energy_error
        allocate(current%z(n), current%p(n), current%g(n))
        current%z = 0.0_dp
        current%p = 0.0_dp
        call evaluate(s, f, context, current)
        if (.not. current%valid) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, 'nuts_sample: log density not finite at x0')
            return
        end if

        s%eps = opts%step_size
        if (s%eps <= 0.0_dp) call initial_step(s, f, context, rng, current)

        ! dual averaging (Hoffman and Gelman 2014; Stan constants)
        mu = log(10.0_dp * s%eps)
        s_bar = 0.0_dp
        x_bar = 0.0_dp
        allocate(result%x(n, opts%n_samples), result%log_density(opts%n_samples), &
            result%accept_stat(opts%n_samples), result%energy(opts%n_samples), &
            result%depth(opts%n_samples), result%n_leapfrog(opts%n_samples), &
            result%divergent(opts%n_samples))
        do it = 1, opts%n_warmup + opts%n_samples
            call transition(s, f, context, rng, opts%max_depth, current, accept, depth, h)
            if (it <= opts%n_warmup) then
                if (opts%adapt_step) then
                    eta = 1.0_dp / (it + 10.0_dp)
                    s_bar = (1.0_dp - eta) * s_bar + eta * (opts%target_accept - accept)
                    log_eps = mu - s_bar * sqrt(real(it, dp)) / 0.05_dp
                    eta = real(it, dp)**(-0.75_dp)
                    x_bar = (1.0_dp - eta) * x_bar + eta * log_eps
                    s%eps = exp(log_eps)
                    if (it == opts%n_warmup) s%eps = exp(x_bar)
                end if
                cycle
            end if
            k = it - opts%n_warmup
            result%x(:, k) = to_x(s, current%z)
            result%log_density(k) = current%log_density
            result%accept_stat(k) = accept
            result%energy(k) = h
            result%depth(k) = depth
            result%n_leapfrog(k) = s%n_leapfrog
            result%divergent(k) = s%divergent
        end do
        result%step_size = s%eps
        result%n_divergent = count(result%divergent)
        result%n_evaluations = s%n_evaluations
    end subroutine nuts_sample

    subroutine transition(s, f, context, rng, max_depth, current, accept, depth, h_out)
        !! One NUTS transition from current (updated in place).
        type(sampler_t), intent(inout) :: s
        procedure(nuts_log_density) :: f
        class(*), intent(inout) :: context
        type(rng_t), intent(inout) :: rng
        integer, intent(in) :: max_depth
        type(point_t), intent(inout) :: current
        real(dp), intent(out) :: accept, h_out
        integer, intent(out) :: depth
        type(point_t) :: z_fwd, z_bck, z_sample, z_propose
        real(dp), allocatable :: rho(:), rho_fwd(:), rho_bck(:), rho_ext(:)
        real(dp), allocatable :: ps_fwd_bck(:), ps_fwd_fwd(:), ps_bck_bck(:), ps_bck_fwd(:)
        real(dp), allocatable :: p_fwd_bck(:), p_fwd_fwd(:), p_bck_bck(:), p_bck_fwd(:)
        real(dp) :: h0, lsw, lsw_sub, u
        integer :: k
        logical :: valid

        do k = 1, s%n
            call rng_normal(rng, current%p(k))
        end do
        h0 = hamiltonian(current)
        z_fwd = current
        z_bck = current
        z_sample = current
        z_propose = current
        ps_fwd_bck = current%p
        ps_fwd_fwd = current%p
        ps_bck_bck = current%p
        ps_bck_fwd = current%p
        p_fwd_bck = current%p
        p_fwd_fwd = current%p
        p_bck_bck = current%p
        p_bck_fwd = current%p
        rho = current%p
        lsw = 0.0_dp
        s%n_leapfrog = 0
        s%sum_metro = 0.0_dp
        s%divergent = .false.
        depth = 0
        do while (depth < max_depth)
            allocate(rho_fwd(s%n), rho_bck(s%n))
            rho_fwd = 0.0_dp
            rho_bck = 0.0_dp
            lsw_sub = NEG_INF
            call rng_uniform(rng, u)
            if (u > 0.5_dp) then
                rho_bck = rho
                p_bck_fwd = p_fwd_bck
                ps_bck_fwd = ps_fwd_bck
                valid = build_tree(s, f, context, rng, depth, z_fwd, z_propose, ps_fwd_bck, &
                    ps_fwd_fwd, rho_fwd, p_fwd_bck, p_fwd_fwd, h0, 1.0_dp, lsw_sub)
            else
                rho_fwd = rho
                p_fwd_bck = p_bck_fwd
                ps_fwd_bck = ps_bck_fwd
                valid = build_tree(s, f, context, rng, depth, z_bck, z_propose, ps_bck_fwd, &
                    ps_bck_bck, rho_bck, p_bck_fwd, p_bck_bck, h0, -1.0_dp, lsw_sub)
            end if
            if (.not. valid) then
                deallocate(rho_fwd, rho_bck)
                exit
            end if
            depth = depth + 1
            if (lsw_sub > lsw) then
                z_sample = z_propose
            else
                call rng_uniform(rng, u)
                if (u < exp(lsw_sub - lsw)) z_sample = z_propose
            end if
            lsw = log_sum_exp(lsw, lsw_sub)
            rho = rho_bck + rho_fwd
            valid = criterion(ps_bck_bck, ps_fwd_fwd, rho)
            rho_ext = rho_bck + p_fwd_bck
            valid = valid .and. criterion(ps_bck_bck, ps_fwd_bck, rho_ext)
            rho_ext = rho_fwd + p_bck_fwd
            valid = valid .and. criterion(ps_bck_fwd, ps_fwd_fwd, rho_ext)
            deallocate(rho_fwd, rho_bck)
            if (.not. valid) exit
        end do
        accept = 0.0_dp
        if (s%n_leapfrog > 0) accept = s%sum_metro / s%n_leapfrog
        current = z_sample
        h_out = hamiltonian(current)
    end subroutine transition

    recursive logical function build_tree(s, f, context, rng, depth, z, z_propose, ps_beg, &
            ps_end, rho, p_beg, p_end, h0, sign, lsw) result(valid)
        !! Extend z by 2**depth leapfrog steps in direction sign; z_propose is
        !! the subtree's multinomial draw, lsw accumulates log sum of weights.
        type(sampler_t), intent(inout) :: s
        procedure(nuts_log_density) :: f
        class(*), intent(inout) :: context
        type(rng_t), intent(inout) :: rng
        integer, intent(in) :: depth
        type(point_t), intent(inout) :: z, z_propose
        real(dp), intent(inout) :: ps_beg(:), ps_end(:), rho(:), p_beg(:), p_end(:), lsw
        real(dp), intent(in) :: h0, sign
        type(point_t) :: z_propose_final
        real(dp), allocatable :: ps_init_end(:), p_init_end(:), rho_init(:)
        real(dp), allocatable :: ps_final_beg(:), p_final_beg(:), rho_final(:), rho_sub(:)
        real(dp) :: h, lsw_init, lsw_final, lsw_sub, u

        if (depth == 0) then
            call leapfrog(s, f, context, z, sign * s%eps)
            s%n_leapfrog = s%n_leapfrog + 1
            h = huge(1.0_dp)
            if (z%valid) h = hamiltonian(z)
            if (.not. ieee_is_finite(h)) h = huge(1.0_dp)
            if (h - h0 > s%max_energy_error) s%divergent = .true.
            lsw = log_sum_exp(lsw, h0 - h)
            if (h0 - h > 0.0_dp) then
                s%sum_metro = s%sum_metro + 1.0_dp
            else
                s%sum_metro = s%sum_metro + exp(h0 - h)
            end if
            z_propose = z
            ps_beg = z%p
            ps_end = z%p
            rho = rho + z%p
            p_beg = z%p
            p_end = z%p
            valid = .not. s%divergent
            return
        end if

        allocate(ps_init_end(s%n), p_init_end(s%n), rho_init(s%n))
        rho_init = 0.0_dp
        lsw_init = NEG_INF
        valid = build_tree(s, f, context, rng, depth - 1, z, z_propose, ps_beg, ps_init_end, &
            rho_init, p_beg, p_init_end, h0, sign, lsw_init)
        if (.not. valid) return

        allocate(ps_final_beg(s%n), p_final_beg(s%n), rho_final(s%n))
        rho_final = 0.0_dp
        lsw_final = NEG_INF
        z_propose_final = z
        valid = build_tree(s, f, context, rng, depth - 1, z, z_propose_final, ps_final_beg, &
            ps_end, rho_final, p_final_beg, p_end, h0, sign, lsw_final)
        if (.not. valid) return

        lsw_sub = log_sum_exp(lsw_init, lsw_final)
        lsw = log_sum_exp(lsw, lsw_sub)
        if (lsw_final > lsw_sub) then
            z_propose = z_propose_final
        else
            call rng_uniform(rng, u)
            if (u < exp(lsw_final - lsw_sub)) z_propose = z_propose_final
        end if
        rho_sub = rho_init + rho_final
        rho = rho + rho_sub
        valid = criterion(ps_beg, ps_end, rho_sub)
        valid = valid .and. criterion(ps_beg, ps_final_beg, rho_init + p_final_beg)
        valid = valid .and. criterion(ps_init_end, ps_end, rho_final + p_init_end)
    end function build_tree

    subroutine leapfrog(s, f, context, z, eps)
        type(sampler_t), intent(inout) :: s
        procedure(nuts_log_density) :: f
        class(*), intent(inout) :: context
        type(point_t), intent(inout) :: z
        real(dp), intent(in) :: eps

        if (.not. z%valid) return
        z%p = z%p + 0.5_dp * eps * z%g
        z%z = z%z + eps * z%p
        call evaluate(s, f, context, z)
        if (.not. z%valid) return
        z%p = z%p + 0.5_dp * eps * z%g
    end subroutine leapfrog

    subroutine evaluate(s, f, context, z)
        !! log density and its z-gradient g_z = L^T g_x at z.
        type(sampler_t), intent(inout) :: s
        procedure(nuts_log_density) :: f
        class(*), intent(inout) :: context
        type(point_t), intent(inout) :: z
        type(fortnum_status_t) :: st
        real(dp) :: x(s%n), gx(s%n), lp

        x = to_x(s, z%z)
        call f(s%n, x, lp, gx, context, st)
        s%n_evaluations = s%n_evaluations + 1
        z%valid = st%code == FORTNUM_OK .and. ieee_is_finite(lp)
        if (z%valid) z%valid = all(ieee_is_finite(gx))
        if (.not. z%valid) then
            z%log_density = NEG_INF
            z%g = 0.0_dp
            return
        end if
        z%log_density = lp
        z%g = matmul(gx, s%chol)
    end subroutine evaluate

    function to_x(s, z) result(x)
        type(sampler_t), intent(in) :: s
        real(dp), intent(in) :: z(:)
        real(dp) :: x(s%n)
        x = s%x_ref + matmul(s%chol, z)
    end function to_x

    real(dp) function hamiltonian(z)
        type(point_t), intent(in) :: z
        if (.not. z%valid) then
            hamiltonian = huge(1.0_dp)
        else
            hamiltonian = -z%log_density + 0.5_dp * dot_product(z%p, z%p)
        end if
    end function hamiltonian

    logical function criterion(ps_minus, ps_plus, rho)
        real(dp), intent(in) :: ps_minus(:), ps_plus(:), rho(:)
        criterion = dot_product(ps_plus, rho) > 0.0_dp .and. dot_product(ps_minus, rho) > 0.0_dp
    end function criterion

    pure real(dp) function log_sum_exp(a, b)
        real(dp), intent(in) :: a, b
        if (a == NEG_INF) then
            log_sum_exp = b
        else if (b == NEG_INF) then
            log_sum_exp = a
        else if (a > b) then
            log_sum_exp = a + log1p_(exp(b - a))
        else
            log_sum_exp = b + log1p_(exp(a - b))
        end if
    end function log_sum_exp

    pure real(dp) function log1p_(x)
        real(dp), intent(in) :: x
        if (abs(x) < 1.0e-4_dp) then
            log1p_ = x * (1.0_dp - x * (0.5_dp - x / 3.0_dp))
        else
            log1p_ = log(1.0_dp + x)
        end if
    end function log1p_

    subroutine initial_step(s, f, context, rng, current)
        !! Heuristic initial step size (Hoffman and Gelman 2014, Alg. 4; Stan's
        !! init_stepsize): double or halve eps until the acceptance probability
        !! of one leapfrog step crosses 0.8.
        type(sampler_t), intent(inout) :: s
        procedure(nuts_log_density) :: f
        class(*), intent(inout) :: context
        type(rng_t), intent(inout) :: rng
        type(point_t), intent(in) :: current
        type(point_t) :: z
        real(dp) :: h0, delta_h, direction
        integer :: k, tries

        s%eps = 1.0_dp
        direction = 0.0_dp
        do tries = 1, 60
            z = current
            do k = 1, s%n
                call rng_normal(rng, z%p(k))
            end do
            h0 = hamiltonian(z)
            call leapfrog(s, f, context, z, s%eps)
            delta_h = h0 - hamiltonian(z)
            if (direction == 0.0_dp) direction = merge(1.0_dp, -1.0_dp, delta_h > log(0.8_dp))
            if (direction > 0.0_dp .and. .not. (delta_h > log(0.8_dp))) exit
            if (direction < 0.0_dp .and. delta_h > log(0.8_dp)) exit
            s%eps = merge(2.0_dp * s%eps, 0.5_dp * s%eps, direction > 0.0_dp)
            if (s%eps > 1.0e7_dp .or. s%eps < 1.0e-12_dp) exit
        end do
    end subroutine initial_step

end module fortnum_nuts
