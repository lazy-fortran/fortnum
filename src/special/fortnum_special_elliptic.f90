module fortnum_special_elliptic
    ! Complete elliptic integrals of the first and second kind,
    !
    !   K(m) = int_0^{pi/2} (1 - m sin^2 t)^(-1/2) dt,
    !   E(m) = int_0^{pi/2} (1 - m sin^2 t)^(+1/2) dt,
    !
    ! in the parameter m = k^2 on [0, 1).  Active argument: m (real scalar).
    !
    ! Algorithm (DLMF 19.2, 19.8, 19.12):
    !   m <= 1/2   Carlson duplication, K = R_F(0,1-m,1) and
    !              E = K - (m/3) R_D(0,1-m,1) (DLMF 19.25.1).
    !   m >  1/2   Logarithmic series in the complement q = 1 - m,
    !                K = -1/2 ln(q) A_K(q) + B_K(q),
    !                E = -1/2 ln(q) A_E(q) + B_E(q),
    !              with A_K(q) = (2/pi) K(q), A_E(q) = q/2 * 2F1(1/2,3/2;2;q)
    !              and B_K, B_E the DLMF 19.12.1 and 19.12.2 series.  These
    !              stay accurate to rounding as q -> 0 where K and E from
    !              m = 1 - q would lose q.
    ! The split pieces A and B are public because singular quadrature of
    ! boundary-element kernels integrates the logarithm analytically.
    !
    ! Derivatives (DLMF 19.4.1):
    !   dK/dm = (E - (1-m) K) / (2 m (1-m)),  dE/dm = (E - K) / (2 m),
    ! evaluated without cancellation near m = 0 through R_D.
    ! Derivative entry points follow the dawson convention: *_jvp (forward
    ! product) and *_vjp (scalar gradient), both with length-1 arrays.
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value, &
        ieee_positive_inf
    use fortnum_kinds, only: dp
    implicit none
    private

    public :: ellip_k, ellip_e
    public :: ellip_k_from_complement, ellip_e_from_complement
    public :: ellip_k_log_a, ellip_k_log_b
    public :: ellip_e_log_a, ellip_e_log_b
    public :: ellip_k_derivative, ellip_e_derivative
    public :: ellip_k_jvp, ellip_k_vjp
    public :: ellip_e_jvp, ellip_e_vjp

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
    real(dp), parameter :: ln4 = 1.38629436111989061883446424291635313_dp
    integer, parameter :: max_series_terms = 80
    real(dp), parameter :: series_switch = 0.5_dp

contains

    elemental function ellip_k(m) result(value)
        !! Complete elliptic integral K(m); NaN outside [0,1], +Inf at 1.
        real(dp), intent(in) :: m
        real(dp) :: value

        if (.not. (m >= 0.0_dp .and. m <= 1.0_dp)) then
            value = ieee_value(m, ieee_quiet_nan)
        else
            value = ellip_k_from_complement(1.0_dp - m)
        end if
    end function ellip_k

    elemental function ellip_e(m) result(value)
        !! Complete elliptic integral E(m); NaN outside [0,1].
        real(dp), intent(in) :: m
        real(dp) :: value

        if (.not. (m >= 0.0_dp .and. m <= 1.0_dp)) then
            value = ieee_value(m, ieee_quiet_nan)
        else
            value = ellip_e_from_complement(1.0_dp - m)
        end if
    end function ellip_e

    elemental function ellip_k_from_complement(q) result(value)
        !! K(1-q) for the complementary parameter q in [0,1], q given exactly.
        real(dp), intent(in) :: q
        real(dp) :: value

        if (.not. (q >= 0.0_dp .and. q <= 1.0_dp)) then
            value = ieee_value(q, ieee_quiet_nan)
        else if (q == 0.0_dp) then
            value = ieee_value(q, ieee_positive_inf)
        else if (q < series_switch) then
            value = -0.5_dp*log(q)*ellip_k_log_a(q) + ellip_k_log_b(q)
        else
            value = carlson_rf(0.0_dp, q, 1.0_dp)
        end if
    end function ellip_k_from_complement

    elemental function ellip_e_from_complement(q) result(value)
        !! E(1-q) for the complementary parameter q in [0,1], q given exactly.
        real(dp), intent(in) :: q
        real(dp) :: value

        if (.not. (q >= 0.0_dp .and. q <= 1.0_dp)) then
            value = ieee_value(q, ieee_quiet_nan)
        else if (q == 0.0_dp) then
            value = 1.0_dp
        else if (q < series_switch) then
            value = -0.5_dp*log(q)*ellip_e_log_a(q) + ellip_e_log_b(q)
        else
            value = carlson_rf(0.0_dp, q, 1.0_dp) &
                - (1.0_dp - q)/3.0_dp*carlson_rd(0.0_dp, q, 1.0_dp)
        end if
    end function ellip_e_from_complement

    elemental function ellip_k_log_a(q) result(value)
        !! A_K(q) in K(1-q) = -1/2 ln(q) A_K(q) + B_K(q); A_K = (2/pi) K(q).
        !! Defined for q in [0,1]; A_K(0) = 1.
        real(dp), intent(in) :: q
        real(dp) :: value

        if (.not. (q >= 0.0_dp .and. q <= 1.0_dp)) then
            value = ieee_value(q, ieee_quiet_nan)
        else if (q < series_switch) then
            value = k_series(q, .false.)
        else
            value = 2.0_dp/pi*carlson_rf(0.0_dp, 1.0_dp - q, 1.0_dp)
        end if
    end function ellip_k_log_a

    elemental function ellip_k_log_b(q) result(value)
        !! B_K(q) in K(1-q) = -1/2 ln(q) A_K(q) + B_K(q); B_K(0) = ln 4.
        real(dp), intent(in) :: q
        real(dp) :: value

        if (.not. (q >= 0.0_dp .and. q <= 1.0_dp)) then
            value = ieee_value(q, ieee_quiet_nan)
        else if (q < series_switch) then
            value = k_series(q, .true.)
        else
            value = carlson_rf(0.0_dp, q, 1.0_dp) &
                + 0.5_dp*log(q)*ellip_k_log_a(q)
        end if
    end function ellip_k_log_b

    elemental function ellip_e_log_a(q) result(value)
        !! A_E(q) in E(1-q) = -1/2 ln(q) A_E(q) + B_E(q); A_E(q) ~ q/2.
        real(dp), intent(in) :: q
        real(dp) :: value

        if (.not. (q >= 0.0_dp .and. q <= 1.0_dp)) then
            value = ieee_value(q, ieee_quiet_nan)
        else if (q < series_switch) then
            value = e_series(q, .false.)
        else
            value = e_log_a_by_rd(q)
        end if
    end function ellip_e_log_a

    elemental function ellip_e_log_b(q) result(value)
        !! B_E(q) in E(1-q) = -1/2 ln(q) A_E(q) + B_E(q); B_E(0) = 1.
        real(dp), intent(in) :: q
        real(dp) :: value

        if (.not. (q >= 0.0_dp .and. q <= 1.0_dp)) then
            value = ieee_value(q, ieee_quiet_nan)
        else if (q < series_switch) then
            value = e_series(q, .true.)
        else
            value = ellip_e_from_complement(q) &
                + 0.5_dp*log(q)*e_log_a_by_rd(q)
        end if
    end function ellip_e_log_b

    elemental function ellip_k_derivative(m) result(value)
        !! dK/dm for m in [0,1); +Inf at m = 1.
        real(dp), intent(in) :: m
        real(dp) :: value, q, k, e

        if (.not. (m >= 0.0_dp .and. m <= 1.0_dp)) then
            value = ieee_value(m, ieee_quiet_nan)
            return
        end if
        q = 1.0_dp - m
        if (q == 0.0_dp) then
            value = ieee_value(m, ieee_positive_inf)
        else if (q < series_switch) then
            k = ellip_k_from_complement(q)
            e = ellip_e_from_complement(q)
            value = (e - q*k)/(2.0_dp*m*q)
        else
            ! E - (1-m) K = m (K - R_D(0,q,1)/3); no cancellation at m -> 0.
            k = carlson_rf(0.0_dp, q, 1.0_dp)
            value = (k - carlson_rd(0.0_dp, q, 1.0_dp)/3.0_dp)/(2.0_dp*q)
        end if
    end function ellip_k_derivative

    elemental function ellip_e_derivative(m) result(value)
        !! dE/dm = (E - K)/(2m) = -R_D(0,1-m,1)/6 for m in [0,1].
        real(dp), intent(in) :: m
        real(dp) :: value

        if (.not. (m >= 0.0_dp .and. m <= 1.0_dp)) then
            value = ieee_value(m, ieee_quiet_nan)
        else if (m == 1.0_dp) then
            value = ieee_value(m, ieee_positive_inf)
            value = -value
        else
            value = -carlson_rd(0.0_dp, 1.0_dp - m, 1.0_dp)/6.0_dp
        end if
    end function ellip_e_derivative

    subroutine ellip_k_jvp(m, v, jv)
        !! Forward product jv = dK/dm * v; m, v, jv have length 1.
        real(dp), intent(in) :: m(:), v(:)
        real(dp), intent(out) :: jv(:)

        jv(1) = ellip_k_derivative(m(1))*v(1)
    end subroutine ellip_k_jvp

    subroutine ellip_k_vjp(m, u, jtu)
        !! Reverse product jtu = dK/dm * u (scalar output cotangent u).
        real(dp), intent(in) :: m(:), u(:)
        real(dp), intent(out) :: jtu(:)

        jtu(1) = u(1)*ellip_k_derivative(m(1))
    end subroutine ellip_k_vjp

    subroutine ellip_e_jvp(m, v, jv)
        !! Forward product jv = dE/dm * v; m, v, jv have length 1.
        real(dp), intent(in) :: m(:), v(:)
        real(dp), intent(out) :: jv(:)

        jv(1) = ellip_e_derivative(m(1))*v(1)
    end subroutine ellip_e_jvp

    subroutine ellip_e_vjp(m, u, jtu)
        !! Reverse product jtu = dE/dm * u (scalar output cotangent u).
        real(dp), intent(in) :: m(:), u(:)
        real(dp), intent(out) :: jtu(:)

        jtu(1) = u(1)*ellip_e_derivative(m(1))
    end subroutine ellip_e_vjp

    pure elemental function k_series(q, constant_part) result(value)
        !! Sum_n c_n q^n (1 or d_n), c_n = ((1/2)_n/n!)^2,
        !! d_n = psi(n+1) - psi(n+1/2) (DLMF 19.12.1, 19.12.3).
        real(dp), intent(in) :: q
        logical, intent(in) :: constant_part
        real(dp) :: value, term, d
        integer :: n

        term = 1.0_dp
        d = ln4
        value = merge(d, 1.0_dp, constant_part)
        do n = 0, max_series_terms
            term = term*((real(n, dp) + 0.5_dp)/real(n + 1, dp))**2*q
            d = d - 1.0_dp/(real(n + 1, dp)*real(2*n + 1, dp))
            if (constant_part) then
                value = value + term*d
            else
                value = value + term
            end if
            if (abs(term) < 1.0e-18_dp*abs(value)) exit
        end do
    end function k_series

    pure elemental function e_series(q, constant_part) result(value)
        !! DLMF 19.12.2 with e_n = (1/2)_n (3/2)_n / ((2)_n n!):
        !! A_E = q/2 sum e_n q^n,
        !! B_E = 1 + q/2 sum e_n q^n (d_n - 1/((2n+1)(2n+2))).
        real(dp), intent(in) :: q
        logical, intent(in) :: constant_part
        real(dp) :: value, term, d, x
        integer :: n

        term = 0.5_dp*q
        d = ln4
        if (constant_part) then
            value = 1.0_dp + term*(d - 0.5_dp)
        else
            value = term
        end if
        do n = 0, max_series_terms
            x = real(n, dp)
            term = term*((x + 0.5_dp)*(x + 1.5_dp)/((x + 2.0_dp)*(x + 1.0_dp)))*q
            d = d - 1.0_dp/((x + 1.0_dp)*(2.0_dp*x + 1.0_dp))
            if (constant_part) then
                value = value + term*(d - 1.0_dp/((2.0_dp*x + 3.0_dp)* &
                    (2.0_dp*x + 4.0_dp)))
            else
                value = value + term
            end if
            if (abs(term) < 1.0e-18_dp*abs(value)) exit
        end do
    end function e_series

    pure elemental function e_log_a_by_rd(q) result(value)
        !! A_E(q) = (2/pi) (K(q) - E(q)) = (2/pi) (q/3) R_D(0,1-q,1).
        real(dp), intent(in) :: q
        real(dp) :: value

        value = 2.0_dp/pi*(q/3.0_dp)*carlson_rd(0.0_dp, 1.0_dp - q, 1.0_dp)
    end function e_log_a_by_rd

    pure elemental function carlson_rf(x0, y0, z0) result(value)
        !! Carlson R_F(x,y,z) for x,y,z >= 0, at most one zero (DLMF 19.36.1).
        real(dp), intent(in) :: x0, y0, z0
        real(dp) :: value, x, y, z, sx, sy, sz, lambda
        real(dp) :: ave, dx, dy, dz, e2, e3
        real(dp), parameter :: errtol = 1.0e-3_dp
        integer :: iter

        x = x0
        y = y0
        z = z0
        do iter = 1, 100
            sx = sqrt(x)
            sy = sqrt(y)
            sz = sqrt(z)
            lambda = sx*(sy + sz) + sy*sz
            x = 0.25_dp*(x + lambda)
            y = 0.25_dp*(y + lambda)
            z = 0.25_dp*(z + lambda)
            ave = (x + y + z)/3.0_dp
            dx = (ave - x)/ave
            dy = (ave - y)/ave
            dz = (ave - z)/ave
            if (max(abs(dx), abs(dy), abs(dz)) < errtol) exit
        end do
        e2 = dx*dy - dz*dz
        e3 = dx*dy*dz
        value = (1.0_dp + (e2/24.0_dp - 0.1_dp - 3.0_dp*e3/44.0_dp)*e2 &
            + e3/14.0_dp)/sqrt(ave)
    end function carlson_rf

    pure elemental function carlson_rd(x0, y0, z0) result(value)
        !! Carlson R_D(x,y,z) for x,y >= 0 (not both 0), z > 0 (DLMF 19.36.1).
        real(dp), intent(in) :: x0, y0, z0
        real(dp) :: value, x, y, z, sx, sy, sz, lambda, total, fac
        real(dp) :: ave, dx, dy, dz, ea, eb, ec, ed, ee
        real(dp), parameter :: errtol = 5.0e-4_dp
        real(dp), parameter :: c1 = 3.0_dp/14.0_dp, c2 = 1.0_dp/6.0_dp
        real(dp), parameter :: c3 = 9.0_dp/22.0_dp, c4 = 3.0_dp/26.0_dp
        real(dp), parameter :: c5 = 0.25_dp*c3, c6 = 1.5_dp*c4
        integer :: iter

        x = x0
        y = y0
        z = z0
        total = 0.0_dp
        fac = 1.0_dp
        do iter = 1, 100
            sx = sqrt(x)
            sy = sqrt(y)
            sz = sqrt(z)
            lambda = sx*(sy + sz) + sy*sz
            total = total + fac/(sz*(z + lambda))
            fac = 0.25_dp*fac
            x = 0.25_dp*(x + lambda)
            y = 0.25_dp*(y + lambda)
            z = 0.25_dp*(z + lambda)
            ave = 0.2_dp*(x + y + 3.0_dp*z)
            dx = (ave - x)/ave
            dy = (ave - y)/ave
            dz = (ave - z)/ave
            if (max(abs(dx), abs(dy), abs(dz)) < errtol) exit
        end do
        ea = dx*dy
        eb = dz*dz
        ec = ea - eb
        ed = ea - 6.0_dp*eb
        ee = ed + ec + ec
        value = 3.0_dp*total + fac*(1.0_dp + ed*(-c1 + c5*ed - c6*dz*ee) &
            + dz*(c2*ee + dz*(-c3*ec + dz*c4*ea)))/(ave*sqrt(ave))
    end function carlson_rd

end module fortnum_special_elliptic
