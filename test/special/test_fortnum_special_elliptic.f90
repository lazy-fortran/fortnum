program test_fortnum_special_elliptic
    ! Oracles: mpmath (50 digits) tables in the complementary parameter q = 1-m,
    ! closed values (DLMF 19.6, Gamma(1/4) lemniscatic value), Legendre's
    ! relation (DLMF 19.7.1), the q -> 0 asymptotic series (DLMF 19.12),
    ! eighth-order finite differences for the derivatives, and the closed
    ! forms of the toroidal functions Q_{-1/2}, P_{-1/2} (DLMF 14.5.15-16)
    ! against fortnum_special_toroidal.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special_elliptic
    use fortnum_special_toroidal, only: toroidal_p, toroidal_q
    implicit none

    real(dp), parameter :: pi = 3.14159265358979323846264338327950288_dp
    integer, parameter :: nq = 11
    real(dp), parameter :: qs(nq) = [ &
        1e-14_dp, &
        1e-10_dp, &
        1e-6_dp, &
        1e-3_dp, &
        0.05_dp, &
        0.3_dp, &
        0.49_dp, &
        0.51_dp, &
        0.7_dp, &
        0.9_dp, &
        0.999984741210937500_dp ]
    real(dp), parameter :: ref_kv(nq) = [ &
        17.5043900120782516679_dp, &
        12.8992198263875995346_dp, &
        8.29405146361543998532_dp, &
        4.84113256055029703034_dp, &
        2.90833724844455210009_dp, &
        2.07536313529246914385_dp, &
        1.86264080233273855203_dp, &
        1.84569399837472351759_dp, &
        1.71388944817879106204_dp, &
        1.61244134872021939823_dp, &
        1.57080231895878055623_dp ]
    real(dp), parameter :: ref_ev(nq) = [ &
        1.00000000000008502195_dp, &
        1.00000000061996099133_dp, &
        1.00000389702617206115_dp, &
        1.00217079083444516587_dp, &
        1.06047372776627824271_dp, &
        1.24167056794582275087_dp, &
        1.34559224536798230385_dp, &
        1.35566113557195546431_dp, &
        1.44536306441266526201_dp, &
        1.53075763689776320247_dp, &
        1.57079033466530026071_dp ]
    real(dp), parameter :: ref_ka(nq) = [ &
        1.0000000000000025_dp, &
        1.000000000025_dp, &
        1.0000002500001406251_dp, &
        1.00025014072273107868_dp, &
        1.01286425658832983978_dp, &
        1.0910959103627815664_dp, &
        1.17500529310552753516_dp, &
        1.18579396358363709292_dp, &
        1.32121720676996162832_dp, &
        1.64126441434237073329_dp, &
        4.41272640782190234782_dp ]
    real(dp), parameter :: ref_kb(nq) = [ &
        1.38629436111989158457_dp, &
        1.38629436112954797786_dp, &
        1.38629445769351178397_dp, &
        1.38639096561026149641_dp, &
        1.39120217734908697413_dp, &
        1.41853823379845003969_dp, &
        1.44354585528661071141_dp, &
        1.44647004503867649336_dp, &
        1.47826691160100873731_dp, &
        1.52597911620724316533_dp, &
        1.57076865227119972201_dp ]
    real(dp), parameter :: ref_ea(nq) = [ &
        5.00000000000001875e-15_dp, &
        5.00000000018750000001e-11_dp, &
        5.00000187500117187585e-7_dp, &
        0.000500187617273016565572_dp, &
        0.025483954429660884736_dp, &
        0.170949205307880796933_dp, &
        0.311964609570132416903_dp, &
        0.329163334637889537754_dp, &
        0.530745172448766505282_dp, &
        0.937942975490848048917_dp, &
        3.77607539748266819214_dp ]
    real(dp), parameter :: ref_eb(nq) = [ &
        1.00000000000000443147_dp, &
        1.00000000004431471806_dp, &
        1.00000044314723736516_dp, &
        1.0004432040075955774_dp, &
        1.02230217539492109167_dp, &
        1.13876147088991303822_dp, &
        1.23432228573868678947_dp, &
        1.24484096631567490885_dp, &
        1.350711312098206941_dp, &
        1.48134655912008735919_dp, &
        1.57076152527651436654_dp ]
    real(dp), parameter :: ref_dk(nq) = [ &
        49999999999995.9989025_dp, &
        4999999997.15019504315_dp, &
        499998.301485655708512_dp, &
        499.163993130077511932_dp, &
        9.63217752993737513379_dp, &
        1.47395625561448097075_dp, &
        0.866142961634534640573_dp, &
        0.829046011206175410855_dp, &
        0.58485821592264647282_dp, &
        0.44200235027536524479_dp, &
        0.392705822932382407164_dp ]
    real(dp), parameter :: ref_de(nq) = [ &
        -8.25219500603916584494_dp, &
        -5.94960991347878026298_dp, &
        -3.64702743032206428415_dp, &
        -1.92140228714507100324_dp, &
        -0.972559747725407293358_dp, &
        -0.595494690961890280701_dp, &
        -0.506910349965447302134_dp, &
        -0.500033533472212299257_dp, &
        -0.447543972943543000046_dp, &
        -0.408418559112280978804_dp, &
        -0.392701328762323611648_dp ]

    integer :: nfail, npass

    nfail = 0
    npass = 0
    call test_closed_values()
    call test_table()
    call test_legendre_relation()
    call test_asymptotics()
    call test_derivatives()
    call test_ad_products()
    call test_toroidal()
    call test_domain()

    write (*, "(a,i0,a,i0,a)") "elliptic: ", npass, " passed, ", nfail, " failed"
    if (nfail > 0) stop 1

contains

    subroutine check(label, got, ref, tol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, ref, tol
        real(dp) :: err

        err = abs(got - ref)/max(abs(ref), tiny(1.0_dp))
        if (err <= tol) then
            npass = npass + 1
        else
            nfail = nfail + 1
            write (error_unit, "(a,a,3es24.15)") "FAIL ", label, got, ref, err
        end if
    end subroutine check

    subroutine check_small(label, got, tol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, tol

        if (abs(got) <= tol) then
            npass = npass + 1
        else
            nfail = nfail + 1
            write (error_unit, "(a,a,2es24.15)") "FAIL ", label, got, tol
        end if
    end subroutine check_small

    subroutine report(label, err)
        character(*), intent(in) :: label
        real(dp), intent(in) :: err

        write (*, "(a,es10.3)") "max rel err "//label//": ", err
    end subroutine report

    subroutine test_closed_values()
        real(dp) :: lemn

        call check("K(0)", ellip_k(0.0_dp), 0.5_dp*pi, 4.0e-16_dp)
        call check("E(0)", ellip_e(0.0_dp), 0.5_dp*pi, 4.0e-16_dp)
        call check("E(1)", ellip_e(1.0_dp), 1.0_dp, 0.0_dp)
        if (ellip_k(1.0_dp) > huge(1.0_dp)) then
            npass = npass + 1
        else
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL K(1) is not +Inf"
        end if
        ! K(1/2) = Gamma(1/4)^2 / (4 sqrt(pi)), E(1/2) via Legendre at m=1/2:
        ! 2 E K - K^2 = pi/2.
        lemn = gamma(0.25_dp)**2/(4.0_dp*sqrt(pi))
        call check("K(1/2)", ellip_k(0.5_dp), lemn, 4.0e-16_dp)
        call check("K(1/2) lit", ellip_k(0.5_dp), &
            1.8540746773013719184338503472_dp, 2.0e-16_dp)
        call check("E(1/2)", ellip_e(0.5_dp), &
            1.35064388104767550252017473534_dp, 2.0e-16_dp)
        ! A&S 17.1 table values (K and E at m = 0.1 and 0.9).
        call check("K(.1)", ellip_k(0.1_dp), 1.6124413487202193982_dp, 4.0e-16_dp)
        call check("E(.1)", ellip_e(0.1_dp), 1.5307576368977632025_dp, 4.0e-16_dp)
        call check("K(.9)", ellip_k(0.9_dp), 2.5780921133481732927_dp, 4.0e-16_dp)
        call check("E(.9)", ellip_e(0.9_dp), 1.1047747327040733079_dp, 4.0e-16_dp)
    end subroutine test_closed_values

    subroutine test_table()
        real(dp) :: ek, ee, eka, ekb, eea, eeb, edk, ede, r
        integer :: i

        ek = 0.0_dp; ee = 0.0_dp; eka = 0.0_dp; ekb = 0.0_dp
        eea = 0.0_dp; eeb = 0.0_dp; edk = 0.0_dp; ede = 0.0_dp
        do i = 1, nq
            r = abs(ellip_k_from_complement(qs(i)) - ref_kv(i))/ref_kv(i)
            ek = max(ek, r)
            r = abs(ellip_e_from_complement(qs(i)) - ref_ev(i))/ref_ev(i)
            ee = max(ee, r)
            r = abs(ellip_k_log_a(qs(i)) - ref_ka(i))/ref_ka(i)
            eka = max(eka, r)
            r = abs(ellip_k_log_b(qs(i)) - ref_kb(i))/ref_kb(i)
            ekb = max(ekb, r)
            r = abs(ellip_e_log_a(qs(i)) - ref_ea(i))/ref_ea(i)
            eea = max(eea, r)
            r = abs(ellip_e_log_b(qs(i)) - ref_eb(i))/ref_eb(i)
            eeb = max(eeb, r)
            ! m = 1 - q is rounded; the derivatives (~1/q) need q >= 1e-6.
            if (qs(i) < 1.0e-6_dp) cycle
            r = abs(ellip_k_derivative(1.0_dp - qs(i)) - ref_dk(i))/ref_dk(i)
            edk = max(edk, r)
            r = abs(ellip_e_derivative(1.0_dp - qs(i)) - ref_de(i))/abs(ref_de(i))
            ede = max(ede, r)
        end do
        call report("K", ek)
        call report("E", ee)
        call report("A_K", eka)
        call report("B_K", ekb)
        call report("A_E", eea)
        call report("B_E", eeb)
        call report("dK/dm", edk)
        call report("dE/dm", ede)
        call check_small("table K", ek, 1.0e-15_dp)
        call check_small("table E", ee, 1.0e-15_dp)
        call check_small("table A_K", eka, 1.0e-15_dp)
        call check_small("table B_K", ekb, 1.0e-15_dp)
        call check_small("table A_E", eea, 1.0e-15_dp)
        call check_small("table B_E", eeb, 1.0e-15_dp)
        call check_small("table dK", edk, 1.0e-9_dp)
        call check_small("table dE", ede, 1.0e-9_dp)
    end subroutine test_table

    subroutine test_legendre_relation()
        ! E K' + E' K - K K' = pi/2 with K' = K(1-m), using exact complements.
        integer, parameter :: n = 400
        real(dp) :: m, q, k, e, kp, ep, worst
        integer :: i

        worst = 0.0_dp
        do i = 1, n - 1
            m = real(i, dp)/real(n, dp)
            q = real(n - i, dp)/real(n, dp)
            k = ellip_k_from_complement(q)
            e = ellip_e_from_complement(q)
            kp = ellip_k(q)
            ep = ellip_e(q)
            worst = max(worst, abs(e*kp + ep*k - k*kp - 0.5_dp*pi)/(0.5_dp*pi))
        end do
        ! Strongly singular end with q = 10^-j given exactly.
        do i = 1, 14
            q = 10.0_dp**(-i)
            k = ellip_k_from_complement(q)
            e = ellip_e_from_complement(q)
            kp = ellip_k(q)
            ep = ellip_e(q)
            worst = max(worst, abs(e*kp + ep*k - k*kp - 0.5_dp*pi)/(0.5_dp*pi))
        end do
        call report("Legendre relation", worst)
        call check_small("Legendre relation", worst, 1.0e-14_dp)
    end subroutine test_legendre_relation

    subroutine test_asymptotics()
        ! DLMF 19.12: K = L + q/4 (L-1) + 9q^2/64 (L-7/6) + ...,
        ! E = 1 + q/2 (L-1/2) + 3q^2/16 (L-13/12) + ..., L = ln(4/sqrt(q)).
        real(dp) :: q, l, k_as, e_as, worst
        integer :: i

        worst = 0.0_dp
        do i = 1, 3
            q = 10.0_dp**(-3 - 2*i)
            l = log(4.0_dp/sqrt(q))
            k_as = l + 0.25_dp*q*(l - 1.0_dp) &
                + 9.0_dp/64.0_dp*q*q*(l - 7.0_dp/6.0_dp)
            e_as = 1.0_dp + 0.5_dp*q*(l - 0.5_dp) &
                + 3.0_dp/16.0_dp*q*q*(l - 13.0_dp/12.0_dp)
            call check("K asymptotic", ellip_k_from_complement(q), k_as, 1.0e-15_dp)
            call check("E asymptotic", ellip_e_from_complement(q), e_as, 1.0e-15_dp)
            worst = max(worst, abs(ellip_k_log_a(q) - 1.0_dp - 0.25_dp*q &
                - 9.0_dp/64.0_dp*q*q))
            worst = max(worst, abs(ellip_k_log_b(q) - log(4.0_dp) &
                - 0.25_dp*q*(log(4.0_dp) - 1.0_dp) &
                - 9.0_dp/64.0_dp*q*q*(log(4.0_dp) - 7.0_dp/6.0_dp)))
        end do
        call check_small("split constants", worst, 2.0e-15_dp)
        call check("A_K(0)", ellip_k_log_a(0.0_dp), 1.0_dp, 0.0_dp)
        call check("B_K(0)", ellip_k_log_b(0.0_dp), log(4.0_dp), 2.0e-16_dp)
        if (ellip_e_log_a(0.0_dp) /= 0.0_dp) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL A_E(0) /= 0"
        else
            npass = npass + 1
        end if
        call check("B_E(0)", ellip_e_log_b(0.0_dp), 1.0_dp, 0.0_dp)
    end subroutine test_asymptotics

    subroutine test_derivatives()
        ! Eighth-order central differences, relative step h = 0.02 min(m,1-m).
        real(dp), parameter :: c(4) = [4.0_dp/5.0_dp, -1.0_dp/5.0_dp, &
            4.0_dp/105.0_dp, -1.0_dp/280.0_dp]
        real(dp) :: m, h, fdk, fde, worst_k, worst_e
        integer :: i, j

        worst_k = 0.0_dp
        worst_e = 0.0_dp
        do i = 1, 19
            m = 0.05_dp*real(i, dp)
            h = 0.02_dp*min(m, 1.0_dp - m)
            fdk = 0.0_dp
            fde = 0.0_dp
            do j = 1, 4
                fdk = fdk + c(j)*(ellip_k(m + j*h) - ellip_k(m - j*h))
                fde = fde + c(j)*(ellip_e(m + j*h) - ellip_e(m - j*h))
            end do
            fdk = fdk/h
            fde = fde/h
            worst_k = max(worst_k, abs(fdk - ellip_k_derivative(m)) &
                /abs(ellip_k_derivative(m)))
            worst_e = max(worst_e, abs(fde - ellip_e_derivative(m)) &
                /abs(ellip_e_derivative(m)))
        end do
        call report("dK/dm vs FD", worst_k)
        call report("dE/dm vs FD", worst_e)
        call check_small("dK/dm FD", worst_k, 1.0e-9_dp)
        call check_small("dE/dm FD", worst_e, 1.0e-9_dp)
        call check("dK/dm(0)", ellip_k_derivative(0.0_dp), pi/8.0_dp, 4.0e-16_dp)
        call check("dE/dm(0)", ellip_e_derivative(0.0_dp), -pi/8.0_dp, 4.0e-16_dp)
        ! Closed forms of DLMF 19.4.1 from the primal values.
        m = 0.3_dp
        call check("dK closed form", ellip_k_derivative(m), &
            (ellip_e(m) - (1.0_dp - m)*ellip_k(m))/(2.0_dp*m*(1.0_dp - m)), &
            4.0e-15_dp)
        call check("dE closed form", ellip_e_derivative(m), &
            (ellip_e(m) - ellip_k(m))/(2.0_dp*m), 4.0e-15_dp)
    end subroutine test_derivatives

    subroutine test_ad_products()
        real(dp) :: jv(1), jtu(1)

        call ellip_k_jvp([0.3_dp], [2.0_dp], jv)
        call check("K jvp", jv(1), 2.0_dp*ellip_k_derivative(0.3_dp), 0.0_dp)
        call ellip_k_vjp([0.3_dp], [-3.0_dp], jtu)
        call check("K vjp", jtu(1), -3.0_dp*ellip_k_derivative(0.3_dp), 0.0_dp)
        call ellip_e_jvp([0.8_dp], [2.0_dp], jv)
        call check("E jvp", jv(1), 2.0_dp*ellip_e_derivative(0.8_dp), 0.0_dp)
        call ellip_e_vjp([0.8_dp], [-3.0_dp], jtu)
        call check("E vjp", jtu(1), -3.0_dp*ellip_e_derivative(0.8_dp), 0.0_dp)
    end subroutine test_ad_products

    subroutine test_toroidal()
        ! Q_{-1/2}(x) = sqrt(2/(x+1)) K(2/(x+1)),
        ! P_{-1/2}(x) = (2/pi) sqrt(2/(x+1)) K((x-1)/(x+1)).
        real(dp), parameter :: xs(5) = [1.05_dp, 1.5_dp, 2.0_dp, 5.0_dp, 20.0_dp]
        real(dp) :: s, worst
        integer :: i

        worst = 0.0_dp
        do i = 1, size(xs)
            s = sqrt(2.0_dp/(xs(i) + 1.0_dp))
            worst = max(worst, abs(toroidal_q(0, 0, xs(i)) &
                - s*ellip_k(2.0_dp/(xs(i) + 1.0_dp)))/toroidal_q(0, 0, xs(i)))
            worst = max(worst, abs(toroidal_p(0, 0, xs(i)) &
                - 2.0_dp/pi*s*ellip_k((xs(i) - 1.0_dp)/(xs(i) + 1.0_dp))) &
                /toroidal_p(0, 0, xs(i)))
        end do
        call report("toroidal cross-check", worst)
        call check_small("toroidal cross-check", worst, 1.0e-10_dp)
    end subroutine test_toroidal

    subroutine test_domain()
        if (ieee_is_nan(ellip_k(-0.1_dp)) .and. ieee_is_nan(ellip_e(1.1_dp)) &
            .and. ieee_is_nan(ellip_k_from_complement(-1.0e-3_dp))) then
            npass = npass + 1
        else
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL domain NaN"
        end if
    end subroutine test_domain

end program test_fortnum_special_elliptic
