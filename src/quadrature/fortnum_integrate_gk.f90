module fortnum_integrate_gk
    ! Adaptive Gauss-Kronrod quadrature (QAG pattern).
    !
    ! Derivative policy: trace_rule (ad.md §1, §4).
    !   integrate_gk is adaptive: the subdivision schedule is data-dependent.
    !   Differentiate at the frozen subdivision the primal chose; the
    !   derivative of the integral w.r.t. integrand values is a fixed weighted
    !   sum once the schedule is frozen (ad.md §4, integrate row).
    !   Active arguments: f (integrand), a, b.
    !   Inactive: epsabs, epsrel, key, limit, abserr, ierr (control/status).
    !
    ! Algorithm: Piessens, de Doncker-Kapenga, Ueberhuber, Kahaner,
    !   "QUADPACK" (Springer 1983). Nodes and weights regenerated to 31
    !   significant digits with mpmath via the Stieltjes polynomial of the
    !   Legendre weight (Kronrod 1965, DLMF 3.5(v)); they agree with the
    !   public-domain Netlib QUADPACK dqk15/21/31/61 tables.
    !
    ! Keys: 15, 21, 31, 61 select the Gauss-Kronrod pair (G7K15, G10K21,
    !   G15K31, G30K61). Default is 21.

    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    ! gk_apply is the single-panel GK pair; the global adaptive driver in
    ! fortnum_integrate reuses it per subinterval instead of re-deriving the
    ! rule or recursing through integrate_gk's own adaptation (ADR integrate.md
    ! section 1).
    public :: integrate_gk, gk_integrand_t, gk_apply

    abstract interface
        function gk_integrand_t(x) result(fx)
            import :: dp
            real(dp), intent(in) :: x
            real(dp) :: fx
        end function gk_integrand_t
    end interface

    ! G7K15 nodes (positive half, outermost first) and weights.
    real(dp), parameter :: xgk15(8) = [ &
        0.9914553711208126392068546975263_dp, 0.9491079123427585245261896840479_dp, &
        0.8648644233597690727897127886409_dp, 0.7415311855993944398638647732808_dp, &
        0.5860872354676911302941448382587_dp, 0.4058451513773971669066064120770_dp, &
        0.2077849550078984676006894037732_dp, 0.0_dp]
    real(dp), parameter :: wgk15(8) = [ &
        0.02293532201052922496373200805897_dp, 0.06309209262997855329070066318920_dp, &
        0.1047900103222501838398763225415_dp,  0.1406532597155259187451895905102_dp,  &
        0.1690047266392679028265834265986_dp,  0.1903505780647854099132564024210_dp,  &
        0.2044329400752988924141619992346_dp,  0.2094821410847278280129991748917_dp]
    real(dp), parameter :: wg15(4) = [ &
        0.1294849661688696932706114326791_dp, 0.2797053914892766679014677714238_dp, &
        0.3818300505051189449503697754890_dp, 0.4179591836734693877551020408163_dp]

    ! G10K21
    real(dp), parameter :: xgk21(11) = [ &
        0.9956571630258080807355272806890_dp, 0.9739065285171717200779640120845_dp, &
        0.9301574913557082260012071800595_dp, 0.8650633666889845107320966884235_dp, &
        0.7808177265864168970637175783450_dp, 0.6794095682990244062343273651149_dp, &
        0.5627571346686046833390000992727_dp, 0.4333953941292471907992659431658_dp, &
        0.2943928627014601981311266031039_dp, 0.1488743389816312108848260011297_dp, &
        0.0_dp]
    real(dp), parameter :: wgk21(11) = [ &
        0.01169463886737187427806439606219_dp, 0.03255816230796472747881897245939_dp, &
        0.05475589657435199603138130024458_dp, 0.07503967481091995276704314091619_dp, &
        0.09312545458369760553506546508337_dp, 0.1093871588022976418992105903258_dp,  &
        0.1234919762620658510779581098311_dp,  0.1347092173114733259280540017717_dp,  &
        0.1427759385770600807970942731387_dp,  0.1477391049013384913748415159721_dp,  &
        0.1494455540029169056649364683898_dp]
    real(dp), parameter :: wg21(5) = [ &
        0.06667134430868813759356880989333_dp, 0.1494513491505805931457763396577_dp, &
        0.2190863625159820439955349342282_dp,  0.2692667193099963550912269215695_dp,  &
        0.2955242247147528701738929946513_dp]

    ! G15K31
    real(dp), parameter :: xgk31(16) = [ &
        0.9980022986933970602851728401523_dp, 0.9879925180204854284895657185866_dp, &
        0.9677390756791391342573479787843_dp, 0.9372733924007059043077589477102_dp, &
        0.8972645323440819008825096564545_dp, 0.8482065834104272162006483207742_dp, &
        0.7904185014424659329676492948179_dp, 0.7244177313601700474161860546139_dp, &
        0.6509967412974169705337358953133_dp, 0.5709721726085388475372267372539_dp, &
        0.4850818636402396806936557402324_dp, 0.3941513470775633698972073709810_dp, &
        0.2991800071531688121667800242664_dp, 0.2011940939974345223006283033946_dp, &
        0.1011420669187174990270742314474_dp, 0.0_dp]
    real(dp), parameter :: wgk31(16) = [ &
        0.005377479872923348987792051430128_dp, 0.01500794732931612253837476307581_dp, &
        0.02546084732671532018687400101965_dp,  0.03534636079137584622203794847836_dp, &
        0.04458975132476487660822729937328_dp,  0.05348152469092808726534314723943_dp, &
        0.06200956780067064028513923096080_dp,  0.06985412131872825870952007709915_dp, &
        0.07684968075772037889443277748266_dp,  0.08308050282313302103828924728610_dp, &
        0.08856444305621177064727544369377_dp,  0.09312659817082532122548687274735_dp, &
        0.09664272698362367850517990762759_dp,  0.09917359872179195933239317348460_dp, &
        0.1007698455238755950449466626176_dp,   0.1013300070147915490173747927675_dp]
    real(dp), parameter :: wg31(8) = [ &
        0.03075324199611726835462839357720_dp, 0.07036604748810812470926741645067_dp, &
        0.1071592204671719350118695466859_dp,  0.1395706779261543144478047945110_dp,  &
        0.1662692058169939335532008604812_dp,  0.1861610000155622110268005618664_dp,  &
        0.1984314853271115764561183264438_dp,  0.2025782419255612728806201999675_dp]

    ! G30K61
    real(dp), parameter :: xgk61(31) = [ &
        0.9994844100504906375713258957058_dp, 0.9968934840746495402716300509187_dp, &
        0.9916309968704045948586283661095_dp, 0.9836681232797472099700325816057_dp, &
        0.9731163225011262683746938684237_dp, 0.9600218649683075122168710255818_dp, &
        0.9443744447485599794158313240374_dp, 0.9262000474292743258793242770805_dp, &
        0.9055733076999077985465225589260_dp, 0.8825605357920526815431164625302_dp, &
        0.8572052335460610989586585106589_dp, 0.8295657623827683974428981197325_dp, &
        0.7997278358218390830136689423227_dp, 0.7677774321048261949179773409745_dp, &
        0.7337900624532268047261711313695_dp, 0.6978504947933157969322923880266_dp, &
        0.6600610641266269613700536681493_dp, 0.6205261829892428611404775564312_dp, &
        0.5793452358263616917560249321725_dp, 0.5366241481420198992641697933111_dp, &
        0.4924804678617785749936930612077_dp, 0.4470337695380891767806099003229_dp, &
        0.4004012548303943925354762115427_dp, 0.3527047255308781134710372070894_dp, &
        0.3040732022736250773726771071993_dp, 0.2546369261678898464398051298178_dp, &
        0.2045251166823098914389576710020_dp, 0.1538699136085835469637946727433_dp, &
        0.1028069379667370301470967513180_dp, 0.05147184255531769583302521316672_dp, &
        0.0_dp]
    real(dp), parameter :: wgk61(31) = [ &
        0.001389013698677007624551591226760_dp, 0.003890461127099884051267201844516_dp, &
        0.006630703915931292173319826369750_dp, 0.009273279659517763428441146892024_dp, &
        0.01182301525349634174223289885325_dp,  0.01436972950704580481245143244358_dp,  &
        0.01692088918905327262757228942032_dp,  0.01941414119394238117340895105013_dp,  &
        0.02182803582160919229716748573834_dp,  0.02419116207808060136568637072523_dp,  &
        0.02650995488233310161060170933508_dp,  0.02875404876504129284397878535433_dp,  &
        0.03090725756238776247288425294309_dp,  0.03298144705748372603181419101685_dp,  &
        0.03497933802806002413749967073147_dp,  0.03688236465182122922391106561714_dp,  &
        0.03867894562472759295034865153228_dp,  0.04037453895153595911199527975247_dp,  &
        0.04196981021516424614714754128597_dp,  0.04345253970135606931683172811707_dp,  &
        0.04481480013316266319235555161672_dp,  0.04605923827100698811627173555937_dp,  &
        0.04718554656929915394526147818110_dp,  0.04818586175708712914077949229830_dp,  &
        0.04905543455502977888752816536724_dp,  0.04979568342707420635781156937994_dp,  &
        0.05040592140278234684089308565359_dp,  0.05088179589874960649229747304980_dp,  &
        0.05122154784925877217065628260494_dp,  0.05142612853745902593386287921578_dp,  &
        0.05149472942945156755834043364710_dp]
    real(dp), parameter :: wg61(15) = [ &
        0.007968192496166605615465883474674_dp, 0.01846646831109095914230213191205_dp, &
        0.02878470788332336934971917961129_dp,  0.03879919256962704959680193644635_dp, &
        0.04840267283059405290293814042281_dp,  0.05749315621761906648172168940206_dp, &
        0.06597422988218049512812851511596_dp,  0.07375597473770520626824385002219_dp, &
        0.08075589522942021535469493846053_dp,  0.08689978720108297980238753071513_dp, &
        0.09212252223778612871763270708762_dp,  0.09636873717464425963946862635181_dp, &
        0.09959342058679526706278028210357_dp,  0.1017623897484055045964289521686_dp,  &
        0.1028526528935588403412856367054_dp]

    real(dp), parameter :: epmach = epsilon(1.0_dp)
    real(dp), parameter :: uflow  = tiny(1.0_dp)

contains

    ! Integrate f over [a,b] to tolerances epsabs+epsrel*|result|.
    ! key selects the GK pair (15/21/31/61); default 21.
    ! limit bounds adaptive subdivisions; default 200.
    ! ierr == 0: converged. 1: max subdivisions reached. 2: roundoff detected.
    ! 3: bad integrand behaviour. 6: invalid input.
    subroutine integrate_gk(f, a, b, epsabs, epsrel, result, abserr, ierr, &
                            key, limit)
        procedure(gk_integrand_t) :: f
        real(dp), intent(in)  :: a, b, epsabs, epsrel
        real(dp), intent(out) :: result, abserr
        integer,  intent(out) :: ierr
        integer,  intent(in), optional :: key, limit

        integer :: key_loc, limit_loc

        key_loc   = 21
        limit_loc = 200
        if (present(key))   key_loc   = key
        if (present(limit)) limit_loc = limit

        result = 0.0_dp
        abserr = 0.0_dp
        if (limit_loc < 1 .or. &
            (key_loc /= 15 .and. key_loc /= 21 .and. &
             key_loc /= 31 .and. key_loc /= 61) .or. &
            (epsabs <= 0.0_dp .and. &
             epsrel < max(50.0_dp*epmach, 0.5e-28_dp))) then
            ierr = 6
            return
        end if

        call qag_adapt(f, a, b, epsabs, epsrel, key_loc, limit_loc, result, &
                       abserr, ierr)
    end subroutine integrate_gk

    ! Single-level attempt; falls through to bisection if not converged.
    subroutine qag_adapt(f, a, b, epsabs, epsrel, key, limit, result, abserr, &
                         ierr)
        procedure(gk_integrand_t) :: f
        real(dp), intent(in)  :: a, b, epsabs, epsrel
        integer,  intent(in)  :: key, limit
        real(dp), intent(out) :: result, abserr
        integer,  intent(out) :: ierr

        real(dp) :: resabs, resasc, errbnd

        ierr = 0
        call gk_apply(f, key, a, b, result, abserr, resabs, resasc)

        errbnd = max(epsabs, epsrel*abs(result))
        if (abserr <= 50.0_dp*epmach*resabs .and. abserr > errbnd) ierr = 2
        if (limit == 1 .and. abserr > errbnd) ierr = 1
        if (ierr /= 0 .or. &
            (abserr <= errbnd .and. abserr /= resasc) .or. &
            abserr == 0.0_dp) return

        call qag_bisect(f, a, b, epsabs, epsrel, key, limit, result, abserr, &
                        ierr)
    end subroutine qag_adapt

    ! Adaptive bisection workspace; only reached when the single-panel rule
    ! did not converge, so the allocation cost is not on the fast path.
    subroutine qag_bisect(f, a, b, epsabs, epsrel, key, limit, result, abserr, &
                          ierr)
        procedure(gk_integrand_t) :: f
        real(dp), intent(in)    :: a, b, epsabs, epsrel
        integer,  intent(in)    :: key, limit
        real(dp), intent(inout) :: result, abserr
        integer,  intent(out)   :: ierr

        real(dp) :: ws(limit, 4)
        real(dp) :: errbnd, errmax, area, errsum
        real(dp) :: a1, b1, a2, b2
        real(dp) :: area1, area2, area12, error1, error2, erro12
        real(dp) :: resabs1, resabs2, defab1, defab2
        integer  :: last, maxerr, iroff1, iroff2

        ierr   = 0
        ws(1,1) = a
        ws(1,2) = b
        ws(1,3) = result
        ws(1,4) = abserr

        area   = result
        errsum = abserr
        iroff1 = 0
        iroff2 = 0
        do last = 2, limit
            maxerr = maxloc(ws(1:last-1, 4), 1)
            errmax = ws(maxerr, 4)
            a1     = ws(maxerr, 1)
            b1     = 0.5_dp*(ws(maxerr, 1) + ws(maxerr, 2))
            a2     = b1
            b2     = ws(maxerr, 2)
            call gk_apply(f, key, a1, b1, area1, error1, resabs1, defab1)
            call gk_apply(f, key, a2, b2, area2, error2, resabs2, defab2)

            area12 = area1 + area2
            erro12 = error1 + error2
            errsum = errsum + erro12 - errmax
            area   = area + area12 - ws(maxerr, 3)
            if (defab1 /= error1 .and. defab2 /= error2) then
                if (abs(ws(maxerr, 3) - area12) <= 1.0e-5_dp*abs(area12) .and. &
                    erro12 >= 0.99_dp*errmax) iroff1 = iroff1 + 1
                if (last > 10 .and. erro12 > errmax) iroff2 = iroff2 + 1
            end if
            ws(maxerr, 3) = area1
            ws(maxerr, 4) = error1
            ws(maxerr, 2) = b1
            ws(last,   1) = a2
            ws(last,   2) = b2
            ws(last,   3) = area2
            ws(last,   4) = error2

            errbnd = max(epsabs, epsrel*abs(area))
            if (errsum <= errbnd) exit
            if (iroff1 >= 6 .or. iroff2 >= 20) ierr = 2
            if (last == limit) ierr = 1
            if (max(abs(a1), abs(b2)) <= &
                (1.0_dp + 100.0_dp*epmach)*(abs(a2) + 1000.0_dp*uflow)) ierr = 3
            if (ierr /= 0) exit
        end do

        ! Sum fresh to avoid accumulated cancellation from the running update.
        result = sum(ws(1:min(last, limit), 3))
        abserr = sum(ws(1:min(last, limit), 4))
    end subroutine qag_bisect

    subroutine gk_apply(f, key, a, b, result, abserr, resabs, resasc)
        procedure(gk_integrand_t) :: f
        integer,  intent(in)  :: key
        real(dp), intent(in)  :: a, b
        real(dp), intent(out) :: result, abserr, resabs, resasc

        select case (key)
        case (15)
            call gk_rule(f, a, b, 8,  xgk15, wgk15, 4,  wg15, .true.,  &
                         result, abserr, resabs, resasc)
        case (21)
            call gk_rule(f, a, b, 11, xgk21, wgk21, 5,  wg21, .false., &
                         result, abserr, resabs, resasc)
        case (31)
            call gk_rule(f, a, b, 16, xgk31, wgk31, 8,  wg31, .true.,  &
                         result, abserr, resabs, resasc)
        case (61)
            call gk_rule(f, a, b, 31, xgk61, wgk61, 15, wg61, .false., &
                         result, abserr, resabs, resasc)
        end select
    end subroutine gk_apply

    ! Apply one GK(nh) panel on [a,b].
    ! Kronrod nodes interlace Gauss: even entries in xgk are Gauss nodes
    ! (QUADPACK layout). center_gauss marks odd-n rules whose center node
    ! also carries a Gauss weight.
    subroutine gk_rule(f, a, b, nh, xgk, wgk, ngw, wg, center_gauss, &
                       result, abserr, resabs, resasc)
        procedure(gk_integrand_t) :: f
        integer,  intent(in)  :: nh, ngw
        real(dp), intent(in)  :: a, b, xgk(nh), wgk(nh), wg(ngw)
        logical,  intent(in)  :: center_gauss
        real(dp), intent(out) :: result, abserr, resabs, resasc

        ! Fixed upper bound at G30K61 pair count; gfortran keeps these on the
        ! stack, avoiding heap traffic in the hot-path single-panel case.
        real(dp) :: fv1(30), fv2(30)
        real(dp) :: centr, hlgth, dhlgth, fc, absc, fval1, fval2, fsum
        real(dp) :: resg, resk, reskh
        integer  :: j, jtw, jtwm1, npair_gauss

        centr  = 0.5_dp*(a + b)
        hlgth  = 0.5_dp*(b - a)
        dhlgth = abs(hlgth)

        fc = f(centr)
        if (center_gauss) then
            resg        = wg(ngw)*fc
            npair_gauss = ngw - 1
        else
            resg        = 0.0_dp
            npair_gauss = ngw
        end if
        resk   = wgk(nh)*fc
        resabs = abs(resk)

        ! Even-indexed Kronrod nodes coincide with Gauss nodes.
        do j = 1, npair_gauss
            jtw  = 2*j
            absc  = hlgth*xgk(jtw)
            fval1 = f(centr - absc)
            fval2 = f(centr + absc)
            fv1(jtw) = fval1
            fv2(jtw) = fval2
            fsum  = fval1 + fval2
            resg   = resg + wg(j)*fsum
            resk   = resk + wgk(jtw)*fsum
            resabs = resabs + wgk(jtw)*(abs(fval1) + abs(fval2))
        end do
        ! Odd-indexed Kronrod-only nodes.
        do j = 1, nh/2
            jtwm1 = 2*j - 1
            absc   = hlgth*xgk(jtwm1)
            fval1  = f(centr - absc)
            fval2  = f(centr + absc)
            fv1(jtwm1) = fval1
            fv2(jtwm1) = fval2
            fsum   = fval1 + fval2
            resk   = resk + wgk(jtwm1)*fsum
            resabs = resabs + wgk(jtwm1)*(abs(fval1) + abs(fval2))
        end do

        reskh  = 0.5_dp*resk
        resasc = wgk(nh)*abs(fc - reskh)
        do j = 1, nh - 1
            resasc = resasc + &
                wgk(j)*(abs(fv1(j) - reskh) + abs(fv2(j) - reskh))
        end do

        result = resk*hlgth
        resabs = resabs*dhlgth
        resasc = resasc*dhlgth
        abserr = abs((resk - resg)*hlgth)
        if (resasc /= 0.0_dp .and. abserr /= 0.0_dp) then
            abserr = resasc*min(1.0_dp, (200.0_dp*abserr/resasc)**1.5_dp)
        end if
        if (resabs > uflow/(50.0_dp*epmach)) then
            abserr = max(abserr, 50.0_dp*epmach*resabs)
        end if
    end subroutine gk_rule

end module fortnum_integrate_gk
