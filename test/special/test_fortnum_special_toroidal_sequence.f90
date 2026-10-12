program test_fortnum_special_toroidal_sequence
    ! Q_(n-1/2)(x) sequences, derivatives, log split and batches.
    ! Independent oracles: closed forms through fortnum_special_elliptic for
    ! n = 0, 1; mpmath legenq/legenp (type=3, 45 digits, exact double x) for
    ! the tables; the Wronskian, degree recurrence and centred finite
    ! differences as behavioral identities.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_special_elliptic, only: ellip_e_from_complement, &
        ellip_k_from_complement
    use fortnum_special_toroidal, only: toroidal_p_derivative, toroidal_q, &
        toroidal_q_log_split, toroidal_q_log_split_batch, &
        toroidal_q_sequence, toroidal_q_sequence_batch
    implicit none

    integer, parameter :: nx = 15, nn = 9, nxs = 7, nns = 5
    real(dp), parameter :: x_tab(nx) = [ &
        1.000000000001_dp, 1.000000001_dp, 1.000001_dp, 1.001_dp, 1.01_dp, &
        1.05_dp, 1.2_dp, 1.5_dp, 2.0_dp, 3.0_dp, 10.0_dp, 100.0_dp, &
        1000.0_dp, 10000.0_dp, 1000000.0_dp]
    integer, parameter :: n_tab(nn) = [0, 1, 2, 3, 5, 8, 12, 16, 20]
    real(dp), parameter :: x_split(nxs) = [ &
        1.000000000001_dp, 1.000001_dp, 1.001_dp, 1.05_dp, 1.5_dp, &
        3.0_dp, 100.0_dp]
    integer, parameter :: n_split(nns) = [0, 1, 2, 5, 10]
    real(dp), parameter :: q_ref_flat(135) = [ &
        1.554833406104685946e+1_dp, &
        1.3548334061054384296e+1_dp, &
        1.2881667394408291959e+1_dp, &
        1.2481667394441249057e+1_dp, &
        1.1973730886606317775e+1_dp, &
        1.1504733217846708391e+1_dp, &
        1.1099628384622330761e+1_dp, &
        1.0812072664915484282e+1_dp, &
        1.0588987642636412953e+1_dp, &
        1.2094500827116072346e+1_dp, &
        1.0094500832913323239e+1_dp, &
        9.4278341816384090937_dp, &
        9.0278342059579965455_dp, &
        8.5198977719584540513_dp, &
        8.0509002736495047303_dp, &
        7.6457957728731391313_dp, &
        7.3582404999188297238_dp, &
        7.1351560340380226756_dp, &
        8.6406222753454926764_dp, &
        6.6406263456563990575_dp, &
        5.9739698899284946646_dp, &
        5.5739855748435751281_dp, &
        5.0660953730522426075_dp, &
        4.5972010434188485477_dp, &
        4.1922909794677207816_dp, &
        3.904989250000913429_dp, &
        3.6822127814365521539_dp, &
        5.1862223889747841037_dp, &
        3.1885653771920873605_dp, &
        2.5269311271007774277_dp, &
        2.1339936668493522666_dp, &
        1.644830365017162476_dp, &
        1.2123681122548624279_dp, &
        8.6610515175120881809e-1_dp, &
        6.4428863737228745141e-1_dp, &
        4.9064820838812618997e-1_dp, &
        4.0316687795887194887_dp, &
        2.049318451462074005_dp, &
        1.4158592547726865215_dp, &
        1.0584374848354170358_dp, &
        6.5142626170330808332e-1_dp, &
        3.49529394741636429e-1_dp, &
        1.6603148605679188104e-1_dp, &
        8.2777050801863850802e-2_dp, &
        4.2422924337499742282e-2_dp, &
        3.2169704206797034663_dp, &
        1.2847427829660803469_dp, &
        7.2631642259261140635e-1_dp, &
        4.4936592017593900612e-1_dp, &
        1.9194076819749848964e-1_dp, &
        6.0315560288532423781e-2_dp, &
        1.4163588904379262426e-2_dp, &
        3.5053760465279088799e-3_dp, &
        8.9357575992726570298e-4_dp, &
        2.5009565105571013434_dp, &
        6.9955482843702442807e-1_dp, &
        2.8563555531353859568e-1_dp, &
        1.2868736913977942657e-1_dp, &
        2.9409355901862177021e-2_dp, &
        3.6463859784294715953e-3_dp, &
        2.4907688838397722224e-4_dp, &
        1.7972082430721266506e-5_dp, &
        1.3370297957803601561e-6_dp, &
        2.0189058199784232156_dp, &
        3.9317514837200473104e-1_dp, &
        1.1338169008453505689e-1_dp, &
        3.6210967179681297904e-2_dp, &
        4.1745653916465122335e-3_dp, &
        1.8607013040078070279e-4_dp, &
        3.2554125878326332115e-6_dp, &
        6.0213614396875759003e-8_dp, &
        1.148738859323242237e-9_dp, &
        1.6566381702365941664_dp, &
        2.2401429283641563704e-1_dp, &
        4.5158724151576976637e-2_dp, &
        1.009934158319694301e-2_dp, &
        5.7191641375056765816e-4_dp, &
        8.7891460374400615386e-6_dp, &
        3.7209455887334494544e-8_dp, &
        1.6659983884258140078e-10_dp, &
        7.6948665465124695227e-13_dp, &
        1.3110287771460599052_dp, &
        1.1288854241046769779e-1_dp, &
        1.4544577259850822759e-2_dp, &
        2.0808454010033305686e-3_dp, &
        4.8268030482453372046e-5_dp, &
        1.9461899387186309928e-7_dp, &
        1.384538402038867365e-10_dp, &
        1.0418849383567112851e-13_dp, &
        8.0886189209108563664e-17_dp, &
        7.0380587894745618557e-1_dp, &
        1.7644903169651926516e-2_dp, &
        6.6341594620695835336e-4_dp, &
        2.7713237520177744395e-5_dp, &
        5.4837831068257536889e-8_dp, &
        5.5115403491693234308e-12_dp, &
        2.855846800808853162e-17_dp, &
        1.5654530474747362175e-22_dp, &
        8.8533117814291858692e-28_dp, &
        2.2214831233847302284e-1_dp, &
        5.5538640149552274399e-4_dp, &
        2.0827532453582534892e-6_dp, &
        8.6783600069118691798e-9_dp, &
        1.7086384492823595123e-13_dp, &
        1.7044762110780853161e-20_dp, &
        8.7443536843966740346e-30_dp, &
        4.7458287729770294853e-39_dp, &
        2.6574117653097128168e-48_dp, &
        7.0248160481942087817e-2_dp, &
        1.756204505981198131e-5_dp, &
        6.5857686124738401832e-9_dp, &
        2.7440709555073857791e-12_dp, &
        5.4023924229891802737e-19_dp, &
        5.3888302077072857396e-29_dp, &
        2.7643191497571314511e-42_dp, &
        1.50013203930609267e-55_dp, &
        8.3991095955771413185e-69_dp, &
        2.2214414732443859008e-2_dp, &
        5.5536036987304751954e-7_dp, &
        2.0826013924473693383e-11_dp, &
        8.6775058240097569539e-16_dp, &
        1.7083839677328190814e-24_dp, &
        1.7040939538680074557e-37_dp, &
        8.7415166085731854976e-55_dp, &
        4.7438140928191699326e-72_dp, &
        2.6560178670717094135e-89_dp, &
        2.2214414690795996438e-3_dp, &
        5.5536036727005610605e-10_dp, &
        2.0826013772632527418e-16_dp, &
        8.6775057385991009959e-23_dp, &
        1.7083839422875610984e-35_dp, &
        1.7040939156480560094e-54_dp, &
        8.7415163249225323555e-80_dp, &
        4.7438138914016231934e-105_dp, &
        2.6560177277237588569e-130_dp &
        ]
    real(dp), parameter :: dq_ref_flat(135) = [ &
        -4.9995555366189103288e+11_dp, &
        -4.9995555365461686585e+11_dp, &
        -4.9995555363479436475e+11_dp, &
        -4.9995555360309019627e+11_dp, &
        -4.999555535055225714e+11_dp, &
        -4.9995555327723411981e+11_dp, &
        -4.9995555282616170101e+11_dp, &
        -4.9995555221335369919e+11_dp, &
        -4.9995555144347358627e+11_dp, &
        -4.9999995995413052625e+8_dp, &
        -4.9999995440688011402e+8_dp, &
        -4.9999993976512886042e+8_dp, &
        -4.9999991669554338746e+8_dp, &
        -4.9999984675858376273e+8_dp, &
        -4.9999968581987631218e+8_dp, &
        -4.9999937290076873644e+8_dp, &
        -4.999989535073808079e+8_dp, &
        -4.9999843230316391428e+8_dp, &
        -5.0000089261852918769e+5_dp, &
        -4.9999707230828409234e+5_dp, &
        -4.9998761136583787489e+5_dp, &
        -4.9997317642872437836e+5_dp, &
        -4.9993087009234002713e+5_dp, &
        -4.9983727955957971403e+5_dp, &
        -4.9966250457007095141e+5_dp, &
        -4.9943650039446499446e+5_dp, &
        -4.9916391922601476233e+5_dp, &
        -5.0046057825384598894e+2_dp, &
        -4.9836792763761238776e+2_dp, &
        -4.9408344749946181422e+2_dp, &
        -4.8826020312920927805e+2_dp, &
        -4.7340434451127724437e+2_dp, &
        -4.4617450411528321344e+2_dp, &
        -4.0559469641179550988e+2_dp, &
        -3.6398389910906921262e+2_dp, &
        -3.2363174331749369094e+2_dp, &
        -5.0315099898570420621e+1_dp, &
        -4.8802416507761765529e+1_dp, &
        -4.6216462995646272611e+1_dp, &
        -4.3138979488671019444e+1_dp, &
        -3.656873258313750861e+1_dp, &
        -2.7381862352488578048e+1_dp, &
        -1.781853774134224086e+1_dp, &
        -1.1257413985798232879e+1_dp, &
        -6.9795367584993330867_dp, &
        -1.0210127603646861128e+1_dp, &
        -9.1121487734893529042_dp, &
        -7.6406420377147004337_dp, &
        -6.2068830831189072788_dp, &
        -3.8863933902574504962_dp, &
        -1.7949637088214783333_dp, &
        -5.9981711419891648448e-1_dp, &
        -1.9243598417548368544e-1_dp, &
        -6.0250539804121322365e-2_dp, &
        -2.6154465729903382165_dp, &
        -1.8880576323098550719_dp, &
        -1.2163369161162893603_dp, &
        -7.4551541105570068917e-1_dp, &
        -2.5989216643269417988e-1_dp, &
        -4.8788355319203469912e-2_dp, &
        -4.8378490357468899966e-3_dp, &
        -4.5757651974682982467e-4_dp, &
        -4.210986548584748157e-5_dp, &
        -1.0540734326382520369_dp, &
        -5.7165723896816644761e-1_dp, &
        -2.6772313589424257485e-1_dp, &
        -1.1813047863002622006e-1_dp, &
        -2.1114525702496988253e-2_dp, &
        -1.4412618483819811011e-3_dp, &
        -3.6872127097178510174e-5_dp, &
        -8.9752427614930659267e-7_dp, &
        -2.1233689889253961678e-8_dp, &
        -5.1487700793946211598e-1_dp, &
        -2.0143493076062714873e-1_dp, &
        -6.6848422266630841886e-2_dp, &
        -2.080003415431924218e-2_dp, &
        -1.8393744527858146815e-3_dp, &
        -4.3501900135634169704e-5_dp, &
        -2.7012971076137724374e-7_dp, &
        -1.5942829035692036705e-9_dp, &
        -9.1408867186525807069e-12_dp, &
        -2.3876236181423200112e-1_dp, &
        -6.0772696869666050741e-2_dp, &
        -1.2985276993296605534e-2_dp, &
        -2.5943878302627597042e-3_dp, &
        -9.4332816701860649459e-5_dp, &
        -5.8683795247720618959e-7_dp, &
        -6.1331140655410834547e-10_dp, &
        -6.0888072260069734207e-13_dp, &
        -5.8709666953391292187e-16_dp, &
        -3.5456635789418737016e-2_dp, &
        -2.6634184204592773758e-3_dp, &
        -1.6682945011488398458e-4_dp, &
        -9.7546356314439623588e-6_dp, &
        -3.0325475401035928857e-8_dp, &
        -4.7097281273365468516e-12_dp, &
        -3.5884876549433909279e-17_dp, &
        -2.5963948201028037565e-22_dp, &
        -1.8242909733308008323e-27_dp, &
        -1.1108248748597763157e-3_dp, &
        -8.3313167411201494368e-6_dp, &
        -5.2071868730827692031e-8_dp, &
        -3.0375968713548019081e-10_dp, &
        -9.3980205429988991296e-15_dp, &
        -1.4488812500448951842e-21_dp, &
        -1.0931009691284412942e-30_dp, &
        -7.8310205526989133163e-40_dp, &
        -5.4479730101787126643e-49_dp, &
        -3.5124106584055098058e-5_dp, &
        -2.6343084054149107403e-8_dp, &
        -1.6464431135438347128e-11_dp, &
        -9.6042537466698826801e-15_dp, &
        -2.9713174421080588067e-21_dp, &
        -4.5805080940421751818e-31_dp, &
        -3.4554007313472407413e-44_dp, &
        -2.4752191388652215139e-57_dp, &
        -1.7218183485005456969e-70_dp, &
        -1.1107207449525985506e-6_dp, &
        -8.330405600160747908e-11_dp, &
        -5.2065035114896939006e-15_dp, &
        -3.0371270554872546621e-19_dp, &
        -9.396111873426111054e-28_dp, &
        -1.4484798684325612175e-40_dp, &
        -1.0926895817452287219e-57_dp, &
        -7.8272932934391694995e-75_dp, &
        -5.4448366553693824841e-92_dp, &
        -1.1107207345406328624e-9_dp, &
        -8.3304055090560480942e-16_dp, &
        -5.2065034431611689815e-22_dp, &
        -3.0371270085113937325e-28_dp, &
        -9.3961116825866756016e-41_dp, &
        -1.4484798283016120834e-59_dp, &
        -1.0926895406158839025e-84_dp, &
        -7.8272929208167070229e-110_dp, &
        -5.4448363418364928943e-135_dp &
        ]
    real(dp), parameter :: p_ref_flat(35) = [ &
        9.9999999999987498889e-1_dp, &
        1.0000000000003750333_dp, &
        1.0000000000018751667_dp, &
        1.0000000000123761001_dp, &
        1.0000000000498794339_dp, &
        9.9999987500003516652e-1_dp, &
        1.0000003749999413754_dp, &
        1.0000018750004100019_dp, &
        1.0000123750351904249_dp, &
        1.0000498756094092271_dp, &
        9.9987503514404765362e-1_dp, &
        1.0003749414233337966_dp, &
        1.0018754101049934758_dp, &
        1.0124102280785706688_dp, &
        1.0504875928851046737_dp, &
        9.9383639337310461411e-1_dp, &
        1.01860561501301738_dp, &
        1.0947690632271895209_dp, &
        1.7114024483248912999_dp, &
        5.4729798318347909312_dp, &
        9.45006330929758054e-1_dp, &
        1.1746724294455385146_dp, &
        2.0343427485811576779_dp, &
        2.1522333339356769294e+1_dp, &
        1.8364101325086368641e+3_dp, &
        8.3462684167407318628e-1_dp, &
        1.5973866051758863743_dp, &
        6.1113374734788544353_dp, &
        7.3296131280272493122e+2_dp, &
        3.4391016721566564163e+6_dp, &
        3.0091748588199264556e-1_dp, &
        9.003746661068960679_dp, &
        1.2003992489805674263e+3_dp, &
        5.852908935099731136e+9_dp, &
        1.3077930651638724956e+21_dp &
        ]
    real(dp), parameter :: r_ref_flat(35) = [ &
        1.3862943611198423277_dp, &
        -6.1370563887971448568e-1_dp, &
        -1.2803723055470517704_dp, &
        -2.188308813497740591_dp, &
        -2.8802166992927533583_dp, &
        1.3862943128331032041_dp, &
        -6.13705244019744675e-1_dp, &
        -1.2803725812443525937_dp, &
        -2.1883232688260409608_dp, &
        -2.8803102258357991006_dp, &
        1.3862460820441701186_dp, &
        -6.1331079917240993987e-1_dp, &
        -1.2806475073101664159_dp, &
        -2.2027853313168774589_dp, &
        -2.9745743986451770167_dp, &
        1.3838990945473332019_dp, &
        -5.940138795672241574e-1_dp, &
        -1.2929191295765582562_dp, &
        -2.9646378964754903006_dp, &
        -1.0065614407742610664e+1_dp, &
        1.36387734608316278_dp, &
        -4.2104573417967157093e-1_dp, &
        -1.2967172503870640685_dp, &
        -1.4913970107854810739e+1_dp, &
        -1.2729024813225157698e+3_dp, &
        1.3110287771460599052_dp, &
        1.1288854241046769779e-1_dp, &
        1.4544577259850822759e-2_dp, &
        4.8268030482453372046e-5_dp, &
        5.1410128069387601202e-9_dp, &
        8.0923421519279766074e-1_dp, &
        1.7566742084034018247e+1_dp, &
        2.3419625331328048338e+3_dp, &
        1.1418945351134194149e+10_dp, &
        2.5514863988643413649e+21_dp &
        ]
    real(dp), parameter :: q_ref(nn, nx) = reshape(q_ref_flat, [nn, nx])
    real(dp), parameter :: dq_ref(nn, nx) = reshape(dq_ref_flat, [nn, nx])
    real(dp), parameter :: p_ref(nns, nxs) = reshape(p_ref_flat, [nns, nxs])
    real(dp), parameter :: r_ref(nns, nxs) = reshape(r_ref_flat, [nns, nxs])
    integer, parameter :: nmax = 20

    integer :: nfail

    nfail = 0
    call test_closed_forms()
    call test_tables()
    call test_recurrence_and_wronskian()
    call test_finite_differences()
    call test_log_split()
    call test_switch_continuity()
    call test_batches_and_domain()

    if (nfail /= 0) then
        write (error_unit, "(a,i0,a)") "toroidal sequence: ", nfail, &
            " failure(s)"
        error stop 1
    end if
    print "(a)", "test_fortnum_special_toroidal_sequence: PASS"

contains

    subroutine expect(label, value, reference, rtol)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, reference, rtol

        if (ieee_is_nan(value) .or. &
            abs(value - reference) > rtol*abs(reference)) then
            nfail = nfail + 1
            write (error_unit, "(a,a,3es24.15)") "FAIL: ", label, value, &
                reference, (value - reference)/reference
        end if
    end subroutine expect

    subroutine test_closed_forms()
        ! Q_(-1/2) = sqrt(s) K(s) and Q_(1/2) = x sqrt(s) K(s) - 2 E(s)/sqrt(s)
        ! with s = 2/(x+1); the complement 1-s = (x-1)/(x+1) is formed
        ! exactly.  Q_(1/2) cancels against its own terms for large x, so the
        ! closed form is used up to x = 10 only.
        integer :: i
        real(dp) :: x, c, s, k, e, q(0:3), dq(0:3), qa, qb

        do i = 1, nx
            x = x_tab(i)
            c = (x - 1.0_dp)/(x + 1.0_dp)
            s = 2.0_dp/(x + 1.0_dp)
            k = ellip_k_from_complement(c)
            e = ellip_e_from_complement(c)
            call toroidal_q_sequence(3, x, q, dq)
            qa = sqrt(s)*k
            call expect("Q_-1/2 closed form", q(0), qa, 2.0e-14_dp)
            if (x <= 10.0_dp) then
                qb = x*sqrt(s)*k - 2.0_dp*e/sqrt(s)
                call expect("Q_1/2 closed form", q(1), qb, 2.0e-13_dp)
            end if
            ! dQ_(-1/2)/dx = -E(s)/(sqrt(s)(x^2-1)) (DLMF 19.4 chain rule)
            call expect("dQ_-1/2 closed form", dq(0), &
                -e/(sqrt(s)*(x - 1.0_dp)*(x + 1.0_dp)), 2.0e-13_dp)
        end do
    end subroutine test_closed_forms

    subroutine test_tables()
        integer :: i, j
        real(dp) :: q(0:nmax), dq(0:nmax)
        character(len=48) :: label

        do i = 1, nx
            call toroidal_q_sequence(nmax, x_tab(i), q, dq)
            do j = 1, nn
                write (label, "(a,i0,a,es9.2)") "Q n=", n_tab(j), " x-1=", &
                    x_tab(i) - 1.0_dp
                call expect(trim(label), q(n_tab(j)), q_ref(j, i), 5.0e-13_dp)
                write (label, "(a,i0,a,es9.2)") "dQ n=", n_tab(j), " x-1=", &
                    x_tab(i) - 1.0_dp
                call expect(trim(label), dq(n_tab(j)), dq_ref(j, i), &
                    1.0e-12_dp)
                write (label, "(a,i0,a,es9.2)") "toroidal_q n=", n_tab(j), &
                    " x-1=", x_tab(i) - 1.0_dp
                call expect(trim(label), toroidal_q(n_tab(j), 0, x_tab(i)), &
                    q_ref(j, i), 5.0e-13_dp)
            end do
        end do
        ! Regression: the order-1 Miller path lost its target scale for large x.
        call expect("Q_39/2^1(1e6)", toroidal_q(20, 1, 1.0e6_dp), &
            -5.44483634183377e-129_dp, 1.0e-12_dp)
    end subroutine test_tables

    subroutine test_recurrence_and_wronskian()
        ! (n+1/2) Q_(n+1) - 2 n x Q_n + (n-1/2) Q_(n-1) = 0 and
        ! P Q' - P' Q = 1/(1-x^2) (DLMF 14.2.4, mu = 0).  P and P' come from
        ! the independent toroidal_p path of the existing module.
        real(dp), parameter :: x_test(6) = [1.001_dp, 1.05_dp, 1.5_dp, &
            3.0_dp, 50.0_dp, 1.0e5_dp]
        integer :: i, n
        real(dp) :: x, q(0:nmax), dq(0:nmax), p, dp_dx, res, scale, w

        do i = 1, size(x_test)
            x = x_test(i)
            call toroidal_q_sequence(nmax, x, q, dq)
            do n = 1, nmax - 1
                res = (real(n, dp) + 0.5_dp)*q(n + 1) - &
                    2.0_dp*real(n, dp)*x*q(n) + (real(n, dp) - 0.5_dp)*q(n - 1)
                scale = 2.0_dp*real(n, dp)*x*q(n)
                if (abs(res) > 4.0e-14_dp*scale) then
                    nfail = nfail + 1
                    write (error_unit, "(a,i0,es10.2,es12.3)") &
                        "FAIL: recurrence n,x,res/scale ", n, x, res/scale
                end if
            end do
            if (x > 10.0_dp) cycle
            do n = 0, 10
                p = toroidal_p_value(n, x)
                dp_dx = toroidal_p_derivative(n, 0, x)
                w = p*dq(n) - dp_dx*q(n)
                call expect("Wronskian", w*(1.0_dp - x*x), 1.0_dp, 2.0e-10_dp)
            end do
        end do
    end subroutine test_recurrence_and_wronskian

    function toroidal_p_value(n, x) result(p)
        use fortnum_special_toroidal, only: toroidal_p
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: p

        p = toroidal_p(n, 0, x)
    end function toroidal_p_value

    subroutine test_finite_differences()
        ! Fourth-order centred differences of the sequence itself.
        real(dp), parameter :: x_test(6) = [1.001_dp, 1.05_dp, 1.5_dp, &
            3.0_dp, 20.0_dp, 1000.0_dp]
        integer :: i, n
        real(dp) :: x, h, qm2(0:nmax), qm1(0:nmax), qp1(0:nmax), &
            qp2(0:nmax), q(0:nmax), dq(0:nmax), fd

        do i = 1, size(x_test)
            x = x_test(i)
            h = 1.0e-3_dp*(x - 1.0_dp)
            call toroidal_q_sequence(nmax, x - 2.0_dp*h, qm2)
            call toroidal_q_sequence(nmax, x - h, qm1)
            call toroidal_q_sequence(nmax, x + h, qp1)
            call toroidal_q_sequence(nmax, x + 2.0_dp*h, qp2)
            call toroidal_q_sequence(nmax, x, q, dq)
            do n = 0, nmax
                fd = (qm2(n) - 8.0_dp*qm1(n) + 8.0_dp*qp1(n) - qp2(n))/ &
                    (12.0_dp*h)
                call expect("dQ vs finite difference", dq(n), fd, 1.0e-8_dp)
            end do
        end do
    end subroutine test_finite_differences

    subroutine test_log_split()
        integer :: i, j, n
        real(dp) :: x, ell, p(0:nmax), r(0:nmax), q(0:nmax)
        character(len=48) :: label

        do i = 1, nxs
            x = x_split(i)
            call toroidal_q_log_split(nmax, x, p, r)
            call toroidal_q_sequence(nmax, x, q)
            ell = -0.5_dp*log(0.5_dp*(x - 1.0_dp))
            do j = 1, nns
                write (label, "(a,i0,a,es9.2)") "split p n=", n_split(j), &
                    " x-1=", x - 1.0_dp
                call expect(trim(label), p(n_split(j)), p_ref(j, i), &
                    1.0e-12_dp)
                write (label, "(a,i0,a,es9.2)") "split r n=", n_split(j), &
                    " x-1=", x - 1.0_dp
                call expect(trim(label), r(n_split(j)), r_ref(j, i), &
                    1.0e-12_dp)
            end do
            do n = 0, nmax
                ! Q is recessive while L p grows, so the sum cancels for
                ! large n x; the budget is relative to the largest term.
                if (abs(ell*p(n) + r(n) - q(n)) > 1.0e-12_dp*(abs(ell*p(n)) &
                    + abs(r(n)) + abs(q(n)))) then
                    nfail = nfail + 1
                    write (error_unit, "(a,i0,es10.2)") &
                        "FAIL: split reassembly n,x-1 ", n, x - 1.0_dp
                end if
            end do
        end do
    end subroutine test_log_split

    subroutine test_switch_continuity()
        ! Forward/Miller switch at (nmax+1) acosh(x) = 1/2: Q is smooth there,
        ! so the two branches must agree after a first-order correction.
        real(dp) :: xs, d, qa(0:nmax), qb(0:nmax), qm(0:nmax), dqm(0:nmax)
        integer :: n

        xs = cosh(0.5_dp/real(nmax + 1, dp))
        d = 1.0e-3_dp*(xs - 1.0_dp)
        call toroidal_q_sequence(nmax, xs - d, qa)
        call toroidal_q_sequence(nmax, xs + d, qb)
        call toroidal_q_sequence(nmax, xs, qm, dqm)
        do n = 0, nmax
            call expect("switch continuity", qb(n) - qa(n), &
                2.0_dp*d*dqm(n), 1.0e-5_dp)
        end do
    end subroutine test_switch_continuity

    subroutine test_batches_and_domain()
        integer, parameter :: np = 200
        real(dp) :: x(np), q(0:nmax, np), dq(0:nmax, np), p(0:nmax, np), &
            r(0:nmax, np), q1(0:nmax), dq1(0:nmax), p1(0:nmax), r1(0:nmax)
        real(dp) :: nanv(0:2)
        integer :: i

        do i = 1, np
            x(i) = 1.0_dp + 10.0_dp**(-12.0_dp + 18.0_dp*real(i - 1, dp)/ &
                real(np - 1, dp))
        end do
        call toroidal_q_sequence_batch(nmax, x, q, dq)
        call toroidal_q_log_split_batch(6, x, p(0:6, :), r(0:6, :))
        do i = 1, np
            call toroidal_q_sequence(nmax, x(i), q1, dq1)
            if (any(abs(q(:, i) - q1) > 1.0e-14_dp*abs(q1)) .or. &
                any(abs(dq(:, i) - dq1) > 1.0e-14_dp*abs(dq1))) then
                nfail = nfail + 1
                write (error_unit, "(a,i0)") "FAIL: batch /= scalar ", i
            end if
            call toroidal_q_log_split(6, x(i), p1(0:6), r1(0:6))
            if (any(abs(p(0:6, i) - p1(0:6)) > 1.0e-14_dp*abs(p1(0:6))) .or. &
                any(abs(r(0:6, i) - r1(0:6)) > 1.0e-14_dp*abs(r1(0:6)))) then
                nfail = nfail + 1
                write (error_unit, "(a,i0)") "FAIL: split batch /= scalar ", i
            end if
        end do

        call toroidal_q_sequence(2, 1.0_dp, nanv)
        if (.not. all(ieee_is_nan(nanv))) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: x = 1 must give NaN"
        end if
        call toroidal_q_sequence(2, 0.5_dp, nanv)
        if (.not. all(ieee_is_nan(nanv))) then
            nfail = nfail + 1
            write (error_unit, "(a)") "FAIL: x < 1 must give NaN"
        end if
    end subroutine test_batches_and_domain

end program test_fortnum_special_toroidal_sequence
