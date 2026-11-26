module fortnum_lapack_ggev
    ! Complex generalized eigenproblem A x = lambda B x with reference LAPACK
    ! arguments and info codes (prefix fortnum_): the driver ZGGEV. B = Q R
    ! by Householder QR, (Q^H A, R) to Hessenberg-triangular form by Givens
    ! rotations (LAPACK's ZGGHRD), the single-shift QZ algorithm with
    ! deflation of infinite eigenvalues (Moler-Stewart, LAPACK's ZHGEQZ) and
    ! eigenvectors by substitution in the triangular pencil (LAPACK's
    ! ZTGEVC). Eigenvalues are pairs (alpha, beta), lambda = alpha/beta,
    ! beta real and nonnegative.
    use fortnum_kinds, only: dp
    use fortnum_lapack_qr, only: fortnum_zgeqrf, fortnum_zunmqr, fortnum_zungqr
    use fortnum_lapack_gees, only: fortnum_zlartg
    implicit none
    private

    public :: fortnum_zggev

    real(dp), parameter :: ulp = epsilon(1.0_dp)
    real(dp), parameter :: safmin = tiny(1.0_dp)
    complex(dp), parameter :: czero = (0.0_dp, 0.0_dp), cone = (1.0_dp, 0.0_dp)

contains

    real(dp) elemental function cabs1(x)
        complex(dp), intent(in) :: x
        cabs1 = abs(real(x, dp)) + abs(aimag(x))
    end function cabs1

    ! [x; y] := [c s; -conj(s) c] [x; y], LAPACK's ZROT
    pure subroutine rot(x, y, c, s)
        complex(dp), intent(inout) :: x(:), y(:)
        real(dp), intent(in) :: c
        complex(dp), intent(in) :: s
        complex(dp) :: t(size(x))
        t = c*x + s*y
        y = c*y - conjg(s)*x
        x = t
    end subroutine rot

    !> Generalized eigenvalues lambda_j = alpha_j/beta_j of the pencil (A, B)
    !> and, with jobvl/jobvr 'V', left (u^H A = lambda u^H B) and right
    !> (A v = lambda B v) eigenvectors, each scaled so that its largest
    !> component has |re| + |im| = 1. beta_j = 0 marks an infinite eigenvalue.
    subroutine fortnum_zggev(jobvl, jobvr, n, a, lda, b, ldb, alpha, beta, vl, ldvl, &
            vr, ldvr, work, lwork, rwork, info)
        character, intent(in) :: jobvl, jobvr
        integer, intent(in) :: n, lda, ldb, ldvl, ldvr, lwork
        complex(dp), intent(inout) :: a(lda, *), b(ldb, *)
        complex(dp), intent(out) :: alpha(*), beta(*), vl(ldvl, *), vr(ldvr, *), work(*)
        real(dp), intent(out) :: rwork(*)
        integer, intent(out) :: info
        logical :: ilvl, ilvr, ilv, ilascl, ilbscl
        complex(dp), allocatable :: tau(:), wk(:), q(:, :), z(:, :)
        real(dp) :: anrm, bnrm, anrmto, bnrmto, smlnum, bignum, temp
        integer :: minwrk, j, ierr
        info = 0
        ilvl = any(jobvl == ['V', 'v'])
        ilvr = any(jobvr == ['V', 'v'])
        ilv = ilvl .or. ilvr
        minwrk = max(1, 2*n)
        if (.not. ilvl .and. .not. any(jobvl == ['N', 'n'])) then
            info = -1
        else if (.not. ilvr .and. .not. any(jobvr == ['N', 'n'])) then
            info = -2
        else if (n < 0) then
            info = -3
        else if (lda < max(1, n)) then
            info = -5
        else if (ldb < max(1, n)) then
            info = -7
        else if (ldvl < 1 .or. (ilvl .and. ldvl < n)) then
            info = -11
        else if (ldvr < 1 .or. (ilvr .and. ldvr < n)) then
            info = -13
        else if (lwork < minwrk .and. lwork /= -1) then
            info = -15
        end if
        if (info /= 0) return
        work(1) = minwrk
        if (lwork == -1 .or. n == 0) return
        rwork(1) = 0
        smlnum = sqrt(safmin)/ulp
        bignum = 1/smlnum
        ! scale A and B into [smlnum, bignum] in max norm
        anrm = 0
        bnrm = 0
        do j = 1, n
            anrm = max(anrm, maxval(abs(a(1:n, j))))
            bnrm = max(bnrm, maxval(abs(b(1:n, j))))
        end do
        call norm_scale(anrm, ilascl, anrmto)
        call norm_scale(bnrm, ilbscl, bnrmto)
        if (ilascl) then
            do j = 1, n
                a(1:n, j) = a(1:n, j)*(anrmto/anrm)
            end do
        end if
        if (ilbscl) then
            do j = 1, n
                b(1:n, j) = b(1:n, j)*(bnrmto/bnrm)
            end do
        end if
        ! B = Q R, A := Q^H A
        allocate (tau(max(1, n)), wk(max(1, n)))
        call fortnum_zgeqrf(n, n, b, ldb, tau, wk, size(wk), ierr)
        call fortnum_zunmqr('L', 'C', n, n, n, b, ldb, tau, a, lda, wk, size(wk), ierr)
        if (ilvl) then
            allocate (q(n, n))
            do j = 1, n
                q(1:n, j) = b(1:n, j)
            end do
            call fortnum_zungqr(n, n, n, q, n, tau, wk, size(wk), ierr)
        else
            allocate (q(1, 1))
        end if
        if (ilvr) then
            allocate (z(n, n))
            z = czero
            do j = 1, n
                z(j, j) = cone
            end do
        else
            allocate (z(1, 1))
        end if
        do j = 1, n - 1
            b(j + 1:n, j) = czero
        end do
        call gghrd(ilvl, ilvr, n, a, lda, b, ldb, q, z)
        call hgeqz(ilv, ilvl, ilvr, n, a, lda, b, ldb, alpha, beta, q, z, info)
        if (info /= 0) then
            if (info > n) return
        else if (ilv) then
            call tgevc(ilvl, ilvr, n, a, lda, b, ldb, q, z)
            if (ilvl) then
                do j = 1, n
                    temp = maxval(cabs1(q(1:n, j)))
                    if (temp < smlnum) then
                        vl(1:n, j) = q(1:n, j)
                    else
                        vl(1:n, j) = q(1:n, j)/temp
                    end if
                end do
            end if
            if (ilvr) then
                do j = 1, n
                    temp = maxval(cabs1(z(1:n, j)))
                    if (temp < smlnum) then
                        vr(1:n, j) = z(1:n, j)
                    else
                        vr(1:n, j) = z(1:n, j)/temp
                    end if
                end do
            end if
        end if
        if (ilascl) alpha(1:n) = alpha(1:n)*(anrm/anrmto)
        if (ilbscl) beta(1:n) = beta(1:n)*(bnrm/bnrmto)
    contains
        subroutine norm_scale(nrm, scaled, target)
            real(dp), intent(in) :: nrm
            logical, intent(out) :: scaled
            real(dp), intent(out) :: target
            scaled = .false.
            target = nrm
            if (nrm > 0 .and. nrm < smlnum) then
                scaled = .true.
                target = smlnum
            else if (nrm > bignum) then
                scaled = .true.
                target = bignum
            end if
        end subroutine norm_scale
    end subroutine fortnum_zggev

    ! (A, B) with B upper triangular to H = Q1^H A Z1 upper Hessenberg and
    ! T = Q1^H B Z1 upper triangular by Givens rotations; Q := Q Q1 (wantq),
    ! Z := Z Z1 (wantz).
    subroutine gghrd(wantq, wantz, n, a, lda, b, ldb, q, z)
        logical, intent(in) :: wantq, wantz
        integer, intent(in) :: n, lda, ldb
        complex(dp), intent(inout) :: a(lda, *), b(ldb, *), q(:, :), z(:, :)
        complex(dp) :: f
        real(dp) :: c
        complex(dp) :: s
        integer :: jcol, jrow
        do jcol = 1, n - 2
            do jrow = n, jcol + 2, -1
                ! rows jrow-1, jrow: zero a(jrow, jcol)
                f = a(jrow - 1, jcol)
                call fortnum_zlartg(f, a(jrow, jcol), c, s, a(jrow - 1, jcol))
                a(jrow, jcol) = czero
                call rot(a(jrow - 1, jcol + 1:n), a(jrow, jcol + 1:n), c, s)
                call rot(b(jrow - 1, jrow - 1:n), b(jrow, jrow - 1:n), c, s)
                if (wantq) call rot(q(:, jrow - 1), q(:, jrow), c, conjg(s))
                ! columns jrow, jrow-1: zero the fill b(jrow, jrow-1)
                f = b(jrow, jrow)
                call fortnum_zlartg(f, b(jrow, jrow - 1), c, s, b(jrow, jrow))
                b(jrow, jrow - 1) = czero
                call rot(a(1:n, jrow), a(1:n, jrow - 1), c, s)
                call rot(b(1:jrow - 1, jrow), b(1:jrow - 1, jrow - 1), c, s)
                if (wantz) call rot(z(:, jrow), z(:, jrow - 1), c, s)
            end do
        end do
    end subroutine gghrd

    ! Single-shift QZ on the Hessenberg-triangular pencil (H, T) = (a, b):
    ! H := Q2^H H Z2 and T := Q2^H T Z2 upper triangular when wants (the
    ! generalized Schur form), only the eigenvalues otherwise. T's diagonal
    ! is made real and nonnegative. info = j > 0: no convergence, the pairs
    ! j+1..n are correct.
    subroutine hgeqz(wants, wantq, wantz, n, h, ldh, t, ldt, alpha, beta, q, z, info)
        logical, intent(in) :: wants, wantq, wantz
        integer, intent(in) :: n, ldh, ldt
        complex(dp), intent(inout) :: h(ldh, *), t(ldt, *), q(:, :), z(:, :)
        complex(dp), intent(out) :: alpha(*), beta(*)
        integer, intent(out) :: info
        real(dp) :: anorm, bnorm, btol, ascale, bscale, c, temp, temp2
        complex(dp) :: s, f, shift, eshift, u12, ad11, ad21, ad12, ad22, abi12, abi22, x, y, ctemp
        integer :: ilast, ifirst, istart, ifrstm, ilastm, iiter, maxit, jiter, j, jch, jr, jc
        logical :: ilazro, handled, chased
        info = 0
        anorm = 0
        bnorm = 0
        do j = 1, n
            anorm = anorm + sum(abs(h(1:min(j + 1, n), j))**2)
            bnorm = bnorm + sum(abs(t(1:j, j))**2)
        end do
        anorm = sqrt(anorm)
        bnorm = sqrt(bnorm)
        btol = max(safmin, ulp*bnorm)
        ascale = 1/max(safmin, anorm)
        bscale = 1/max(safmin, bnorm)
        ilast = n
        ifrstm = 1
        ilastm = n
        iiter = 0
        eshift = czero
        maxit = 30*n
        do jiter = 1, maxit + 1
            if (ilast < 1) exit
            if (jiter > maxit) then
                info = ilast
                return
            end if
            if (.not. wants) ilastm = ilast
            ! is the bottom eigenvalue split off, or t(ilast, ilast) zero?
            if (ilast == 1) then
                call deflate()
                cycle
            end if
            if (cabs1(h(ilast, ilast - 1)) <= max(safmin, ulp*(cabs1(h(ilast, ilast)) &
                    + cabs1(h(ilast - 1, ilast - 1))))) then
                h(ilast, ilast - 1) = czero
                call deflate()
                cycle
            end if
            if (abs(t(ilast, ilast)) <= btol) then
                t(ilast, ilast) = czero
                call infinite_at_bottom()
                call deflate()
                cycle
            end if
            ! find ifirst: the top of the unreduced block ending at ilast
            handled = .false.
            ifirst = 1
            do j = ilast - 1, 1, -1
                if (j == 1) then
                    ilazro = .true.
                else if (cabs1(h(j, j - 1)) <= max(safmin, ulp*(cabs1(h(j, j)) &
                        + cabs1(h(j - 1, j - 1))))) then
                    h(j, j - 1) = czero
                    ilazro = .true.
                else
                    ilazro = .false.
                end if
                if (abs(t(j, j)) < btol) then
                    t(j, j) = czero
                    chased = .true.
                    if (ilazro) then
                        ! the zero of T at the top of the block moves down by
                        ! rotations of H's rows
                        do jch = j, ilast - 1
                            f = h(jch, jch)
                            call fortnum_zlartg(f, h(jch + 1, jch), c, s, h(jch, jch))
                            h(jch + 1, jch) = czero
                            call rot(h(jch, jch + 1:ilastm), h(jch + 1, jch + 1:ilastm), c, s)
                            call rot(t(jch, jch + 1:ilastm), t(jch + 1, jch + 1:ilastm), c, s)
                            if (wantq) call rot(q(:, jch), q(:, jch + 1), c, conjg(s))
                            if (abs(t(jch + 1, jch + 1)) >= btol) then
                                if (jch + 1 >= ilast) then
                                    call deflate()
                                else
                                    ifirst = jch + 1
                                    call sweep()
                                end if
                                chased = .false.
                                exit
                            end if
                            t(jch + 1, jch + 1) = czero
                        end do
                    else
                        ! chase the zero of T to the bottom
                        do jch = j, ilast - 1
                            f = t(jch, jch + 1)
                            call fortnum_zlartg(f, t(jch + 1, jch + 1), c, s, t(jch, jch + 1))
                            t(jch + 1, jch + 1) = czero
                            if (jch < ilastm - 1) call rot(t(jch, jch + 2:ilastm), &
                                t(jch + 1, jch + 2:ilastm), c, s)
                            call rot(h(jch, jch - 1:ilastm), h(jch + 1, jch - 1:ilastm), c, s)
                            if (wantq) call rot(q(:, jch), q(:, jch + 1), c, conjg(s))
                            f = h(jch + 1, jch)
                            call fortnum_zlartg(f, h(jch + 1, jch - 1), c, s, h(jch + 1, jch))
                            h(jch + 1, jch - 1) = czero
                            call rot(h(ifrstm:jch, jch), h(ifrstm:jch, jch - 1), c, s)
                            call rot(t(ifrstm:jch - 1, jch), t(ifrstm:jch - 1, jch - 1), c, s)
                            if (wantz) call rot(z(:, jch), z(:, jch - 1), c, s)
                        end do
                    end if
                    if (chased) then
                        call infinite_at_bottom()
                        call deflate()
                    end if
                    handled = .true.
                    exit
                end if
                if (ilazro) then
                    ifirst = j
                    exit
                end if
            end do
            if (.not. handled) call sweep()
        end do
    contains
        ! t(ilast, ilast) = 0: a column rotation zeroes h(ilast, ilast-1)
        subroutine infinite_at_bottom()
            if (ilast == 1) return
            f = h(ilast, ilast)
            call fortnum_zlartg(f, h(ilast, ilast - 1), c, s, h(ilast, ilast))
            h(ilast, ilast - 1) = czero
            call rot(h(ifrstm:ilast - 1, ilast), h(ifrstm:ilast - 1, ilast - 1), c, s)
            call rot(t(ifrstm:ilast - 1, ilast), t(ifrstm:ilast - 1, ilast - 1), c, s)
            if (wantz) call rot(z(:, ilast), z(:, ilast - 1), c, s)
        end subroutine infinite_at_bottom

        ! (h(ilast,ilast), t(ilast,ilast)) is an eigenvalue: make t's entry
        ! real and nonnegative, store it and shrink the active block
        subroutine deflate()
            complex(dp) :: signbc
            real(dp) :: absb
            absb = abs(t(ilast, ilast))
            if (absb > safmin) then
                signbc = conjg(t(ilast, ilast)/absb)
                t(ilast, ilast) = absb
                if (wants) then
                    t(ifrstm:ilast - 1, ilast) = t(ifrstm:ilast - 1, ilast)*signbc
                    h(ifrstm:ilast, ilast) = h(ifrstm:ilast, ilast)*signbc
                else
                    h(ilast, ilast) = h(ilast, ilast)*signbc
                end if
                if (wantz) z(:, ilast) = z(:, ilast)*signbc
            else
                t(ilast, ilast) = czero
            end if
            alpha(ilast) = h(ilast, ilast)
            beta(ilast) = t(ilast, ilast)
            ilast = ilast - 1
            if (ilast < 1) return
            iiter = 0
            eshift = czero
            if (.not. wants) then
                ilastm = ilast
                if (ifrstm > ilast) ifrstm = 1
            end if
        end subroutine deflate

        ! one implicit single-shift QZ sweep over rows ifirst..ilast
        subroutine sweep()
            integer :: k
            iiter = iiter + 1
            if (.not. wants) ifrstm = ifirst
            if (mod(iiter, 10) /= 0) then
                ! the eigenvalue of the trailing 2 x 2 of H T^-1 nearer to
                ! h(ilast,ilast)/t(ilast,ilast) (Wilkinson)
                u12 = (bscale*t(ilast - 1, ilast))/(bscale*t(ilast, ilast))
                ad11 = (ascale*h(ilast - 1, ilast - 1))/(bscale*t(ilast - 1, ilast - 1))
                ad21 = (ascale*h(ilast, ilast - 1))/(bscale*t(ilast - 1, ilast - 1))
                ad12 = (ascale*h(ilast - 1, ilast))/(bscale*t(ilast, ilast))
                ad22 = (ascale*h(ilast, ilast))/(bscale*t(ilast, ilast))
                abi22 = ad22 - u12*ad21
                abi12 = ad12 - u12*ad11
                shift = abi22
                ctemp = sqrt(abi12)*sqrt(ad21)
                temp = cabs1(ctemp)
                if (ctemp /= czero) then
                    x = (ad11 - shift)/2
                    temp2 = cabs1(x)
                    temp = max(temp, temp2)
                    y = temp*sqrt((x/temp)**2 + (ctemp/temp)**2)
                    if (temp2 > 0) then
                        if (real(x/temp2, dp)*real(y, dp) + aimag(x/temp2)*aimag(y) < 0) y = -y
                    end if
                    shift = shift - ctemp*(ctemp/(x + y))
                end if
            else
                ! exceptional shift
                if (mod(iiter, 20) == 0 .and. bscale*abs(t(ilast, ilast)) > safmin) then
                    eshift = eshift + (ascale*h(ilast, ilast))/(bscale*t(ilast, ilast))
                else
                    eshift = eshift + (ascale*h(ilast, ilast - 1))/(bscale*t(ilast - 1, ilast - 1))
                end if
                shift = eshift
            end if
            istart = ifirst
            call fortnum_zlartg(ascale*h(istart, istart) - shift*(bscale*t(istart, istart)), &
                ascale*h(istart + 1, istart), c, s, ctemp)
            do k = istart, ilast - 1
                if (k > istart) then
                    f = h(k, k - 1)
                    call fortnum_zlartg(f, h(k + 1, k - 1), c, s, h(k, k - 1))
                    h(k + 1, k - 1) = czero
                end if
                jc = k
                call rot(h(k, jc:ilastm), h(k + 1, jc:ilastm), c, s)
                call rot(t(k, jc:ilastm), t(k + 1, jc:ilastm), c, s)
                if (wantq) call rot(q(:, k), q(:, k + 1), c, conjg(s))
                f = t(k + 1, k + 1)
                call fortnum_zlartg(f, t(k + 1, k), c, s, t(k + 1, k + 1))
                t(k + 1, k) = czero
                jr = min(k + 2, ilast)
                call rot(h(ifrstm:jr, k + 1), h(ifrstm:jr, k), c, s)
                call rot(t(ifrstm:k, k + 1), t(ifrstm:k, k), c, s)
                if (wantz) call rot(z(:, k + 1), z(:, k), c, s)
            end do
        end subroutine sweep
    end subroutine hgeqz

    ! Eigenvectors of the triangular pencil (S, P) = (a, b): right x with
    ! (beta S - alpha P) x = 0 and left y with y^H (beta S - alpha P) = 0,
    ! (alpha, beta) = (S(k,k), P(k,k)), multiplied into z and q.
    subroutine tgevc(wantl, wantr, n, s, lds, p, ldp, q, z)
        logical, intent(in) :: wantl, wantr
        integer, intent(in) :: n, lds, ldp
        complex(dp), intent(in) :: s(lds, *), p(ldp, *)
        complex(dp), intent(inout) :: q(:, :), z(:, :)
        complex(dp), allocatable :: x(:), vout(:, :)
        complex(dp) :: acoef, bcoef, d, salpha
        real(dp) :: anorm, bnorm, ascale, bscale, temp, dmin, sbeta, bignum
        integer :: je, j
        anorm = 0
        bnorm = 0
        do j = 1, n
            anorm = max(anorm, sum(cabs1(s(1:j, j))))
            bnorm = max(bnorm, sum(cabs1(p(1:j, j))))
        end do
        ascale = 1/max(anorm, safmin)
        bscale = 1/max(bnorm, safmin)
        bignum = 1/(safmin*n/ulp)
        allocate (x(n))
        if (wantr) then
            allocate (vout(n, n))
            do je = n, 1, -1
                x = czero
                if (cabs1(s(je, je)) <= safmin .and. abs(real(p(je, je), dp)) <= safmin) then
                    ! singular pencil: a unit vector
                    x(je) = cone
                    vout(:, je) = matmul(z(:, 1:n), x)
                    cycle
                end if
                temp = 1/max(cabs1(s(je, je))*ascale, abs(real(p(je, je), dp))*bscale, safmin)
                salpha = (temp*s(je, je))*ascale
                sbeta = (temp*real(p(je, je), dp))*bscale
                acoef = sbeta*ascale
                bcoef = salpha*bscale
                dmin = max(ulp*cabs1(acoef)*anorm, ulp*cabs1(bcoef)*bnorm, safmin)
                x(je) = cone
                do j = je - 1, 1, -1
                    x(j) = -sum((acoef*s(j, j + 1:je) - bcoef*p(j, j + 1:je))*x(j + 1:je))
                    d = acoef*s(j, j) - bcoef*p(j, j)
                    if (cabs1(d) <= dmin) d = dmin
                    x(j) = x(j)/d
                    if (cabs1(x(j)) > bignum) x(j:je) = x(j:je)/cabs1(x(j))
                end do
                vout(:, je) = matmul(z(:, 1:je), x(1:je))
            end do
            z = vout
            deallocate (vout)
        end if
        if (wantl) then
            allocate (vout(n, n))
            do je = 1, n
                x = czero
                if (cabs1(s(je, je)) <= safmin .and. abs(real(p(je, je), dp)) <= safmin) then
                    x(je) = cone
                    vout(:, je) = matmul(q(:, 1:n), x)
                    cycle
                end if
                temp = 1/max(cabs1(s(je, je))*ascale, abs(real(p(je, je), dp))*bscale, safmin)
                salpha = (temp*s(je, je))*ascale
                sbeta = (temp*real(p(je, je), dp))*bscale
                acoef = sbeta*ascale
                bcoef = salpha*bscale
                dmin = max(ulp*cabs1(acoef)*anorm, ulp*cabs1(bcoef)*bnorm, safmin)
                x(je) = cone
                do j = je + 1, n
                    x(j) = -sum(conjg(acoef*s(je:j - 1, j) - bcoef*p(je:j - 1, j))*x(je:j - 1))
                    d = conjg(acoef*s(j, j) - bcoef*p(j, j))
                    if (cabs1(d) <= dmin) d = dmin
                    x(j) = x(j)/d
                    if (cabs1(x(j)) > bignum) x(je:j) = x(je:j)/cabs1(x(j))
                end do
                vout(:, je) = matmul(q(:, je:n), x(je:n))
            end do
            q = vout
        end if
    end subroutine tgevc

end module fortnum_lapack_ggev
