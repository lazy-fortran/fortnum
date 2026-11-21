    ! Divide-and-conquer eigensolver of a symmetric tridiagonal matrix
    ! (Cuppen 1981; Gu and Eisenstat, SIAM J. Matrix Anal. Appl. 16 (1995)
    ! 172), ported from reference LAPACK 3.12.1: DSTEDC with DLAED0 (COMPZ =
    ! 'I' path), DLAED1, DLAED2, DLAED3, DLAED4, DLAED5, DLAED6 and DLAMRG.
    ! Leaves of at most STEDC_LEAF rows use the implicit QL iteration.
    !
    ! Copyright (c) 1992-2023 The University of Tennessee and The University of
    ! Tennessee Research Foundation. All rights reserved.
    ! Copyright (c) 2000-2023 The University of California Berkeley. All
    ! rights reserved.
    ! Copyright (c) 2006-2023 The University of Colorado Denver. All rights
    ! reserved.
    !
    ! Redistribution and use in source and binary forms, with or without
    ! modification, are permitted provided that the following conditions are
    ! met:
    ! - Redistributions of source code must retain the above copyright
    !   notice, this list of conditions and the following disclaimer.
    ! - Redistributions in binary form must reproduce the above copyright
    !   notice, this list of conditions and the following disclaimer listed
    !   in this license in the documentation and/or other materials provided
    !   with the distribution.
    ! - Neither the name of the copyright holders nor the names of its
    !   contributors may be used to endorse or promote products derived from
    !   this software without specific prior written permission.
    ! The copyright holders provide no reassurances that the source code
    ! provided does not infringe any patent, copyright, or any other
    ! intellectual property rights of third parties. The copyright holders
    ! disclaim any liability to any recipient for claims brought against
    ! recipient by any third party for infringement of that recipients
    ! intellectual property rights.
    ! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    ! IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    ! TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    ! PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
    ! OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
    ! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
    ! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
    ! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
    ! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
    ! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    ! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

    !> DSTEDC: eigenvalues and optionally eigenvectors of a symmetric
    !> tridiagonal matrix by divide and conquer. COMPZ 'N', 'I' (vectors of
    !> the tridiagonal matrix) or 'V' (Z times them); LAPACK's workspace
    !> contract, though the routine keeps its own scratch.
    subroutine fortnum_dstedc(compz, n, d, e, z, ldz, work, lwork, iwork, liwork, info)
        character, intent(in) :: compz
        integer, intent(in) :: n, ldz, lwork, liwork
        real(dp), intent(inout) :: d(*), e(*), z(ldz, *), work(*)
        integer, intent(inout) :: iwork(*)
        integer, intent(out) :: info
        integer :: icompz, lwmin, liwmin, lgn
        real(dp), allocatable :: v(:, :), zt(:, :)
        logical :: query

        info = 0
        query = lwork == -1 .or. liwork == -1
        icompz = index('NnVvIi', compz)
        if (icompz == 0) then
            info = -1
        else if (n < 0) then
            info = -2
        else if (ldz < 1 .or. (icompz > 2 .and. ldz < max(1, n))) then
            info = -6
        end if
        icompz = (icompz - 1)/2
        if (info == 0) then
            if (n <= 1 .or. icompz == 0) then
                lwmin = 1
                liwmin = 1
            else
                lgn = ceiling(log(real(n, dp))/log(2.0_dp))
                if (2**lgn < n) lgn = lgn + 1
                if (icompz == 1) then
                    lwmin = 1 + 3*n + 2*n*lgn + 4*n**2
                    liwmin = 6 + 6*n + 5*n*lgn
                else
                    lwmin = 1 + 4*n + n**2
                    liwmin = 3 + 5*n
                end if
            end if
            work(1) = lwmin
            iwork(1) = liwmin
            if (lwork < lwmin .and. .not. query) then
                info = -8
            else if (liwork < liwmin .and. .not. query) then
                info = -10
            end if
        end if
        if (info /= 0 .or. query) return
        if (n == 0) return
        if (n == 1) then
            if (icompz /= 0) z(1, 1) = 1
            return
        end if
        select case (icompz)
        case (0)
            call fortnum_dsterf(n, d, e, info)
        case (2)
            call stedc_vectors(n, d, e, z, ldz, info)
        case default
            ! Z times the vectors of the tridiagonal matrix.
            allocate (v(n, n), zt(n, n))
            call stedc_vectors(n, d, e, v, n, info)
            if (info == 0) then
                zt = z(1:n, 1:n)
                call fortnum_dgemm('N', 'N', n, n, n, 1.0_dp, zt, n, v, n, 0.0_dp, z, ldz)
            end if
        end select
        work(1) = lwmin
        iwork(1) = liwmin
    end subroutine fortnum_dstedc

    !> Eigenpairs of the tridiagonal matrix (d, e) into z (COMPZ = 'I'),
    !> ascending; DSTEDC's splitting into unreduced blocks with scaling.
    subroutine stedc_vectors(n, d, e, z, ldz, info)
        integer, intent(in) :: n, ldz
        real(dp), intent(inout) :: d(*), e(*), z(ldz, *)
        integer, intent(out) :: info
        real(dp), allocatable :: work(:)
        integer, allocatable :: iwork(:)
        real(dp) :: orgnrm, eps, tiny_e, p
        integer :: start, finish, m, i, j, k

        info = 0
        z(1:n, 1:n) = 0
        do j = 1, n
            z(j, j) = 1
        end do
        if (n <= STEDC_LEAF) then
            allocate (work(max(1, 2*n - 2)))
            call fortnum_dsteqr('I', n, d, e, z, ldz, work, info)
            return
        end if
        orgnrm = max(maxval(abs(d(1:n))), maxval(abs(e(1:n - 1))))
        if (orgnrm == 0) return
        eps = epsilon(1.0_dp)
        allocate (work(1 + 4*n + n**2), iwork(3 + 5*n))
        start = 1
        do while (start <= n)
            ! The unreduced block start..finish.
            finish = start
            do while (finish < n)
                tiny_e = eps*sqrt(abs(d(finish)))*sqrt(abs(d(finish + 1)))
                if (abs(e(finish)) <= tiny_e) exit
                finish = finish + 1
            end do
            m = finish - start + 1
            if (m > 1) then
                if (m > STEDC_LEAF) then
                    orgnrm = max(maxval(abs(d(start:finish))), &
                        maxval(abs(e(start:finish - 1))))
                    d(start:finish) = d(start:finish)/orgnrm
                    e(start:finish - 1) = e(start:finish - 1)/orgnrm
                    call laed0(m, d(start), e(start), z(start, start), ldz, work, iwork, info)
                    if (info /= 0) then
                        info = (info/(m + 1) + start - 1)*(n + 1) + mod(info, m + 1) + start - 1
                        return
                    end if
                    d(start:finish) = d(start:finish)*orgnrm
                else
                    call fortnum_dsteqr('I', m, d(start), e(start), z(start, start), ldz, &
                        work, info)
                    if (info /= 0) then
                        info = start*(n + 1) + finish
                        return
                    end if
                end if
            end if
            start = finish + 1
        end do
        ! Selection sort of the eigenvalues with their vectors.
        do i = 1, n - 1
            k = i
            p = d(i)
            do j = i + 1, n
                if (d(j) < p) then
                    k = j
                    p = d(j)
                end if
            end do
            if (k /= i) then
                d(k) = d(i)
                d(i) = p
                work(1:n) = z(1:n, i)
                z(1:n, i) = z(1:n, k)
                z(1:n, k) = work(1:n)
            end if
        end do
    end subroutine stedc_vectors

    !> DLAED0 with ICOMPQ = 2: eigenpairs of the unreduced tridiagonal block
    !> into q (set to the identity by the caller).
    subroutine laed0(n, d, e, q, ldq, work, iwork, info)
        integer, intent(in) :: n, ldq
        real(dp), intent(inout) :: d(*), e(*), q(ldq, *), work(*)
        integer, intent(inout) :: iwork(*)
        integer, intent(out) :: info
        integer :: subpbs, tlvls, j, i, spm1, submat, smm1, indxq, matsiz, k
        integer :: spm2, msd2

        info = 0
        ! The subproblem sizes, halving until each is at most STEDC_LEAF.
        iwork(1) = n
        subpbs = 1
        tlvls = 0
        do while (iwork(subpbs) > STEDC_LEAF)
            do j = subpbs, 1, -1
                iwork(2*j) = (iwork(j) + 1)/2
                iwork(2*j - 1) = iwork(j)/2
            end do
            tlvls = tlvls + 1
            subpbs = 2*subpbs
        end do
        do j = 2, subpbs
            iwork(j) = iwork(j) + iwork(j - 1)
        end do
        ! Rank-one modifications that split the blocks.
        spm1 = subpbs - 1
        do i = 1, spm1
            submat = iwork(i) + 1
            smm1 = submat - 1
            d(smm1) = d(smm1) - abs(e(smm1))
            d(submat) = d(submat) - abs(e(smm1))
        end do
        indxq = 4*n + 3
        ! The leaves.
        do i = 0, spm1
            if (i == 0) then
                submat = 1
                matsiz = iwork(1)
            else
                submat = iwork(i) + 1
                matsiz = iwork(i + 1) - iwork(i)
            end if
            call fortnum_dsteqr('I', matsiz, d(submat), e(submat), q(submat, submat), ldq, &
                work, info)
            if (info /= 0) then
                info = submat*(n + 1) + submat + matsiz - 1
                return
            end if
            k = 1
            do j = submat, iwork(i + 1)
                iwork(indxq + j) = k
                k = k + 1
            end do
        end do
        ! Merge pairs of subproblems up the tree.
        do while (subpbs > 1)
            spm2 = subpbs - 2
            do i = 0, spm2, 2
                if (i == 0) then
                    submat = 1
                    matsiz = iwork(2)
                    msd2 = iwork(1)
                else
                    submat = iwork(i) + 1
                    matsiz = iwork(i + 2) - iwork(i)
                    msd2 = matsiz/2
                end if
                call laed1(matsiz, d(submat), q(submat, submat), ldq, iwork(indxq + submat), &
                    e(submat + msd2 - 1), msd2, work, iwork(subpbs + 1), info)
                if (info /= 0) then
                    info = submat*(n + 1) + submat + matsiz - 1
                    return
                end if
                iwork(i/2 + 1) = iwork(i + 2)
            end do
            subpbs = subpbs/2
        end do
        ! Sort the eigenvalues by the merge permutation.
        do i = 1, n
            j = iwork(indxq + i)
            work(i) = d(j)
            work(n*i + 1:n*i + n) = q(1:n, j)
        end do
        d(1:n) = work(1:n)
        do i = 1, n
            q(1:n, i) = work(n*i + 1:n*i + n)
        end do
    end subroutine laed0

    !> DLAED1: eigenpairs of the merged block from those of its halves and
    !> the rank-one modification rho z z^T, z built from q.
    subroutine laed1(n, d, q, ldq, indxq, rho, cutpnt, work, iwork, info)
        integer, intent(in) :: n, ldq, cutpnt
        real(dp), intent(inout) :: d(*), q(ldq, *), work(*), rho
        integer, intent(inout) :: indxq(*), iwork(*)
        integer, intent(out) :: info
        integer :: iz, idlmda, iw, iq2, indx, indxc, coltyp, indxp, k, is, zpp1, i
        real(dp) :: rho_local

        info = 0
        if (n == 0) return
        iz = 1
        idlmda = iz + n
        iw = idlmda + n
        iq2 = iw + n
        indx = 1
        indxc = indx + n
        coltyp = indxc + n
        indxp = coltyp + n
        ! z: the last row of the first eigenvector block and the first row of
        ! the second.
        do i = 1, cutpnt
            work(iz + i - 1) = q(cutpnt, i)
        end do
        zpp1 = cutpnt + 1
        do i = 1, n - cutpnt
            work(iz + cutpnt + i - 1) = q(zpp1, zpp1 + i - 1)
        end do
        rho_local = rho
        call laed2(k, n, cutpnt, d, q, ldq, indxq, rho_local, work(iz), work(idlmda), &
            work(iw), work(iq2), iwork(indx), iwork(indxc), iwork(indxp), iwork(coltyp), info)
        if (info /= 0) return
        if (k /= 0) then
            is = (iwork(coltyp) + iwork(coltyp + 1))*cutpnt + &
                (iwork(coltyp + 1) + iwork(coltyp + 2))*(n - cutpnt) + iq2
            call laed3(k, n, cutpnt, d, q, ldq, rho_local, work(idlmda), work(iq2), &
                iwork(indxc), iwork(coltyp), work(iw), work(is), info)
            if (info /= 0) return
            call lamrg(k, n - k, d, 1, -1, indxq)
        else
            do i = 1, n
                indxq(i) = i
            end do
        end if
    end subroutine laed1

    !> DLAED2: deflation of the merged problem; k is the size of the secular
    !> equation that remains.
    subroutine laed2(k, n, n1, d, q, ldq, indxq, rho, z, dlambda, w, q2, indx, indxc, &
            indxp, coltyp, info)
        integer, intent(out) :: k, info
        integer, intent(in) :: n, n1, ldq
        real(dp), intent(inout) :: d(*), q(ldq, *), rho, z(*), dlambda(*), w(*), q2(*)
        integer, intent(inout) :: indxq(*), indx(*), indxc(*), indxp(*), coltyp(*)
        integer :: ctot(4), psm(4), ct, i, imax, iq1, iq2, j, jmax, js, k2, n1p1, n2, nj, pj
        real(dp) :: c, s, t, tau, tol

        info = 0
        k = 0
        if (n == 0) return
        n2 = n - n1
        n1p1 = n1 + 1
        if (rho < 0) z(n1p1:n) = -z(n1p1:n)
        ! z is the concatenation of two unit vectors: normalize it.
        z(1:n) = z(1:n)/sqrt(2.0_dp)
        rho = abs(2*rho)
        do i = n1p1, n
            indxq(i) = indxq(i) + n1
        end do
        do i = 1, n
            dlambda(i) = d(indxq(i))
        end do
        call lamrg(n1, n2, dlambda, 1, 1, indxc)
        do i = 1, n
            indx(i) = indxq(indxc(i))
        end do
        imax = maxloc(abs(z(1:n)), 1)
        jmax = maxloc(abs(d(1:n)), 1)
        tol = 8*epsilon(1.0_dp)*max(abs(d(jmax)), abs(z(imax)))
        if (rho*abs(z(imax)) <= tol) then
            ! Everything deflates: sort and return.
            k = 0
            iq2 = 1
            do j = 1, n
                i = indx(j)
                q2(iq2:iq2 + n - 1) = q(1:n, i)
                dlambda(j) = d(i)
                iq2 = iq2 + n
            end do
            do j = 1, n
                q(1:n, j) = q2((j - 1)*n + 1:j*n)
            end do
            d(1:n) = dlambda(1:n)
            return
        end if
        coltyp(1:n1) = 1
        coltyp(n1p1:n) = 3
        k = 0
        k2 = n + 1
        ! The first component that does not deflate.
        pj = 0
        do j = 1, n
            nj = indx(j)
            if (rho*abs(z(nj)) <= tol) then
                k2 = k2 - 1
                coltyp(nj) = 4
                indxp(k2) = nj
            else
                pj = nj
                exit
            end if
        end do
        if (pj /= 0) then
            do
                j = j + 1
                if (j > n) exit
                nj = indx(j)
                if (rho*abs(z(nj)) <= tol) then
                    k2 = k2 - 1
                    coltyp(nj) = 4
                    indxp(k2) = nj
                else
                    s = z(pj)
                    c = z(nj)
                    tau = hypot(c, s)
                    t = d(nj) - d(pj)
                    c = c/tau
                    s = -s/tau
                    if (abs(t*c*s) <= tol) then
                        ! Close eigenvalues: rotate one component of z away.
                        z(nj) = tau
                        z(pj) = 0
                        if (coltyp(nj) /= coltyp(pj)) coltyp(nj) = 2
                        coltyp(pj) = 4
                        call rotate(n, q(1, pj), q(1, nj), c, s)
                        t = d(pj)*c**2 + d(nj)*s**2
                        d(nj) = d(pj)*s**2 + d(nj)*c**2
                        d(pj) = t
                        k2 = k2 - 1
                        i = 1
                        do
                            if (k2 + i <= n) then
                                if (d(pj) < d(indxp(k2 + i))) then
                                    indxp(k2 + i - 1) = indxp(k2 + i)
                                    indxp(k2 + i) = pj
                                    i = i + 1
                                    cycle
                                end if
                            end if
                            indxp(k2 + i - 1) = pj
                            exit
                        end do
                        pj = nj
                    else
                        k = k + 1
                        dlambda(k) = d(pj)
                        w(k) = z(pj)
                        indxp(k) = pj
                        pj = nj
                    end if
                end if
            end do
            k = k + 1
            dlambda(k) = d(pj)
            w(k) = z(pj)
            indxp(k) = pj
        end if
        ! Columns by type: 1 first half only, 2 both, 3 second half, 4 deflated.
        ctot = 0
        do j = 1, n
            ct = coltyp(j)
            ctot(ct) = ctot(ct) + 1
        end do
        psm(1) = 1
        psm(2) = 1 + ctot(1)
        psm(3) = psm(2) + ctot(2)
        psm(4) = psm(3) + ctot(3)
        k = n - ctot(4)
        do j = 1, n
            js = indxp(j)
            ct = coltyp(js)
            indx(psm(ct)) = js
            indxc(psm(ct)) = j
            psm(ct) = psm(ct) + 1
        end do
        i = 1
        iq1 = 1
        iq2 = 1 + (ctot(1) + ctot(2))*n1
        do j = 1, ctot(1)
            js = indx(i)
            q2(iq1:iq1 + n1 - 1) = q(1:n1, js)
            z(i) = d(js)
            i = i + 1
            iq1 = iq1 + n1
        end do
        do j = 1, ctot(2)
            js = indx(i)
            q2(iq1:iq1 + n1 - 1) = q(1:n1, js)
            q2(iq2:iq2 + n2 - 1) = q(n1 + 1:n, js)
            z(i) = d(js)
            i = i + 1
            iq1 = iq1 + n1
            iq2 = iq2 + n2
        end do
        do j = 1, ctot(3)
            js = indx(i)
            q2(iq2:iq2 + n2 - 1) = q(n1 + 1:n, js)
            z(i) = d(js)
            i = i + 1
            iq2 = iq2 + n2
        end do
        iq1 = iq2
        do j = 1, ctot(4)
            js = indx(i)
            q2(iq2:iq2 + n - 1) = q(1:n, js)
            iq2 = iq2 + n
            z(i) = d(js)
            i = i + 1
        end do
        if (k < n) then
            do j = 1, ctot(4)
                q(1:n, k + j) = q2(iq1 + (j - 1)*n:iq1 + j*n - 1)
            end do
            d(k + 1:n) = z(k + 1:n)
        end if
        coltyp(1:4) = ctot
    end subroutine laed2

    !> DLAED3: roots of the secular equation and the eigenvectors of the
    !> merged block, by GEMM with the stored halves.
    subroutine laed3(k, n, n1, d, q, ldq, rho, dlambda, q2, indx, ctot, w, s, info)
        integer, intent(in) :: k, n, n1, ldq
        real(dp), intent(inout) :: d(*), q(ldq, *), dlambda(*), q2(*), w(*), s(*)
        real(dp), intent(in) :: rho
        integer, intent(in) :: indx(*), ctot(*)
        integer, intent(out) :: info
        integer :: i, j, ii, iq2, n12, n2, n23
        real(dp) :: temp

        info = 0
        if (k == 0) return
        do j = 1, k
            call laed4(k, j, dlambda, w, q(1, j), rho, d(j), info)
            if (info /= 0) return
        end do
        if (k == 2) then
            do j = 1, k
                w(1) = q(1, j)
                w(2) = q(2, j)
                q(1, j) = w(indx(1))
                q(2, j) = w(indx(2))
            end do
        else if (k > 2) then
            ! Recompute w so that the eigenvectors are orthogonal (Gu and
            ! Eisenstat).
            s(1:k) = w(1:k)
            do i = 1, k
                w(i) = q(i, i)
            end do
            do j = 1, k
                do i = 1, j - 1
                    w(i) = w(i)*(q(i, j)/(dlambda(i) - dlambda(j)))
                end do
                do i = j + 1, k
                    w(i) = w(i)*(q(i, j)/(dlambda(i) - dlambda(j)))
                end do
            end do
            do i = 1, k
                w(i) = sign(sqrt(-w(i)), s(i))
            end do
            do j = 1, k
                do i = 1, k
                    s(i) = w(i)/q(i, j)
                end do
                temp = norm2(s(1:k))
                do i = 1, k
                    ii = indx(i)
                    q(i, j) = s(ii)/temp
                end do
            end do
        end if
        ! The eigenvectors of the merged block: stored halves times the
        ! secular-equation vectors.
        n2 = n - n1
        n12 = ctot(1) + ctot(2)
        n23 = ctot(2) + ctot(3)
        do j = 1, k
            s((j - 1)*n23 + 1:j*n23) = q(ctot(1) + 1:ctot(1) + n23, j)
        end do
        iq2 = n1*n12 + 1
        if (n23 /= 0) then
            call fortnum_dgemm('N', 'N', n2, k, n23, 1.0_dp, q2(iq2), n2, s, n23, 0.0_dp, &
                q(n1 + 1, 1), ldq)
        else
            q(n1 + 1:n, 1:k) = 0
        end if
        do j = 1, k
            s((j - 1)*n12 + 1:j*n12) = q(1:n12, j)
        end do
        if (n12 /= 0) then
            call fortnum_dgemm('N', 'N', n1, k, n12, 1.0_dp, q2, n1, s, n12, 0.0_dp, q, ldq)
        else
            q(1:n1, 1:k) = 0
        end if
    end subroutine laed3

    !> DLAED4: the i-th root dlam of the secular equation
    !> 1 + rho sum_j z(j)^2/(d(j) - lambda) = 0 for d ascending, and
    !> delta(j) = d(j) - dlam.
    subroutine laed4(n, i, d, z, delta, rho, dlam, info)
        integer, intent(in) :: n, i
        real(dp), intent(in) :: d(*), z(*), rho
        real(dp), intent(out) :: delta(*), dlam
        integer, intent(out) :: info
        integer, parameter :: maxit = 30
        logical :: orgati, swtch, swtch3
        integer :: ii, iim1, iip1, ip1, iter, j, niter
        real(dp) :: a, b, c, del, dltlb, dltub, dphi, dpsi, dw, eps, erretm, eta, midpt, &
            phi, prew, psi, rhoinv, tau, temp, temp1, w, zz(3)

        info = 0
        if (n == 1) then
            dlam = d(1) + rho*z(1)*z(1)
            delta(1) = 1
            return
        end if
        if (n == 2) then
            call laed5(i, d, z, delta, rho, dlam)
            return
        end if
        eps = epsilon(1.0_dp)/2
        rhoinv = 1/rho
        if (i == n) then
            ! The largest root, between d(n) and d(n) + rho.
            ii = n - 1
            niter = 1
            midpt = rho/2
            do j = 1, n
                delta(j) = (d(j) - d(i)) - midpt
            end do
            psi = 0
            do j = 1, n - 2
                psi = psi + z(j)*z(j)/delta(j)
            end do
            c = rhoinv + psi
            w = c + z(ii)*z(ii)/delta(ii) + z(n)*z(n)/delta(n)
            if (w <= 0) then
                temp = z(n - 1)*z(n - 1)/(d(n) - d(n - 1) + rho) + z(n)*z(n)/rho
                if (c <= temp) then
                    tau = rho
                else
                    del = d(n) - d(n - 1)
                    a = -c*del + z(n - 1)*z(n - 1) + z(n)*z(n)
                    b = z(n)*z(n)*del
                    if (a < 0) then
                        tau = 2*b/(sqrt(a*a + 4*b*c) - a)
                    else
                        tau = (a + sqrt(a*a + 4*b*c))/(2*c)
                    end if
                end if
                dltlb = midpt
                dltub = rho
            else
                del = d(n) - d(n - 1)
                a = -c*del + z(n - 1)*z(n - 1) + z(n)*z(n)
                b = z(n)*z(n)*del
                if (a < 0) then
                    tau = 2*b/(sqrt(a*a + 4*b*c) - a)
                else
                    tau = (a + sqrt(a*a + 4*b*c))/(2*c)
                end if
                dltlb = 0
                dltub = midpt
            end if
            do j = 1, n
                delta(j) = (d(j) - d(i)) - tau
            end do
            call last_sums()
            if (abs(w) <= eps*erretm) then
                dlam = d(i) + tau
                return
            end if
            if (w <= 0) then
                dltlb = max(dltlb, tau)
            else
                dltub = min(dltub, tau)
            end if
            niter = niter + 1
            c = w - delta(n - 1)*dpsi - delta(n)*dphi
            a = (delta(n - 1) + delta(n))*w - delta(n - 1)*delta(n)*(dpsi + dphi)
            b = delta(n - 1)*delta(n)*w
            if (c < 0) c = abs(c)
            if (c == 0) then
                eta = -w/(dpsi + dphi)
            else if (a >= 0) then
                eta = (a + sqrt(abs(a*a - 4*b*c)))/(2*c)
            else
                eta = 2*b/(a - sqrt(abs(a*a - 4*b*c)))
            end if
            if (w*eta > 0) eta = -w/(dpsi + dphi)
            temp = tau + eta
            if (temp > dltub .or. temp < dltlb) then
                if (w < 0) then
                    eta = (dltub - tau)/2
                else
                    eta = (dltlb - tau)/2
                end if
            end if
            delta(1:n) = delta(1:n) - eta
            tau = tau + eta
            call last_sums()
            iter = niter + 1
            do niter = iter, maxit
                if (abs(w) <= eps*erretm) then
                    dlam = d(i) + tau
                    return
                end if
                if (w <= 0) then
                    dltlb = max(dltlb, tau)
                else
                    dltub = min(dltub, tau)
                end if
                c = w - delta(n - 1)*dpsi - delta(n)*dphi
                a = (delta(n - 1) + delta(n))*w - delta(n - 1)*delta(n)*(dpsi + dphi)
                b = delta(n - 1)*delta(n)*w
                if (a >= 0) then
                    eta = (a + sqrt(abs(a*a - 4*b*c)))/(2*c)
                else
                    eta = 2*b/(a - sqrt(abs(a*a - 4*b*c)))
                end if
                if (w*eta > 0) eta = -w/(dpsi + dphi)
                temp = tau + eta
                if (temp > dltub .or. temp < dltlb) then
                    if (w < 0) then
                        eta = (dltub - tau)/2
                    else
                        eta = (dltlb - tau)/2
                    end if
                end if
                delta(1:n) = delta(1:n) - eta
                tau = tau + eta
                call last_sums()
            end do
            info = 1
            dlam = d(i) + tau
            return
        end if
        ! A root between d(i) and d(i + 1).
        niter = 1
        ip1 = i + 1
        del = d(ip1) - d(i)
        midpt = del/2
        do j = 1, n
            delta(j) = (d(j) - d(i)) - midpt
        end do
        psi = 0
        do j = 1, i - 1
            psi = psi + z(j)*z(j)/delta(j)
        end do
        phi = 0
        do j = n, i + 2, -1
            phi = phi + z(j)*z(j)/delta(j)
        end do
        c = rhoinv + psi + phi
        w = c + z(i)*z(i)/delta(i) + z(ip1)*z(ip1)/delta(ip1)
        if (w > 0) then
            ! The root is nearer d(i).
            orgati = .true.
            a = c*del + z(i)*z(i) + z(ip1)*z(ip1)
            b = z(i)*z(i)*del
            if (a > 0) then
                tau = 2*b/(a + sqrt(abs(a*a - 4*b*c)))
            else
                tau = (a - sqrt(abs(a*a - 4*b*c)))/(2*c)
            end if
            dltlb = 0
            dltub = midpt
        else
            orgati = .false.
            a = c*del - z(i)*z(i) - z(ip1)*z(ip1)
            b = z(ip1)*z(ip1)*del
            if (a < 0) then
                tau = 2*b/(a - sqrt(abs(a*a + 4*b*c)))
            else
                tau = -(a + sqrt(abs(a*a + 4*b*c)))/(2*c)
            end if
            dltlb = -midpt
            dltub = 0
        end if
        if (orgati) then
            do j = 1, n
                delta(j) = (d(j) - d(i)) - tau
            end do
            ii = i
        else
            do j = 1, n
                delta(j) = (d(j) - d(ip1)) - tau
            end do
            ii = i + 1
        end if
        iim1 = ii - 1
        iip1 = ii + 1
        call middle_sums()
        w = rhoinv + phi + psi
        swtch3 = .false.
        if (orgati) then
            if (w < 0) swtch3 = .true.
        else
            if (w > 0) swtch3 = .true.
        end if
        if (ii == 1 .or. ii == n) swtch3 = .false.
        temp = z(ii)/delta(ii)
        dw = dpsi + dphi + temp*temp
        temp = z(ii)*temp
        w = w + temp
        erretm = 8*(phi - psi) + erretm + 2*rhoinv + 3*abs(temp) + abs(tau)*dw
        if (abs(w) <= eps*erretm) then
            call finish_root()
            return
        end if
        if (w <= 0) then
            dltlb = max(dltlb, tau)
        else
            dltub = min(dltub, tau)
        end if
        niter = niter + 1
        if (.not. swtch3) then
            if (orgati) then
                c = w - delta(ip1)*dw - (d(i) - d(ip1))*(z(i)/delta(i))**2
            else
                c = w - delta(i)*dw - (d(ip1) - d(i))*(z(ip1)/delta(ip1))**2
            end if
            a = (delta(i) + delta(ip1))*w - delta(i)*delta(ip1)*dw
            b = delta(i)*delta(ip1)*w
            if (c == 0) then
                if (a == 0) then
                    if (orgati) then
                        a = z(i)*z(i) + delta(ip1)*delta(ip1)*(dpsi + dphi)
                    else
                        a = z(ip1)*z(ip1) + delta(i)*delta(i)*(dpsi + dphi)
                    end if
                end if
                eta = b/a
            else if (a <= 0) then
                eta = (a - sqrt(abs(a*a - 4*b*c)))/(2*c)
            else
                eta = 2*b/(a + sqrt(abs(a*a - 4*b*c)))
            end if
        else
            ! Interpolation with three poles.
            temp = rhoinv + psi + phi
            if (orgati) then
                temp1 = z(iim1)/delta(iim1)
                temp1 = temp1*temp1
                c = temp - delta(iip1)*(dpsi + dphi) - (d(iim1) - d(iip1))*temp1
                zz(1) = z(iim1)*z(iim1)
                zz(3) = delta(iip1)*delta(iip1)*((dpsi - temp1) + dphi)
            else
                temp1 = z(iip1)/delta(iip1)
                temp1 = temp1*temp1
                c = temp - delta(iim1)*(dpsi + dphi) - (d(iip1) - d(iim1))*temp1
                zz(1) = delta(iim1)*delta(iim1)*(dpsi + (dphi - temp1))
                zz(3) = z(iip1)*z(iip1)
            end if
            zz(2) = z(ii)*z(ii)
            call laed6(niter, orgati, c, delta(iim1), zz, w, eta, info)
            if (info /= 0) return
        end if
        if (w*eta >= 0) eta = -w/dw
        temp = tau + eta
        if (temp > dltub .or. temp < dltlb) then
            if (w < 0) then
                eta = (dltub - tau)/2
            else
                eta = (dltlb - tau)/2
            end if
        end if
        prew = w
        delta(1:n) = delta(1:n) - eta
        call middle_sums()
        temp = z(ii)/delta(ii)
        dw = dpsi + dphi + temp*temp
        temp = z(ii)*temp
        w = rhoinv + phi + psi + temp
        erretm = 8*(phi - psi) + erretm + 2*rhoinv + 3*abs(temp) + abs(tau + eta)*dw
        swtch = .false.
        if (orgati) then
            if (-w > abs(prew)/10) swtch = .true.
        else
            if (w > abs(prew)/10) swtch = .true.
        end if
        tau = tau + eta
        iter = niter + 1
        do niter = iter, maxit
            if (abs(w) <= eps*erretm) then
                call finish_root()
                return
            end if
            if (w <= 0) then
                dltlb = max(dltlb, tau)
            else
                dltub = min(dltub, tau)
            end if
            if (.not. swtch3) then
                if (.not. swtch) then
                    if (orgati) then
                        c = w - delta(ip1)*dw - (d(i) - d(ip1))*(z(i)/delta(i))**2
                    else
                        c = w - delta(i)*dw - (d(ip1) - d(i))*(z(ip1)/delta(ip1))**2
                    end if
                else
                    temp = z(ii)/delta(ii)
                    if (orgati) then
                        dpsi = dpsi + temp*temp
                    else
                        dphi = dphi + temp*temp
                    end if
                    c = w - delta(i)*dpsi - delta(ip1)*dphi
                end if
                a = (delta(i) + delta(ip1))*w - delta(i)*delta(ip1)*dw
                b = delta(i)*delta(ip1)*w
                if (c == 0) then
                    if (a == 0) then
                        if (.not. swtch) then
                            if (orgati) then
                                a = z(i)*z(i) + delta(ip1)*delta(ip1)*(dpsi + dphi)
                            else
                                a = z(ip1)*z(ip1) + delta(i)*delta(i)*(dpsi + dphi)
                            end if
                        else
                            a = delta(i)*delta(i)*dpsi + delta(ip1)*delta(ip1)*dphi
                        end if
                    end if
                    eta = b/a
                else if (a <= 0) then
                    eta = (a - sqrt(abs(a*a - 4*b*c)))/(2*c)
                else
                    eta = 2*b/(a + sqrt(abs(a*a - 4*b*c)))
                end if
            else
                temp = rhoinv + psi + phi
                if (swtch) then
                    c = temp - delta(iim1)*dpsi - delta(iip1)*dphi
                    zz(1) = delta(iim1)*delta(iim1)*dpsi
                    zz(3) = delta(iip1)*delta(iip1)*dphi
                else
                    if (orgati) then
                        temp1 = z(iim1)/delta(iim1)
                        temp1 = temp1*temp1
                        c = temp - delta(iip1)*(dpsi + dphi) - (d(iim1) - d(iip1))*temp1
                        zz(1) = z(iim1)*z(iim1)
                        zz(3) = delta(iip1)*delta(iip1)*((dpsi - temp1) + dphi)
                    else
                        temp1 = z(iip1)/delta(iip1)
                        temp1 = temp1*temp1
                        c = temp - delta(iim1)*(dpsi + dphi) - (d(iip1) - d(iim1))*temp1
                        zz(1) = delta(iim1)*delta(iim1)*(dpsi + (dphi - temp1))
                        zz(3) = z(iip1)*z(iip1)
                    end if
                end if
                call laed6(niter, orgati, c, delta(iim1), zz, w, eta, info)
                if (info /= 0) return
            end if
            if (w*eta >= 0) eta = -w/dw
            temp = tau + eta
            if (temp > dltub .or. temp < dltlb) then
                if (w < 0) then
                    eta = (dltub - tau)/2
                else
                    eta = (dltlb - tau)/2
                end if
            end if
            delta(1:n) = delta(1:n) - eta
            tau = tau + eta
            prew = w
            call middle_sums()
            temp = z(ii)/delta(ii)
            dw = dpsi + dphi + temp*temp
            temp = z(ii)*temp
            w = rhoinv + phi + psi + temp
            erretm = 8*(phi - psi) + erretm + 2*rhoinv + 3*abs(temp) + abs(tau)*dw
            if (w*prew > 0 .and. abs(w) > abs(prew)/10) swtch = .not. swtch
        end do
        info = 1
        call finish_root()
    contains
        subroutine last_sums()
            ! psi, phi, their derivatives and the error bound for i = n.
            integer :: jj
            real(dp) :: t
            dpsi = 0
            psi = 0
            erretm = 0
            do jj = 1, ii
                t = z(jj)/delta(jj)
                psi = psi + z(jj)*t
                dpsi = dpsi + t*t
                erretm = erretm + psi
            end do
            erretm = abs(erretm)
            t = z(n)/delta(n)
            phi = z(n)*t
            dphi = t*t
            erretm = 8*(-phi - psi) + erretm - phi + rhoinv + abs(tau)*(dpsi + dphi)
            w = rhoinv + phi + psi
        end subroutine last_sums

        subroutine middle_sums()
            ! psi over 1..iim1 and phi over iip1..n, with derivatives.
            integer :: jj
            real(dp) :: t
            dpsi = 0
            psi = 0
            erretm = 0
            do jj = 1, iim1
                t = z(jj)/delta(jj)
                psi = psi + z(jj)*t
                dpsi = dpsi + t*t
                erretm = erretm + psi
            end do
            erretm = abs(erretm)
            dphi = 0
            phi = 0
            do jj = n, iip1, -1
                t = z(jj)/delta(jj)
                phi = phi + z(jj)*t
                dphi = dphi + t*t
                erretm = erretm + phi
            end do
        end subroutine middle_sums

        subroutine finish_root()
            if (orgati) then
                dlam = d(i) + tau
            else
                dlam = d(ip1) + tau
            end if
        end subroutine finish_root
    end subroutine laed4

    !> DLAED5: the i-th eigenvalue of diag(d) + rho z z^T for n = 2.
    subroutine laed5(i, d, z, delta, rho, dlam)
        integer, intent(in) :: i
        real(dp), intent(in) :: d(2), z(2), rho
        real(dp), intent(out) :: delta(2), dlam
        real(dp) :: b, c, del, tau, temp, w

        del = d(2) - d(1)
        if (i == 1) then
            w = 1 + 2*rho*(z(2)*z(2) - z(1)*z(1))/del
            if (w > 0) then
                b = del + rho*(z(1)*z(1) + z(2)*z(2))
                c = rho*z(1)*z(1)*del
                tau = 2*c/(b + sqrt(abs(b*b - 4*c)))
                dlam = d(1) + tau
                delta(1) = -z(1)/tau
                delta(2) = z(2)/(del - tau)
            else
                b = -del + rho*(z(1)*z(1) + z(2)*z(2))
                c = rho*z(2)*z(2)*del
                if (b > 0) then
                    tau = -2*c/(b + sqrt(b*b + 4*c))
                else
                    tau = (b - sqrt(b*b + 4*c))/2
                end if
                dlam = d(2) + tau
                delta(1) = -z(1)/(del + tau)
                delta(2) = -z(2)/tau
            end if
        else
            b = -del + rho*(z(1)*z(1) + z(2)*z(2))
            c = rho*z(2)*z(2)*del
            if (b > 0) then
                tau = (b + sqrt(b*b + 4*c))/2
            else
                tau = 2*c/(-b + sqrt(b*b + 4*c))
            end if
            dlam = d(2) + tau
            delta(1) = -z(1)/(del + tau)
            delta(2) = -z(2)/tau
        end if
        temp = sqrt(delta(1)*delta(1) + delta(2)*delta(2))
        delta = delta/temp
    end subroutine laed5

    !> DLAED6: one step of the root of the three-pole secular function
    !> finit + tau*sum_j z(j)/(d(j)*(d(j) - tau)) near the origin.
    subroutine laed6(kniter, orgati, rho, d, z, finit, tau, info)
        integer, intent(in) :: kniter
        logical, intent(in) :: orgati
        real(dp), intent(in) :: rho, d(3), z(3), finit
        real(dp), intent(out) :: tau
        integer, intent(out) :: info
        integer, parameter :: maxit = 40
        real(dp) :: a, b, c, ddf, df, eps, erretm, eta, f, fc, sclfac, sclinv, small1, &
            small2, sminv1, sminv2, temp, temp1, temp2, temp3, temp4, lbd, ubd, &
            dscale(3), zscale(3)
        logical :: scale
        integer :: j, niter

        info = 0
        if (orgati) then
            lbd = d(2)
            ubd = d(3)
        else
            lbd = d(1)
            ubd = d(2)
        end if
        if (finit < 0) then
            lbd = 0
        else
            ubd = 0
        end if
        niter = 1
        tau = 0
        if (kniter == 2) then
            if (orgati) then
                temp = (d(3) - d(2))/2
                c = rho + z(1)/((d(1) - d(2)) - temp)
                a = c*(d(2) + d(3)) + z(2) + z(3)
                b = c*d(2)*d(3) + z(2)*d(3) + z(3)*d(2)
            else
                temp = (d(1) - d(2))/2
                c = rho + z(3)/((d(3) - d(2)) - temp)
                a = c*(d(1) + d(2)) + z(1) + z(2)
                b = c*d(1)*d(2) + z(1)*d(2) + z(2)*d(1)
            end if
            temp = max(abs(a), abs(b), abs(c))
            a = a/temp
            b = b/temp
            c = c/temp
            if (c == 0) then
                tau = b/a
            else if (a <= 0) then
                tau = (a - sqrt(abs(a*a - 4*b*c)))/(2*c)
            else
                tau = 2*b/(a + sqrt(abs(a*a - 4*b*c)))
            end if
            if (tau < lbd .or. tau > ubd) tau = (lbd + ubd)/2
            if (d(1) == tau .or. d(2) == tau .or. d(3) == tau) then
                tau = 0
            else
                temp = finit + tau*z(1)/(d(1)*(d(1) - tau)) + tau*z(2)/(d(2)*(d(2) - tau)) + &
                    tau*z(3)/(d(3)*(d(3) - tau))
                if (temp <= 0) then
                    lbd = tau
                else
                    ubd = tau
                end if
                if (abs(finit) <= abs(temp)) tau = 0
            end if
        end if
        eps = epsilon(1.0_dp)/2
        small1 = radix(1.0_dp)**(int(log(tiny(1.0_dp))/log(real(radix(1.0_dp), dp))/3))
        sminv1 = 1/small1
        small2 = small1*small1
        sminv2 = sminv1*sminv1
        ! Scale when tau is near a pole.
        if (orgati) then
            temp = min(abs(d(2) - tau), abs(d(3) - tau))
        else
            temp = min(abs(d(1) - tau), abs(d(2) - tau))
        end if
        scale = temp <= small1
        if (scale) then
            if (temp <= small2) then
                sclfac = sminv2
                sclinv = small2
            else
                sclfac = sminv1
                sclinv = small1
            end if
            dscale = d*sclfac
            zscale = z*sclfac
            tau = tau*sclfac
            lbd = lbd*sclfac
            ubd = ubd*sclfac
        else
            dscale = d
            zscale = z
        end if
        fc = 0
        df = 0
        ddf = 0
        do j = 1, 3
            temp = 1/(dscale(j) - tau)
            temp1 = zscale(j)*temp
            temp2 = temp1*temp
            temp3 = temp2*temp
            fc = fc + temp1/dscale(j)
            df = df + temp2
            ddf = ddf + temp3
        end do
        f = finit + tau*fc
        converge: block
            if (abs(f) <= 0) exit converge
            if (f <= 0) then
                lbd = tau
            else
                ubd = tau
            end if
            do niter = niter + 1, maxit
                if (orgati) then
                    temp1 = dscale(2) - tau
                    temp2 = dscale(3) - tau
                else
                    temp1 = dscale(1) - tau
                    temp2 = dscale(2) - tau
                end if
                a = (temp1 + temp2)*f - temp1*temp2*df
                b = temp1*temp2*f
                c = f - (temp1 + temp2)*df + temp1*temp2*ddf
                temp = max(abs(a), abs(b), abs(c))
                a = a/temp
                b = b/temp
                c = c/temp
                if (c == 0) then
                    eta = b/a
                else if (a <= 0) then
                    eta = (a - sqrt(abs(a*a - 4*b*c)))/(2*c)
                else
                    eta = 2*b/(a + sqrt(abs(a*a - 4*b*c)))
                end if
                if (f*eta >= 0) eta = -f/df
                tau = tau + eta
                if (tau < lbd .or. tau > ubd) tau = (lbd + ubd)/2
                fc = 0
                erretm = 0
                df = 0
                ddf = 0
                do j = 1, 3
                    if (dscale(j) - tau == 0) exit converge
                    temp = 1/(dscale(j) - tau)
                    temp1 = zscale(j)*temp
                    temp2 = temp1*temp
                    temp3 = temp2*temp
                    temp4 = temp1/dscale(j)
                    fc = fc + temp4
                    erretm = erretm + abs(temp4)
                    df = df + temp2
                    ddf = ddf + temp3
                end do
                f = finit + tau*fc
                erretm = 8*(abs(finit) + abs(tau)*erretm) + abs(tau)*df
                if (abs(f) <= 4*eps*erretm .or. (ubd - lbd) <= 4*eps*abs(tau)) exit converge
                if (f <= 0) then
                    lbd = tau
                else
                    ubd = tau
                end if
            end do
            info = 1
        end block converge
        if (scale) tau = tau*sclinv
    end subroutine laed6

    !> DLAMRG: the permutation that merges two sorted runs of a into one
    !> ascending list; strides dtrd1 and dtrd2 are +1 or -1.
    subroutine lamrg(n1, n2, a, dtrd1, dtrd2, index)
        integer, intent(in) :: n1, n2, dtrd1, dtrd2
        real(dp), intent(in) :: a(*)
        integer, intent(out) :: index(*)
        integer :: i, ind1, ind2, n1sv, n2sv

        n1sv = n1
        n2sv = n2
        if (dtrd1 > 0) then
            ind1 = 1
        else
            ind1 = n1
        end if
        if (dtrd2 > 0) then
            ind2 = 1 + n1
        else
            ind2 = n1 + n2
        end if
        i = 1
        do while (n1sv > 0 .and. n2sv > 0)
            if (a(ind1) <= a(ind2)) then
                index(i) = ind1
                ind1 = ind1 + dtrd1
                n1sv = n1sv - 1
            else
                index(i) = ind2
                ind2 = ind2 + dtrd2
                n2sv = n2sv - 1
            end if
            i = i + 1
        end do
        do while (n2sv > 0)
            index(i) = ind2
            i = i + 1
            ind2 = ind2 + dtrd2
            n2sv = n2sv - 1
        end do
        do while (n1sv > 0)
            index(i) = ind1
            i = i + 1
            ind1 = ind1 + dtrd1
            n1sv = n1sv - 1
        end do
    end subroutine lamrg

    !> DROT on two columns.
    subroutine rotate(n, x, y, c, s)
        integer, intent(in) :: n
        real(dp), intent(inout) :: x(*), y(*)
        real(dp), intent(in) :: c, s
        real(dp) :: t
        integer :: j
        do j = 1, n
            t = c*x(j) + s*y(j)
            y(j) = c*y(j) - s*x(j)
            x(j) = t
        end do
    end subroutine rotate
