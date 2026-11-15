module fortnum_blas12
    ! Level-1 and level-2 BLAS used by our codes and by fortnum's LAPACK
    ! layers, with the arguments and semantics of reference BLAS, including
    ! increments of either sign (a negative increment walks the vector from
    ! its far end). These routines are memory bound; they are written so the
    ! compiler vectorises the unit-stride paths.
    use fortnum_kinds, only: dp, sp
    implicit none
    private

    public :: fortnum_ddot, fortnum_sdot, fortnum_zdotc, fortnum_zdotu
    public :: fortnum_daxpy, fortnum_zaxpy, fortnum_dscal, fortnum_zscal, fortnum_zdscal
    public :: fortnum_dcopy, fortnum_zcopy, fortnum_dswap, fortnum_zswap
    public :: fortnum_dnrm2, fortnum_dznrm2, fortnum_idamax, fortnum_izamax, fortnum_drot
    public :: fortnum_dasum, fortnum_dzasum
    public :: fortnum_dgemv, fortnum_zgemv, fortnum_dsymv, fortnum_zhemv
    public :: fortnum_dtrsv, fortnum_ztrsv, fortnum_dger, fortnum_zgeru, fortnum_zgerc

contains

    ! Position of element i (1-based) of a vector of n elements with increment inc.
    pure integer function at(i, n, inc)
        integer, intent(in) :: i, n, inc
        if (inc >= 0) then
            at = 1 + (i - 1)*inc
        else
            at = 1 + (n - i)*(-inc)
        end if
    end function at

    real(dp) function fortnum_ddot(n, x, incx, y, incy) result(s)
        integer, intent(in) :: n, incx, incy
        real(dp), intent(in) :: x(*), y(*)
        integer :: i
        s = 0
        if (n <= 0) return
        if (incx == 1 .and. incy == 1) then
            s = dot_product(x(1:n), y(1:n))
        else
            do i = 1, n
                s = s + x(at(i, n, incx))*y(at(i, n, incy))
            end do
        end if
    end function fortnum_ddot

    real(sp) function fortnum_sdot(n, x, incx, y, incy) result(s)
        integer, intent(in) :: n, incx, incy
        real(sp), intent(in) :: x(*), y(*)
        integer :: i
        s = 0
        do i = 1, n
            s = s + x(at(i, n, incx))*y(at(i, n, incy))
        end do
    end function fortnum_sdot

    !> sum conjg(x) y
    complex(dp) function fortnum_zdotc(n, x, incx, y, incy) result(s)
        integer, intent(in) :: n, incx, incy
        complex(dp), intent(in) :: x(*), y(*)
        integer :: i
        s = 0
        do i = 1, n
            s = s + conjg(x(at(i, n, incx)))*y(at(i, n, incy))
        end do
    end function fortnum_zdotc

    complex(dp) function fortnum_zdotu(n, x, incx, y, incy) result(s)
        integer, intent(in) :: n, incx, incy
        complex(dp), intent(in) :: x(*), y(*)
        integer :: i
        s = 0
        do i = 1, n
            s = s + x(at(i, n, incx))*y(at(i, n, incy))
        end do
    end function fortnum_zdotu

    subroutine fortnum_daxpy(n, alpha, x, incx, y, incy)
        integer, intent(in) :: n, incx, incy
        real(dp), intent(in) :: alpha, x(*)
        real(dp), intent(inout) :: y(*)
        integer :: i
        if (n <= 0 .or. alpha == 0) return
        if (incx == 1 .and. incy == 1) then
            y(1:n) = y(1:n) + alpha*x(1:n)
        else
            do i = 1, n
                y(at(i, n, incy)) = y(at(i, n, incy)) + alpha*x(at(i, n, incx))
            end do
        end if
    end subroutine fortnum_daxpy

    subroutine fortnum_zaxpy(n, alpha, x, incx, y, incy)
        integer, intent(in) :: n, incx, incy
        complex(dp), intent(in) :: alpha, x(*)
        complex(dp), intent(inout) :: y(*)
        integer :: i
        if (n <= 0 .or. alpha == (0.0_dp, 0.0_dp)) return
        do i = 1, n
            y(at(i, n, incy)) = y(at(i, n, incy)) + alpha*x(at(i, n, incx))
        end do
    end subroutine fortnum_zaxpy

    ! Scaling uses the increment as given (reference BLAS ignores incx <= 0).
    subroutine fortnum_dscal(n, alpha, x, incx)
        integer, intent(in) :: n, incx
        real(dp), intent(in) :: alpha
        real(dp), intent(inout) :: x(*)
        integer :: i
        if (n <= 0 .or. incx <= 0) return
        do i = 1, 1 + (n - 1)*incx, incx
            x(i) = alpha*x(i)
        end do
    end subroutine fortnum_dscal

    subroutine fortnum_zscal(n, alpha, x, incx)
        integer, intent(in) :: n, incx
        complex(dp), intent(in) :: alpha
        complex(dp), intent(inout) :: x(*)
        integer :: i
        if (n <= 0 .or. incx <= 0) return
        do i = 1, 1 + (n - 1)*incx, incx
            x(i) = alpha*x(i)
        end do
    end subroutine fortnum_zscal

    subroutine fortnum_zdscal(n, alpha, x, incx)
        integer, intent(in) :: n, incx
        real(dp), intent(in) :: alpha
        complex(dp), intent(inout) :: x(*)
        integer :: i
        if (n <= 0 .or. incx <= 0) return
        do i = 1, 1 + (n - 1)*incx, incx
            x(i) = alpha*x(i)
        end do
    end subroutine fortnum_zdscal

    subroutine fortnum_dcopy(n, x, incx, y, incy)
        integer, intent(in) :: n, incx, incy
        real(dp), intent(in) :: x(*)
        real(dp), intent(inout) :: y(*)
        integer :: i
        do i = 1, n
            y(at(i, n, incy)) = x(at(i, n, incx))
        end do
    end subroutine fortnum_dcopy

    subroutine fortnum_zcopy(n, x, incx, y, incy)
        integer, intent(in) :: n, incx, incy
        complex(dp), intent(in) :: x(*)
        complex(dp), intent(inout) :: y(*)
        integer :: i
        do i = 1, n
            y(at(i, n, incy)) = x(at(i, n, incx))
        end do
    end subroutine fortnum_zcopy

    subroutine fortnum_dswap(n, x, incx, y, incy)
        integer, intent(in) :: n, incx, incy
        real(dp), intent(inout) :: x(*), y(*)
        real(dp) :: t
        integer :: i
        do i = 1, n
            t = x(at(i, n, incx))
            x(at(i, n, incx)) = y(at(i, n, incy))
            y(at(i, n, incy)) = t
        end do
    end subroutine fortnum_dswap

    subroutine fortnum_zswap(n, x, incx, y, incy)
        integer, intent(in) :: n, incx, incy
        complex(dp), intent(inout) :: x(*), y(*)
        complex(dp) :: t
        integer :: i
        do i = 1, n
            t = x(at(i, n, incx))
            x(at(i, n, incx)) = y(at(i, n, incy))
            y(at(i, n, incy)) = t
        end do
    end subroutine fortnum_zswap

    !> Euclidean norm without overflow or harmful underflow: one pass keeping
    !> the largest magnitude seen and the sum of squares scaled by it.
    !> sum |x(i)|; zero for n <= 0 or incx <= 0, as reference BLAS
    real(dp) function fortnum_dasum(n, x, incx) result(s)
        integer, intent(in) :: n, incx
        real(dp), intent(in) :: x(*)
        integer :: i
        s = 0
        if (n <= 0 .or. incx <= 0) return
        do i = 1, 1 + (n - 1)*incx, incx
            s = s + abs(x(i))
        end do
    end function fortnum_dasum

    !> sum |re x(i)| + |im x(i)|; zero for n <= 0 or incx <= 0
    real(dp) function fortnum_dzasum(n, x, incx) result(s)
        integer, intent(in) :: n, incx
        complex(dp), intent(in) :: x(*)
        integer :: i
        s = 0
        if (n <= 0 .or. incx <= 0) return
        do i = 1, 1 + (n - 1)*incx, incx
            s = s + abs(real(x(i), dp)) + abs(aimag(x(i)))
        end do
    end function fortnum_dzasum

    real(dp) function fortnum_dnrm2(n, x, incx) result(norm)
        integer, intent(in) :: n, incx
        real(dp), intent(in) :: x(*)
        real(dp) :: scale, ssq, v
        integer :: i
        norm = 0
        if (n <= 0 .or. incx <= 0) return
        scale = 0
        ssq = 1
        do i = 1, 1 + (n - 1)*incx, incx
            v = abs(x(i))
            if (v /= 0) call accumulate(v, scale, ssq)
            if (v /= v) then
                norm = v
                return
            end if
        end do
        norm = scale*sqrt(ssq)
    end function fortnum_dnrm2

    real(dp) function fortnum_dznrm2(n, x, incx) result(norm)
        integer, intent(in) :: n, incx
        complex(dp), intent(in) :: x(*)
        real(dp) :: scale, ssq, v
        integer :: i, part
        norm = 0
        if (n <= 0 .or. incx <= 0) return
        scale = 0
        ssq = 1
        do i = 1, 1 + (n - 1)*incx, incx
            do part = 1, 2
                v = abs(merge(real(x(i), dp), aimag(x(i)), part == 1))
                if (v /= 0) call accumulate(v, scale, ssq)
                if (v /= v) then
                    norm = v
                    return
                end if
            end do
        end do
        norm = scale*sqrt(ssq)
    end function fortnum_dznrm2

    pure subroutine accumulate(v, scale, ssq)
        real(dp), intent(in) :: v
        real(dp), intent(inout) :: scale, ssq
        if (scale < v) then
            ssq = 1 + ssq*(scale/v)**2
            scale = v
        else
            ssq = ssq + (v/scale)**2
        end if
    end subroutine accumulate

    !> First index of the largest |x(i)| (0 for n < 1).
    integer function fortnum_idamax(n, x, incx) result(k)
        integer, intent(in) :: n, incx
        real(dp), intent(in) :: x(*)
        real(dp) :: big
        integer :: i
        k = 0
        if (n < 1 .or. incx <= 0) return
        k = 1
        big = abs(x(1))
        do i = 2, n
            if (abs(x(1 + (i - 1)*incx)) > big) then
                k = i
                big = abs(x(1 + (i - 1)*incx))
            end if
        end do
    end function fortnum_idamax

    !> First index of the largest |re| + |im| (the reference BLAS measure).
    integer function fortnum_izamax(n, x, incx) result(k)
        integer, intent(in) :: n, incx
        complex(dp), intent(in) :: x(*)
        real(dp) :: big, v
        integer :: i
        k = 0
        if (n < 1 .or. incx <= 0) return
        k = 1
        big = abs(real(x(1), dp)) + abs(aimag(x(1)))
        do i = 2, n
            v = abs(real(x(1 + (i - 1)*incx), dp)) + abs(aimag(x(1 + (i - 1)*incx)))
            if (v > big) then
                k = i
                big = v
            end if
        end do
    end function fortnum_izamax

    !> Plane rotation (x, y) := (c x + s y, c y - s x).
    subroutine fortnum_drot(n, x, incx, y, incy, c, s)
        integer, intent(in) :: n, incx, incy
        real(dp), intent(inout) :: x(*), y(*)
        real(dp), intent(in) :: c, s
        real(dp) :: t
        integer :: i, ix, iy
        do i = 1, n
            ix = at(i, n, incx)
            iy = at(i, n, incy)
            t = c*x(ix) + s*y(iy)
            y(iy) = c*y(iy) - s*x(ix)
            x(ix) = t
        end do
    end subroutine fortnum_drot

    !> y := alpha op(A) x + beta y, op(A) = A or A^T; reference DGEMV.
    subroutine fortnum_dgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        character, intent(in) :: trans
        integer, intent(in) :: m, n, lda, incx, incy
        real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
        real(dp), intent(inout) :: y(*)
        logical :: notrans
        integer :: leny, lenx, i, j
        real(dp) :: t
        select case (trans)
        case ('N', 'n')
            notrans = .true.
        case ('T', 't', 'C', 'c')
            notrans = .false.
        case default
            call argument_error('DGEMV', 1)
        end select
        call check_mv('DGEMV', m, n, lda, incx, incy)
        if (m == 0 .or. n == 0 .or. (alpha == 0 .and. beta == 1)) return
        leny = merge(m, n, notrans)
        lenx = merge(n, m, notrans)
        if (beta /= 1) then
            do i = 1, leny
                if (beta == 0) then
                    y(at(i, leny, incy)) = 0
                else
                    y(at(i, leny, incy)) = beta*y(at(i, leny, incy))
                end if
            end do
        end if
        if (alpha == 0) return
        if (notrans) then
            do j = 1, n
                t = alpha*x(at(j, lenx, incx))
                if (incy == 1) then
                    y(1:m) = y(1:m) + t*a(1:m, j)
                else
                    do i = 1, m
                        y(at(i, m, incy)) = y(at(i, m, incy)) + t*a(i, j)
                    end do
                end if
            end do
        else
            do j = 1, n
                if (incx == 1) then
                    t = dot_product(a(1:m, j), x(1:m))
                else
                    t = 0
                    do i = 1, m
                        t = t + a(i, j)*x(at(i, m, incx))
                    end do
                end if
                y(at(j, n, incy)) = y(at(j, n, incy)) + alpha*t
            end do
        end if
    end subroutine fortnum_dgemv

    !> y := alpha op(A) x + beta y, op(A) = A, A^T or A^H; reference ZGEMV.
    subroutine fortnum_zgemv(trans, m, n, alpha, a, lda, x, incx, beta, y, incy)
        character, intent(in) :: trans
        integer, intent(in) :: m, n, lda, incx, incy
        complex(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
        complex(dp), intent(inout) :: y(*)
        logical :: notrans, conj
        integer :: leny, lenx, i, j
        complex(dp) :: t
        select case (trans)
        case ('N', 'n')
            notrans = .true.
            conj = .false.
        case ('T', 't')
            notrans = .false.
            conj = .false.
        case ('C', 'c')
            notrans = .false.
            conj = .true.
        case default
            call argument_error('ZGEMV', 1)
        end select
        call check_mv('ZGEMV', m, n, lda, incx, incy)
        if (m == 0 .or. n == 0 .or. (alpha == (0.0_dp, 0.0_dp) .and. beta == (1.0_dp, 0.0_dp))) return
        leny = merge(m, n, notrans)
        lenx = merge(n, m, notrans)
        if (beta /= (1.0_dp, 0.0_dp)) then
            do i = 1, leny
                if (beta == (0.0_dp, 0.0_dp)) then
                    y(at(i, leny, incy)) = 0
                else
                    y(at(i, leny, incy)) = beta*y(at(i, leny, incy))
                end if
            end do
        end if
        if (alpha == (0.0_dp, 0.0_dp)) return
        if (notrans) then
            do j = 1, n
                t = alpha*x(at(j, lenx, incx))
                do i = 1, m
                    y(at(i, m, incy)) = y(at(i, m, incy)) + t*a(i, j)
                end do
            end do
        else
            do j = 1, n
                t = 0
                if (conj) then
                    do i = 1, m
                        t = t + conjg(a(i, j))*x(at(i, m, incx))
                    end do
                else
                    do i = 1, m
                        t = t + a(i, j)*x(at(i, m, incx))
                    end do
                end if
                y(at(j, n, incy)) = y(at(j, n, incy)) + alpha*t
            end do
        end if
    end subroutine fortnum_zgemv

    !> y := alpha A x + beta y, A symmetric, only its uplo triangle read.
    subroutine fortnum_dsymv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda, incx, incy
        real(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
        real(dp), intent(inout) :: y(*)
        logical :: upper
        integer :: i, j
        real(dp) :: t1, t2
        upper = upper_flag(uplo, 'DSYMV')
        call check_mv('DSYMV', n, n, lda, incx, incy)
        if (n == 0 .or. (alpha == 0 .and. beta == 1)) return
        do i = 1, n
            if (beta == 0) then
                y(at(i, n, incy)) = 0
            else if (beta /= 1) then
                y(at(i, n, incy)) = beta*y(at(i, n, incy))
            end if
        end do
        if (alpha == 0) return
        do j = 1, n
            t1 = alpha*x(at(j, n, incx))
            t2 = 0
            if (upper) then
                do i = 1, j - 1
                    y(at(i, n, incy)) = y(at(i, n, incy)) + t1*a(i, j)
                    t2 = t2 + a(i, j)*x(at(i, n, incx))
                end do
            else
                do i = j + 1, n
                    y(at(i, n, incy)) = y(at(i, n, incy)) + t1*a(i, j)
                    t2 = t2 + a(i, j)*x(at(i, n, incx))
                end do
            end if
            y(at(j, n, incy)) = y(at(j, n, incy)) + t1*a(j, j) + alpha*t2
        end do
    end subroutine fortnum_dsymv

    !> y := alpha A x + beta y, A Hermitian (diagonal taken as real).
    subroutine fortnum_zhemv(uplo, n, alpha, a, lda, x, incx, beta, y, incy)
        character, intent(in) :: uplo
        integer, intent(in) :: n, lda, incx, incy
        complex(dp), intent(in) :: alpha, beta, a(lda, *), x(*)
        complex(dp), intent(inout) :: y(*)
        logical :: upper
        integer :: i, j
        complex(dp) :: t1, t2
        upper = upper_flag(uplo, 'ZHEMV')
        call check_mv('ZHEMV', n, n, lda, incx, incy)
        if (n == 0 .or. (alpha == (0.0_dp, 0.0_dp) .and. beta == (1.0_dp, 0.0_dp))) return
        do i = 1, n
            if (beta == (0.0_dp, 0.0_dp)) then
                y(at(i, n, incy)) = 0
            else if (beta /= (1.0_dp, 0.0_dp)) then
                y(at(i, n, incy)) = beta*y(at(i, n, incy))
            end if
        end do
        if (alpha == (0.0_dp, 0.0_dp)) return
        do j = 1, n
            t1 = alpha*x(at(j, n, incx))
            t2 = 0
            if (upper) then
                do i = 1, j - 1
                    y(at(i, n, incy)) = y(at(i, n, incy)) + t1*a(i, j)
                    t2 = t2 + conjg(a(i, j))*x(at(i, n, incx))
                end do
            else
                do i = j + 1, n
                    y(at(i, n, incy)) = y(at(i, n, incy)) + t1*a(i, j)
                    t2 = t2 + conjg(a(i, j))*x(at(i, n, incx))
                end do
            end if
            y(at(j, n, incy)) = y(at(j, n, incy)) + t1*real(a(j, j), dp) + alpha*t2
        end do
    end subroutine fortnum_zhemv

    !> x := inv(op(A)) x, A triangular; reference DTRSV.
    subroutine fortnum_dtrsv(uplo, trans, diag, n, a, lda, x, incx)
        character, intent(in) :: uplo, trans, diag
        integer, intent(in) :: n, lda, incx
        real(dp), intent(in) :: a(lda, *)
        real(dp), intent(inout) :: x(*)
        logical :: upper, notrans, unit
        integer :: i, l
        real(dp) :: t
        upper = upper_flag(uplo, 'DTRSV')
        notrans = trans == 'N' .or. trans == 'n'
        if (.not. (notrans .or. index('TtCc', trans) > 0)) call argument_error('DTRSV', 2)
        unit = diag == 'U' .or. diag == 'u'
        if (.not. (unit .or. diag == 'N' .or. diag == 'n')) call argument_error('DTRSV', 3)
        if (n < 0) call argument_error('DTRSV', 4)
        if (lda < max(1, n)) call argument_error('DTRSV', 6)
        if (incx == 0) call argument_error('DTRSV', 8)
        if (n == 0) return
        if (upper .eqv. notrans) then
            ! op(A) upper: back substitution
            do i = n, 1, -1
                t = x(at(i, n, incx))
                do l = i + 1, n
                    t = t - opa(i, l)*x(at(l, n, incx))
                end do
                if (.not. unit) t = t/opa(i, i)
                x(at(i, n, incx)) = t
            end do
        else
            do i = 1, n
                t = x(at(i, n, incx))
                do l = 1, i - 1
                    t = t - opa(i, l)*x(at(l, n, incx))
                end do
                if (.not. unit) t = t/opa(i, i)
                x(at(i, n, incx)) = t
            end do
        end if
    contains
        pure real(dp) function opa(r, c)
            integer, intent(in) :: r, c
            opa = merge(a(r, c), a(c, r), notrans)
        end function opa
    end subroutine fortnum_dtrsv

    !> x := inv(op(A)) x, A triangular, op = A, A^T or A^H; reference ZTRSV.
    subroutine fortnum_ztrsv(uplo, trans, diag, n, a, lda, x, incx)
        character, intent(in) :: uplo, trans, diag
        integer, intent(in) :: n, lda, incx
        complex(dp), intent(in) :: a(lda, *)
        complex(dp), intent(inout) :: x(*)
        logical :: upper, notrans, conj, unit
        integer :: i, l
        complex(dp) :: t
        upper = upper_flag(uplo, 'ZTRSV')
        notrans = trans == 'N' .or. trans == 'n'
        conj = trans == 'C' .or. trans == 'c'
        if (.not. (notrans .or. conj .or. trans == 'T' .or. trans == 't')) call argument_error('ZTRSV', 2)
        unit = diag == 'U' .or. diag == 'u'
        if (.not. (unit .or. diag == 'N' .or. diag == 'n')) call argument_error('ZTRSV', 3)
        if (n < 0) call argument_error('ZTRSV', 4)
        if (lda < max(1, n)) call argument_error('ZTRSV', 6)
        if (incx == 0) call argument_error('ZTRSV', 8)
        if (n == 0) return
        if (upper .eqv. notrans) then
            do i = n, 1, -1
                t = x(at(i, n, incx))
                do l = i + 1, n
                    t = t - opa(i, l)*x(at(l, n, incx))
                end do
                if (.not. unit) t = t/opa(i, i)
                x(at(i, n, incx)) = t
            end do
        else
            do i = 1, n
                t = x(at(i, n, incx))
                do l = 1, i - 1
                    t = t - opa(i, l)*x(at(l, n, incx))
                end do
                if (.not. unit) t = t/opa(i, i)
                x(at(i, n, incx)) = t
            end do
        end if
    contains
        pure complex(dp) function opa(r, c)
            integer, intent(in) :: r, c
            if (notrans) then
                opa = a(r, c)
            else if (conj) then
                opa = conjg(a(c, r))
            else
                opa = a(c, r)
            end if
        end function opa
    end subroutine fortnum_ztrsv

    !> A := A + alpha x y^T; reference DGER.
    subroutine fortnum_dger(m, n, alpha, x, incx, y, incy, a, lda)
        integer, intent(in) :: m, n, incx, incy, lda
        real(dp), intent(in) :: alpha, x(*), y(*)
        real(dp), intent(inout) :: a(lda, *)
        integer :: i, j
        real(dp) :: t
        if (m < 0) call argument_error('DGER', 1)
        if (n < 0) call argument_error('DGER', 2)
        if (incx == 0) call argument_error('DGER', 5)
        if (incy == 0) call argument_error('DGER', 7)
        if (lda < max(1, m)) call argument_error('DGER', 9)
        if (m == 0 .or. n == 0 .or. alpha == 0) return
        do j = 1, n
            t = alpha*y(at(j, n, incy))
            if (incx == 1) then
                a(1:m, j) = a(1:m, j) + t*x(1:m)
            else
                do i = 1, m
                    a(i, j) = a(i, j) + t*x(at(i, m, incx))
                end do
            end if
        end do
    end subroutine fortnum_dger

    !> A := A + alpha x y^T (complex, unconjugated); reference ZGERU.
    subroutine fortnum_zgeru(m, n, alpha, x, incx, y, incy, a, lda)
        integer, intent(in) :: m, n, incx, incy, lda
        complex(dp), intent(in) :: alpha, x(*), y(*)
        complex(dp), intent(inout) :: a(lda, *)
        call zger('ZGERU', m, n, alpha, x, incx, y, incy, a, lda, .false.)
    end subroutine fortnum_zgeru

    !> A := A + alpha x y^H; reference ZGERC.
    subroutine fortnum_zgerc(m, n, alpha, x, incx, y, incy, a, lda)
        integer, intent(in) :: m, n, incx, incy, lda
        complex(dp), intent(in) :: alpha, x(*), y(*)
        complex(dp), intent(inout) :: a(lda, *)
        call zger('ZGERC', m, n, alpha, x, incx, y, incy, a, lda, .true.)
    end subroutine fortnum_zgerc

    subroutine zger(name, m, n, alpha, x, incx, y, incy, a, lda, conj)
        character(*), intent(in) :: name
        integer, intent(in) :: m, n, incx, incy, lda
        complex(dp), intent(in) :: alpha, x(*), y(*)
        complex(dp), intent(inout) :: a(lda, *)
        logical, intent(in) :: conj
        integer :: i, j
        complex(dp) :: t
        if (m < 0) call argument_error(name, 1)
        if (n < 0) call argument_error(name, 2)
        if (incx == 0) call argument_error(name, 5)
        if (incy == 0) call argument_error(name, 7)
        if (lda < max(1, m)) call argument_error(name, 9)
        if (m == 0 .or. n == 0 .or. alpha == (0.0_dp, 0.0_dp)) return
        do j = 1, n
            t = y(at(j, n, incy))
            if (conj) t = conjg(t)
            t = alpha*t
            do i = 1, m
                a(i, j) = a(i, j) + t*x(at(i, m, incx))
            end do
        end do
    end subroutine zger

    subroutine check_mv(name, m, n, lda, incx, incy)
        character(*), intent(in) :: name
        integer, intent(in) :: m, n, lda, incx, incy
        if (m < 0) call argument_error(name, 2)
        if (n < 0) call argument_error(name, 3)
        if (lda < max(1, m)) call argument_error(name, 6)
        if (incx == 0) call argument_error(name, 8)
        if (incy == 0) call argument_error(name, 11)
    end subroutine check_mv

    logical function upper_flag(uplo, name)
        character, intent(in) :: uplo
        character(*), intent(in) :: name
        select case (uplo)
        case ('U', 'u')
            upper_flag = .true.
        case ('L', 'l')
            upper_flag = .false.
        case default
            upper_flag = .false.
            call argument_error(name, 1)
        end select
    end function upper_flag

    subroutine argument_error(routine, position)
        character(*), intent(in) :: routine
        integer, intent(in) :: position
        character(len=16) :: number
        write (number, '(i0)') position
        error stop routine//': argument '//trim(number)//' has an illegal value'
    end subroutine argument_error

end module fortnum_blas12
