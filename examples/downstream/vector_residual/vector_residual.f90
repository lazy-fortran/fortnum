program vector_residual
    ! Downstream example: vector residual through the #41 optimizer-facing API.
    !
    ! Residual r: R^3 -> R^3
    !   r1(x) = x1^2 + sin(x2)
    !   r2(x) = x2 * x3
    !   r3(x) = x1 * x3 - x2^2
    !
    ! Analytic Jacobian J:
    !   J = [ 2x1   cos(x2)   0  ]
    !       [ 0     x3        x2 ]
    !       [ x3   -2x2       x1 ]
    !
    ! Checks performed:
    !   (1) value matches analytic formula
    !   (2) JVP matches central finite difference
    !   (3) dot-product adjoint identity: u.(Jv) = v.(J^T u)
    !   (4) VJP row check: J^T e_i matches row i of J from analytic J
    !
    ! The active vector has a single block "params" (size 3). The context
    ! is the layout. This shows how the abstract interfaces are wired: the
    ! kernel unpacks x from the layout and calls the residual.
    !
    ! Build standalone (from repo root):
    !   gfortran -I<build>/include -o vector_residual \
    !       examples/downstream/vector_residual/vector_residual.f90 \
    !       -L<build>/lib -lfortnum
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use fortnum_kinds,         only: dp
    use fortnum_status,        only: fortnum_status_t, status_ok, FORTNUM_OK
    use fortnum_active_vector, only: fortnum_active_layout_t, layout_init, &
                                     layout_add, pack_block, unpack_block
    use fortnum_ad_interfaces, only: fortnum_ad_status_t, ad_status_set, &
                                     ad_status_ok, &
                                     FORTNUM_AD_BACKEND_ANALYTIC, &
                                     FORTNUM_AD_QUALITY_EXACT
    implicit none

    integer :: nfail
    nfail = 0

    call run_checks(nfail)

    if (nfail > 0) then
        write (error_unit, '(i0,a)') nfail, " test(s) failed"
        stop 1
    end if
    write (*, '(a)') "PASS"
    stop 0

contains

    subroutine run_checks(nfail)
        integer, intent(inout) :: nfail

        type(fortnum_active_layout_t) :: layout
        type(fortnum_status_t)        :: us
        type(fortnum_ad_status_t)     :: st
        real(dp) :: x(3), y(3), y_ref(3)
        real(dp) :: v(3), jv_ad(3), jv_fd(3)
        real(dp) :: u(3), jtu_ad(3)
        real(dp) :: err, h, lhs, rhs
        integer  :: i

        ! Build layout: one block "params" of length 3.
        call layout_init(layout, 1)
        call layout_add(layout, "params", 3, us)
        if (.not. status_ok(us)) error stop "layout_add failed"

        ! Base point.
        x = [0.5_dp, 1.3_dp, -0.8_dp]
        call pack_block(layout, x, "params", x, us)

        ! ------------------------------------------------------------ (1) value
        call res_value(3, x, y, layout, st)
        if (.not. ad_status_ok(st)) then
            write (error_unit, '(a)') "FAIL [value] ad_status not ok"
            nfail = nfail + 1
        end if
        y_ref(1) = x(1)**2 + sin(x(2))
        y_ref(2) = x(2)*x(3)
        y_ref(3) = x(1)*x(3) - x(2)**2
        err = maxval(abs(y - y_ref)) / max(maxval(abs(y_ref)), 1.0_dp)
        if (err > 1.0e-14_dp) then
            write (error_unit, '(a,es12.4)') "FAIL [value] rel_err=", err
            nfail = nfail + 1
        end if

        ! ------------------------------------------------------------ (2) JVP vs FD
        v  = [0.2_dp, -0.4_dp, 1.1_dp]
        call res_jvp(3, x, v, y, jv_ad, layout, st)

        h = 1.0e-6_dp
        call fd_jvp_vec(layout, x, v, h, jv_fd)
        err = maxval(abs(jv_ad - jv_fd)) / max(maxval(abs(jv_fd)), 1.0_dp)
        if (err > 1.0e-6_dp) then
            write (error_unit, '(a,es12.4)') "FAIL [JVP-vs-FD] rel_err=", err
            nfail = nfail + 1
        end if

        ! ------------------------------------------------------------ (3) adjoint identity
        u = [1.0_dp, -0.5_dp, 0.3_dp]
        call res_vjp(3, x, u, jtu_ad, layout, st)
        lhs = dot_product(u, jv_ad)        ! u . (J v)
        rhs = dot_product(v, jtu_ad)       ! v . (J^T u)
        err = abs(lhs - rhs) / max(abs(lhs), abs(rhs), 1.0_dp)
        if (err > 1.0e-13_dp) then
            write (error_unit, '(a,es12.4,a,es24.16,a,es24.16)') &
                "FAIL [adjoint] rel_err=", err, "  lhs=", lhs, "  rhs=", rhs
            nfail = nfail + 1
        end if

        ! ------------------------------------------------------------ (4) VJP row check
        ! VJP with u = e_i returns J^T e_i = row i of J (a vector of length n).
        ! J rows:
        !   row 1: [2x1, cos(x2), 0]
        !   row 2: [0, x3, x2]
        !   row 3: [x3, -2x2, x1]
        block
            real(dp) :: ei(3), jtu(3), jrow_ref(3)
            real(dp) :: worst
            worst = 0.0_dp
            do i = 1, 3
                ei = 0.0_dp; ei(i) = 1.0_dp
                call res_vjp(3, x, ei, jtu, layout, st)
                select case (i)
                case (1)
                    jrow_ref = [2.0_dp*x(1), cos(x(2)), 0.0_dp]
                case (2)
                    jrow_ref = [0.0_dp, x(3), x(2)]
                case (3)
                    jrow_ref = [x(3), -2.0_dp*x(2), x(1)]
                end select
                worst = max(worst, maxval(abs(jtu - jrow_ref)) &
                    / max(maxval(abs(jrow_ref)), 1.0_dp))
            end do
            if (worst > 1.0e-14_dp) then
                write (error_unit, '(a,es12.4)') &
                    "FAIL [VJP-rows] rel_err=", worst
                nfail = nfail + 1
            end if
        end block

    end subroutine run_checks

    ! value_fn: y = r(x).
    subroutine res_value(n, x, y, layout, st)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x(n)
        real(dp), intent(out) :: y(:)
        type(fortnum_active_layout_t), intent(in)  :: layout
        type(fortnum_ad_status_t),     intent(out) :: st
        real(dp) :: p(3)
        type(fortnum_status_t) :: us
        call unpack_block(layout, x, "params", p, us)
        y(1) = p(1)**2 + sin(p(2))
        y(2) = p(2)*p(3)
        y(3) = p(1)*p(3) - p(2)**2
        call ad_status_set(st, FORTNUM_OK, "", &
            FORTNUM_AD_BACKEND_ANALYTIC, FORTNUM_AD_QUALITY_EXACT)
    end subroutine res_value

    ! jvp_fn: (y, y_dot) = (r(x), J(x) v).
    subroutine res_jvp(n, x, v, y, y_dot, layout, st)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x(n), v(n)
        real(dp), intent(out) :: y(:), y_dot(:)
        type(fortnum_active_layout_t), intent(in)  :: layout
        type(fortnum_ad_status_t),     intent(out) :: st
        real(dp) :: p(3)
        type(fortnum_status_t) :: us
        call unpack_block(layout, x, "params", p, us)
        y(1) = p(1)**2 + sin(p(2))
        y(2) = p(2)*p(3)
        y(3) = p(1)*p(3) - p(2)**2
        ! J v: J row i times v, using analytic Jacobian.
        y_dot(1) = 2.0_dp*p(1)*v(1) + cos(p(2))*v(2)
        y_dot(2) = p(3)*v(2) + p(2)*v(3)
        y_dot(3) = p(3)*v(1) - 2.0_dp*p(2)*v(2) + p(1)*v(3)
        call ad_status_set(st, FORTNUM_OK, "", &
            FORTNUM_AD_BACKEND_ANALYTIC, FORTNUM_AD_QUALITY_EXACT)
    end subroutine res_jvp

    ! vjp_fn: x_bar = J(x)^T y_bar.
    subroutine res_vjp(n, x, y_bar, x_bar, layout, st)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: x(n), y_bar(:)
        real(dp), intent(out) :: x_bar(n)
        type(fortnum_active_layout_t), intent(in)  :: layout
        type(fortnum_ad_status_t),     intent(out) :: st
        real(dp) :: p(3), u(3)
        type(fortnum_status_t) :: us
        call unpack_block(layout, x, "params", p, us)
        u = y_bar(1:3)
        ! J^T u: column j of J dotted with u.
        x_bar(1) = 2.0_dp*p(1)*u(1) + p(3)*u(3)
        x_bar(2) = cos(p(2))*u(1) + p(3)*u(2) - 2.0_dp*p(2)*u(3)
        x_bar(3) = p(2)*u(2) + p(1)*u(3)
        call ad_status_set(st, FORTNUM_OK, "", &
            FORTNUM_AD_BACKEND_ANALYTIC, FORTNUM_AD_QUALITY_EXACT)
    end subroutine res_vjp

    ! Central FD of J(x) v.
    subroutine fd_jvp_vec(layout, x, v, h, jv)
        type(fortnum_active_layout_t), intent(in)  :: layout
        real(dp), intent(in)  :: x(:), v(:), h
        real(dp), intent(out) :: jv(:)
        type(fortnum_ad_status_t) :: st
        real(dp) :: xp(size(x)), xm(size(x)), yp(size(jv)), ym(size(jv))
        xp = x + h*v
        xm = x - h*v
        call res_value(size(x), xp, yp, layout, st)
        call res_value(size(x), xm, ym, layout, st)
        jv = (yp - ym) / (2.0_dp*h)
    end subroutine fd_jvp_vec

end program vector_residual
