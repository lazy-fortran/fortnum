! Fixed bicubic coefficient polynomial and analytical derivative products.
! c(i+1,j+1) multiplies x**i*y**j, i,j=0..3. jet=[P,P_x,P_y].
! No allocations or mutable state. Inputs must not alias intent(out) results.
! Floating-point centers only: see docs/bicubic.md for rounding/domain limits.
module fortnum_bicubic
    use fortnum_kinds, only: dp
    use fortnum_bicubic_kernel, only: bicubic_kernel, bicubic_jvp_kernel, &
                                      bicubic_vjp_kernel
    implicit none
    private
    public :: bicubic_jet, bicubic_jet_jvp, bicubic_jet_vjp
contains
    pure subroutine bicubic_jet(c, x, y, jet)
        real(dp), intent(in) :: c(4, 4), x, y
        real(dp), intent(out) :: jet(3)
        call bicubic_kernel( &
            x, y, c(1, 1), c(2, 1), &
            c(3, 1), c(4, 1), c(1, 2), c(2, 2), &
            c(3, 2), c(4, 2), c(1, 3), c(2, 3), &
            c(3, 3), c(4, 3), c(1, 4), c(2, 4), &
            c(3, 4), c(4, 4), jet(1), jet(2), &
            jet(3))
    end subroutine bicubic_jet

    ! All coefficients and both coordinates are active; returns primal and JVP.
    pure subroutine bicubic_jet_jvp(c, x, y, dc, dx, dy, jet, djet)
        real(dp), intent(in) :: c(4, 4), x, y, dc(4, 4), dx, dy
        real(dp), intent(out) :: jet(3), djet(3)
        call bicubic_jvp_kernel( &
            x, y, c(1, 1), c(2, 1), &
            c(3, 1), c(4, 1), c(1, 2), c(2, 2), &
            c(3, 2), c(4, 2), c(1, 3), c(2, 3), &
            c(3, 3), c(4, 3), c(1, 4), c(2, 4), &
            c(3, 4), c(4, 4), dx, dy, &
            dc(1, 1), dc(2, 1), dc(3, 1), dc(4, 1), &
            dc(1, 2), dc(2, 2), dc(3, 2), dc(4, 2), &
            dc(1, 3), dc(2, 3), dc(3, 3), dc(4, 3), &
            dc(1, 4), dc(2, 4), dc(3, 4), dc(4, 4), &
            jet(1), jet(2), jet(3), djet(1), &
            djet(2), djet(3))
    end subroutine bicubic_jet_jvp

    ! Euclidean pairing, overwrites (rather than accumulates) input cotangents.
    pure subroutine bicubic_jet_vjp(c, x, y, bar_jet, jet, bar_c, bar_x, bar_y)
        real(dp), intent(in) :: c(4, 4), x, y, bar_jet(3)
        real(dp), intent(out) :: jet(3), bar_c(4, 4), bar_x, bar_y
        call bicubic_vjp_kernel( &
            x, y, c(1, 1), c(2, 1), &
            c(3, 1), c(4, 1), c(1, 2), c(2, 2), &
            c(3, 2), c(4, 2), c(1, 3), c(2, 3), &
            c(3, 3), c(4, 3), c(1, 4), c(2, 4), &
            c(3, 4), c(4, 4), bar_jet(1), bar_jet(2), &
            bar_jet(3), jet(1), jet(2), jet(3), &
            bar_x, bar_y, bar_c(1, 1), bar_c(2, 1), &
            bar_c(3, 1), bar_c(4, 1), bar_c(1, 2), bar_c(2, 2), &
            bar_c(3, 2), bar_c(4, 2), bar_c(1, 3), bar_c(2, 3), &
            bar_c(3, 3), bar_c(4, 3), bar_c(1, 4), bar_c(2, 4), &
            bar_c(3, 4), bar_c(4, 4))
    end subroutine bicubic_jet_vjp
end module fortnum_bicubic
