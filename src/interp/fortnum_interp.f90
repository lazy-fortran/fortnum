module fortnum_interp
    ! Binary search / cell-index location into a sorted grid.
    !
    ! DERIVATIVE POLICY (ad.md §1, §4):
    !   grid_search: primal_only for the index output.
    !   The index i returned by grid_search is the result of control flow
    !   branching on the grid values; it is an integer and carries no derivative.
    !   Smooth derivatives are valid only INSIDE a fixed cell [xp(i-1), xp(i)].
    !   Crossing a cell boundary is a non-smooth event: the index jumps
    !   discontinuously, so any function that calls grid_search and then
    !   evaluates a piecewise polynomial must treat the cell boundary as a
    !   derivative discontinuity. The AD contract (ad.md §3) classifies integer
    !   outputs as inactive; callers must hold the index fixed when computing
    !   derivatives of the interpolant with respect to the evaluation point x.
    !
    ! grid_search(p, nmin, nmax, xi, i):
    !   Finds the index i in [nmin+1, nmax] such that
    !     p(i-1) < xi <= p(i)   (interior)
    !   Clamping at boundaries:
    !     xi <= p(nmin)  =>  i = nmin + 1  (clamp low)
    !     xi >= p(nmax)  =>  i = nmax      (clamp high)
    !   p(nmin:nmax) must be strictly increasing; no check is performed.
    !   The algorithm is a binary bisection identical in behavior to the libneo
    !   binsrc routine (clean-room reimplementation).
    !
    ! Complexity: O(log2(nmax - nmin)) comparisons.

    use fortnum_kinds, only: dp
    implicit none
    private

    public :: grid_search

contains

    ! Binary search in a strictly increasing array p(nmin:nmax).
    ! Returns the cell index i such that p(i-1) < xi <= p(i), clamped.
    pure subroutine grid_search(p, nmin, nmax, xi, i)
        integer,  intent(in)  :: nmin, nmax
        real(dp), intent(in)  :: p(nmin:nmax)
        real(dp), intent(in)  :: xi
        integer,  intent(out) :: i

        integer :: imin, imax, imid, k, n

        ! Clamp-low: return the first interior cell regardless of how far below
        ! the grid xi falls; mirrors binsrc behavior.
        if (xi <= p(nmin)) then
            i = nmin + 1
            return
        end if

        ! Clamp-high.
        if (xi >= p(nmax)) then
            i = nmax
            return
        end if

        imin = nmin
        imax = nmax
        n    = nmax - nmin

        ! Bisection: at each step narrow [imin, imax] by one half.
        ! Invariant: p(imin) < xi <= p(imax).
        do k = 1, n
            imid = (imax - imin)/2 + imin
            if (p(imid) >= xi) then
                imax = imid
            else
                imin = imid
            end if
            if (imax == imin + 1) exit
        end do

        i = imax
    end subroutine grid_search

end module fortnum_interp
