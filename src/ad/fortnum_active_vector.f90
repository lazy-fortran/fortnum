module fortnum_active_vector
    ! Flat active-vector layout for optimizer-facing kernels (issue #41).
    !
    ! A downstream optimization code packs heterogeneous active inputs --
    ! boundary modes, radial profiles, kernel parameters -- into a single
    ! contiguous x(:) and hands that vector to a value/JVP/VJP/grad/HVP kernel.
    ! The layout records the named blocks (offset and size) so the kernel can
    ! unpack the slice it needs inside its own context without the optimizer
    ! knowing the internal ordering.
    !
    ! Pure, no module-level state. The layout is data the caller owns and
    ! threads through every call; nothing here mutates shared memory.
    use fortnum_kinds,  only: dp
    use fortnum_status, only: fortnum_status_t, status_set, &
        FORTNUM_OK, FORTNUM_DOMAIN_ERROR
    implicit none
    private

    ! Maximum length of a component name.
    integer, parameter, public :: FORTNUM_AV_NAME_LEN = 32

    ! One named block within the flat vector: a contiguous slice
    ! x(offset : offset + size - 1) carrying one logical quantity.
    type, public :: fortnum_av_block_t
        character(FORTNUM_AV_NAME_LEN) :: name = ""
        integer :: offset = 0 ! 1-based start index into x(:)
        integer :: size   = 0 ! number of entries in this block
    end type fortnum_av_block_t

    ! Layout descriptor: total length n and the named blocks that tile it.
    ! Built incrementally with layout_add, then consumed by pack/unpack. The
    ! blocks are stored in declaration order; offsets are assigned so the
    ! blocks tile [1, n] without gaps or overlap.
    type, public :: fortnum_active_layout_t
        integer :: n      = 0 ! total flat length
        integer :: nblock = 0 ! number of named blocks
        type(fortnum_av_block_t), allocatable :: blocks(:)
    end type fortnum_active_layout_t

    public :: layout_init
    public :: layout_add
    public :: layout_index
    public :: layout_block
    public :: pack_block
    public :: unpack_block

contains

    ! Initializes an empty layout with capacity for cap named blocks. The
    ! caller adds blocks with layout_add; n grows as blocks are appended.
    pure subroutine layout_init(layout, cap)
        type(fortnum_active_layout_t), intent(out) :: layout
        integer,                       intent(in)  :: cap
        layout%n      = 0
        layout%nblock = 0
        allocate (layout%blocks(max(cap, 0)))
    end subroutine layout_init

    ! Appends a named block of length blk_size to the layout. The block is
    ! placed immediately after the previous one, so n advances by blk_size and
    ! the new block's offset is the old n + 1. Reports FORTNUM_DOMAIN_ERROR on
    ! a non-positive size, a name collision, or exhausted capacity.
    pure subroutine layout_add(layout, name, blk_size, status)
        type(fortnum_active_layout_t), intent(inout) :: layout
        character(*),                  intent(in)    :: name
        integer,                       intent(in)    :: blk_size
        type(fortnum_status_t),        intent(out)   :: status
        integer :: i

        call status_set(status, FORTNUM_OK, "")
        if (blk_size <= 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "layout_add: block size must be positive")
            return
        end if
        if (.not. allocated(layout%blocks) .or. &
            layout%nblock >= size(layout%blocks)) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "layout_add: layout capacity exhausted")
            return
        end if
        do i = 1, layout%nblock
            if (trim(layout%blocks(i)%name) == trim(name)) then
                call status_set(status, FORTNUM_DOMAIN_ERROR, &
                    "layout_add: duplicate block name")
                return
            end if
        end do

        layout%nblock = layout%nblock + 1
        layout%blocks(layout%nblock)%name   = name
        layout%blocks(layout%nblock)%offset = layout%n + 1
        layout%blocks(layout%nblock)%size   = blk_size
        layout%n = layout%n + blk_size
    end subroutine layout_add

    ! Returns the 1-based position of the named block, or 0 if absent.
    pure function layout_index(layout, name) result(idx)
        type(fortnum_active_layout_t), intent(in) :: layout
        character(*),                  intent(in) :: name
        integer                                   :: idx
        integer :: i
        idx = 0
        do i = 1, layout%nblock
            if (trim(layout%blocks(i)%name) == trim(name)) then
                idx = i
                return
            end if
        end do
    end function layout_index

    ! Returns the block descriptor for a name. Reports FORTNUM_DOMAIN_ERROR and
    ! leaves blk zeroed when the name is unknown.
    pure subroutine layout_block(layout, name, blk, status)
        type(fortnum_active_layout_t), intent(in)  :: layout
        character(*),                  intent(in)  :: name
        type(fortnum_av_block_t),      intent(out) :: blk
        type(fortnum_status_t),        intent(out) :: status
        integer :: idx
        call status_set(status, FORTNUM_OK, "")
        idx = layout_index(layout, name)
        if (idx == 0) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "layout_block: unknown block name")
            return
        end if
        blk = layout%blocks(idx)
    end subroutine layout_block

    ! Writes vals into the named block's slice of the flat vector x. Validates
    ! that x is long enough and that size(vals) matches the block size.
    pure subroutine pack_block(layout, x, name, vals, status)
        type(fortnum_active_layout_t), intent(in)    :: layout
        real(dp),                      intent(inout) :: x(:)
        character(*),                  intent(in)    :: name
        real(dp),                      intent(in)    :: vals(:)
        type(fortnum_status_t),        intent(out)   :: status
        type(fortnum_av_block_t) :: blk
        integer :: lo, hi

        call layout_block(layout, name, blk, status)
        if (status%code /= FORTNUM_OK) return
        if (size(vals) /= blk%size .or. size(x) < layout%n) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "pack_block: size mismatch")
            return
        end if
        lo = blk%offset
        hi = blk%offset + blk%size - 1
        x(lo:hi) = vals
    end subroutine pack_block

    ! Reads the named block's slice of the flat vector x into vals. Validates
    ! that x is long enough and that size(vals) matches the block size.
    pure subroutine unpack_block(layout, x, name, vals, status)
        type(fortnum_active_layout_t), intent(in)  :: layout
        real(dp),                      intent(in)  :: x(:)
        character(*),                  intent(in)  :: name
        real(dp),                      intent(out) :: vals(:)
        type(fortnum_status_t),        intent(out) :: status
        type(fortnum_av_block_t) :: blk
        integer :: lo, hi

        call layout_block(layout, name, blk, status)
        if (status%code /= FORTNUM_OK) return
        if (size(vals) /= blk%size .or. size(x) < layout%n) then
            call status_set(status, FORTNUM_DOMAIN_ERROR, &
                "unpack_block: size mismatch")
            return
        end if
        lo = blk%offset
        hi = blk%offset + blk%size - 1
        vals = x(lo:hi)
    end subroutine unpack_block

end module fortnum_active_vector
