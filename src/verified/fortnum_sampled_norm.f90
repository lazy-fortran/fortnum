module fortnum_sampled_norm
    use, intrinsic :: iso_fortran_env, only: dp => real64, i64 => int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_rng, only: rng_t, rng_next_u64
    use fortnum_interval, only: interval_t, ipoint, mag, operator(+), &
        operator(*), operator(/), log, sqrt
    implicit none
    private
    public :: sampled_norm_upper
    abstract interface
        subroutine coordinate_enclosure(index, value, ok)
            import :: interval_t
            integer, intent(in) :: index
            type(interval_t), intent(out) :: value
            logical, intent(out) :: ok
        end subroutine
    end interface
contains
    ! Probability law: iid UNIFORM coordinates, independent of a fixed vector.
    ! A caller must prove |v_i| <= envelope for ALL coordinates, not a sampled max.
    ! Threefry is the reproducible implementation; confidence is conditional on
    ! the ideal iid-bit law (a deterministic PRNG seed alone is not random).
    subroutine sampled_norm_upper(n, samples, delta, envelope, draw, coordinate, &
                                   estimate, upper, ok)
        integer, intent(in) :: n, samples
        real(dp), intent(in) :: delta, envelope
        type(rng_t), intent(inout) :: draw
        procedure(coordinate_enclosure) :: coordinate
        real(dp), intent(out) :: estimate, upper
        logical, intent(out) :: ok
        type(interval_t) :: value, scaled, sumsq, tail, radius
        integer(i64) :: bits, limit
        integer :: j,index
        logical :: valid
        ok = .false.; estimate = 0; upper = huge(0._dp)
        if (n < 1 .or. samples < 1) return
        if (.not. ieee_is_finite(delta) .or. .not. ieee_is_finite(envelope)) return
        if (delta <= 0 .or. delta >= 1 .or. envelope <= 0) return
        limit = huge(bits) - modulo(huge(bits),int(n,i64))
        sumsq = ipoint(0._dp)
        do j = 1,samples
            ! Rejection avoids modulo bias for every n (including non-powers of 2).
            do
                call rng_next_u64(draw,bits)
                bits = shiftr(bits,1)
                if (bits < limit) exit
            end do
            index = int(modulo(bits,int(n,i64))) + 1
            call coordinate(index,value,valid)
            if (.not. valid) return
            if (.not. ieee_is_finite(value%lo)) return
            if (.not. ieee_is_finite(value%hi)) return
            if (value%lo > value%hi) return
            if (mag(value) > envelope) return
            scaled = ipoint(mag(value))/ipoint(envelope)
            sumsq = sumsq + scaled*scaled
        end do
        ! Hoeffding: P(E X > mean(X)+sqrt(log(1/delta)/(2m))) <= delta,
        ! X=(v_I/envelope)^2 in [0,1]. Every arithmetic bridge is outward.
        tail = sqrt(log(ipoint(1._dp)/ipoint(delta))/ &
            (ipoint(2._dp)*ipoint(real(samples,dp))))
        scaled = sumsq/ipoint(real(samples,dp))
        estimate = envelope*sqrt(real(n,dp)*scaled%hi)
        scaled = scaled + tail
        radius = ipoint(envelope)*sqrt(ipoint(real(n,dp))* &
            ipoint(min(1._dp,scaled%hi)))
        upper = radius%hi
        ok = ieee_is_finite(upper)
    end subroutine
end module
