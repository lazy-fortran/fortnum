program test_gram_factor
    use fortnum_kinds, only: dp
    use fortnum_status, only: fortnum_status_t, status_ok
    use fortnum_gram_factor, only: gram_factor_workspace_t
    use fortnum_psd_factor, only: psd_factor
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    implicit none
    integer, parameter :: qp = selected_real_kind(33, 4931)
    type(gram_factor_workspace_t) :: ws
    type(fortnum_status_t) :: status
    real(dp) :: w(3, 5), b(3, 3), a(3, 3), small(2, 2), root(2, 2), err
    real(dp) :: empty(3, 0), bad(2, 3), zero(0, 0), zero_out(0, 0), scales(3)
    integer :: i, j, k, rank

    do j = 1, 5
        do i = 1, 3
            w(i, j) = sin(real(i+3*j, dp))
        end do
    end do
    call ws%factor(w, b, status)
    call require(status_ok(status), 'wide Gram factor')
    call gram_oracle(w, b)
    call ws%factor(w(:, :2), b, status)
    call require(status_ok(status), 'tall Gram factor')
    call gram_oracle(w(:, :2), b)
    a = matmul(w(:, :2), transpose(w(:, :2)))
    call psd_factor(a, b, status, rank)
    call require(status_ok(status), 'dense rank deficient PSD')
    call require(rank >= 2, 'dense rank deficient resolved rank')
    call psd_oracle(a, b)
    call ws%factor(empty, b, status)
    call require(status_ok(status), 'empty column Gram factor')
    call require(all(b == 0), 'empty column zero root')
    call ws%factor(w, bad, status)
    call require(.not. status_ok(status), 'wrong output shape')
    call ws%factor(zero, zero_out, status)
    call require(.not. status_ok(status), 'empty row Gram rejected')
    w(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call ws%factor(w, b, status)
    call require(.not. status_ok(status), 'nonfinite Gram input')
    w = 1.0e-320_dp
    call ws%factor(w, b, status)
    call require(status_ok(status), 'subnormal root QR')
    call gram_oracle(w, b, 2.0e-3_qp)
    scales = [1.0e-150_dp, 1.0_dp, 1.0e150_dp]
    do k = 1, 3
        w = reshape([1.0_dp, 2.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, &
            1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, &
            0.0_dp, -1.0_dp, 0.0_dp, 2.0_dp], [3, 5])*scales(k)
        call ws%factor(w, b, status)
        call require(status_ok(status), 'scaled Gram factor')
        call gram_oracle(w, b)
        a = matmul(w, transpose(w))
        call psd_factor(a, b, status, rank, err)
        call require(status_ok(status), 'scaled dense PSD factor')
        call require(rank == 3, 'dense positive rank')
        call psd_oracle(a, b)
    end do
    small = reshape([0.25_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2, 2])
    call psd_factor(small, root, status, rank, err)
    call require(status_ok(status), 'binary rank one PSD')
    call require(rank == 1, 'rank one numerical rank')
    call psd_oracle(small, root)
    small = 0
    call psd_factor(small, root, status, rank, err)
    call require(status_ok(status), 'zero PSD')
    call require(rank == 0 .and. all(root == 0), 'zero PSD rank and root')
    small(1, 1) = 1.0e12_dp
    small(2, 2) = 1.0e-20_dp
    call psd_factor(small, root, status, rank)
    call require(status_ok(status), 'mixed scale positive PSD')
    call require(rank == 2, 'tiny positive pivot retained')
    call require(abs(norm2(root(2, :))/1.0e-10_dp-1) < 1.0e-14_dp, &
        'tiny positive diagonal reconstructed')
    call psd_oracle(small, root)
    small = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2, 2])
    call psd_factor(small, root, status)
    call require(.not. status_ok(status), 'indefinite PSD rejected')
    small = reshape([1.0_dp, 0.0_dp, 0.1_dp, 1.0_dp], [2, 2])
    call psd_factor(small, root, status)
    call require(.not. status_ok(status), 'asymmetric PSD rejected')
    small(1, 1) = ieee_value(0.0_dp, ieee_quiet_nan)
    call psd_factor(small, root, status)
    call require(.not. status_ok(status), 'nonfinite PSD rejected')
    call psd_factor(zero, zero_out, status)
    call require(.not. status_ok(status), 'empty PSD rejected')
    call psd_factor(bad, root, status)
    call require(.not. status_ok(status), 'nonsquare PSD rejected')
contains
    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if (.not. ok) then
            print *, message
            stop 1
        end if
    end subroutine
    subroutine gram_oracle(matrix, factor, tolerance)
        real(dp), intent(in) :: matrix(:, :), factor(:, :)
        real(qp), intent(in), optional :: tolerance
        real(qp) :: scale, expected, actual, tol
        integer :: i, j
        scale = maxval(abs(real(matrix, qp)))
        tol = 2.0e-13_qp
        if (present(tolerance)) tol = tolerance
        do j = 1, size(matrix, 1)
            do i = 1, size(matrix, 1)
                expected = sum((real(matrix(i, :), qp)/scale)* &
                    (real(matrix(j, :), qp)/scale))
                actual = sum((real(factor(i, :), qp)/scale)* &
                    (real(factor(j, :), qp)/scale))
                call require(abs(actual-expected) <= tol*max(1.0_qp, abs(expected)), &
                    'quad precision Gram oracle')
            end do
        end do
    end subroutine
    subroutine psd_oracle(matrix, factor)
        real(dp), intent(in) :: matrix(:, :), factor(:, :)
        real(qp) :: scale, actual, expected
        integer :: i, j
        scale = maxval(abs(real(matrix, qp)))
        do j = 1, size(matrix, 1)
            do i = 1, size(matrix, 1)
                expected = real(matrix(i, j), qp)/scale
                actual = sum((real(factor(i, :), qp)/sqrt(scale))* &
                    (real(factor(j, :), qp)/sqrt(scale)))
                call require(abs(actual-expected) < 2.0e-13_qp, &
                    'quad precision PSD oracle')
            end do
        end do
    end subroutine
end program
