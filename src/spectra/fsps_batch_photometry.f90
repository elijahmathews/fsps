module fsps_batch_photometry
    !> @brief Magnitudes of many spectra at once on the device (batched CSP, compute_csp_scenario).
    !>
    !> @details
    !> Same quantities as compute_magnitudes, for the columns of a spectrum array resident on the
    !> device: observed-frame interpolation, then the trapezoid filter integrals over each band's
    !> support, compensated at FP32 kernel precision (or, with -Dcublas=true, a split-K matrix
    !> product over 256-pixel blocks, the blocks summed in FP64). The grid-only data for a
    !> redshift (observed-frame map, distance term, band weights and supports, trapezoid weights,
    !> IGM transmission) are cached in ctx%state%batch and on the device (prepare_batch_grid).
    !> Light-weighted ages and redshifted colors are not handled (callers use compute_magnitudes).

    use fsps_precision, only: WP, KP
    use fsps_constants, only: SAFE_FLOOR, C_LIGHT, ABS_MAG_ZEROPOINT_LOG
    use fsps_context_types, only: fsps_context_t, batch_grid_release
    use fsps_cosmology, only: compute_igm_transmission
    use fsps_interpolation, only: find_interval, interpolate_linear
    use fsps_special_functions, only: mag_from_flux
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
#ifdef FSPS_CUBLAS
    use cublas
    use openacc, only: acc_get_cuda_stream, acc_async_sync
#endif

    implicit none
    private

    public :: compute_magnitudes_columns, prepare_batch_grid
#ifdef FSPS_CUBLAS
    public :: blas, blas_init
#endif

    real(WP), parameter :: AB_ZEROPOINT = 48.60_wp
    integer, parameter  :: IDX_V_BAND = 1
    !> Pixels per block of the split-K band product (cuBLAS path): each block is a plain FP32 sum
    !> of at most 256 positive terms; the blocks are added in FP64
    integer, parameter :: BAND_BLOCK = 256
#ifdef FSPS_CUBLAS
    type(cublasHandle), save :: blas
    logical, save :: blas_ready = .false.
#endif

contains

    !> @brief Magnitudes of spec(:, 1:n_col) (present on the device, L_sun/Hz) at redshift zred,
    !> as compute_magnitudes(ctx, zred, spec(:, i), mags(:, i), mag_compute) for each column.
    subroutine compute_magnitudes_columns(ctx, zred, igm_factor, spec, n_col, mags, mag_compute)
        type(fsps_context_t), intent(inout) :: ctx
        real(WP), intent(in)                :: zred, igm_factor
        real(KP), intent(in)                :: spec(:,:)
        integer, intent(in)                 :: n_col
        real(WP), intent(out)               :: mags(:,:)
        integer, intent(in), optional       :: mag_compute(:)

        integer  :: nspec, nbands, b, m
        integer, allocatable  :: do_band(:)
        real(WP), allocatable :: fband(:,:)
        real(KP), allocatable :: fobs(:,:)
        real(WP) :: flux
        logical  :: done
#ifdef FSPS_CUBLAS
        integer :: npad
        integer, parameter :: GEMM_MIN_COLS = 16
#endif

        nspec = ctx%state%nspec
        nbands = ctx%state%nbands
        call prepare_batch_grid(ctx, zred, igm_factor)

        allocate(do_band(nbands), fband(nbands, n_col))
        do_band = 1
        if (present(mag_compute)) then
            b = min(nbands, size(mag_compute))
            do_band(1:b) = mag_compute(1:b)
            if (ctx%compute_vega_mags_val == 1) do_band(IDX_V_BAND) = 1
        end if

        associate(g => ctx%state%batch)
        done = .false.
#ifdef FSPS_CUBLAS
        ! cuBLAS for many spectra; a few (single CSP calls) take the compensated kernels below,
        ! which are more accurate (no plain FP32 block sums) and as fast there
        if (n_col >= GEMM_MIN_COLS) then
            npad = BAND_BLOCK * ((nspec + BAND_BLOCK - 1) / BAND_BLOCK)
            allocate(fobs(npad, n_col))
            !$acc data present(spec) create(fobs) copyout(fband)
            call observed_frame(n_col, nspec, npad, g%i_lo, g%t_obs, g%wl, spec, fobs)
            call band_fluxes_gemm(n_col, npad, nbands, g%cmat, fobs, fband)
            !$acc end data
            done = .true.
        end if
#endif
        if (.not. done) then
            allocate(fobs(n_col, nspec))
            !$acc data present(spec) copyin(do_band) create(fobs) copyout(fband)
            call observed_frame_t(n_col, nspec, g%i_lo, g%t_obs, g%wl, spec, fobs)
            call band_fluxes(n_col, nbands, g%bands_t, g%j_lo, g%j_hi, do_band, g%wtl, fobs, fband)
            !$acc end data
        end if
        end associate

        ! The cache's device copies are released with the context's (remove_from_device); a context
        ! not moved there by fsps_context_move_to_device would leave them behind, so drop them now
        if (.not. ctx%on_device) call batch_grid_release(ctx%state%batch)

        do m = 1, n_col
            mags(:, m) = ieee_value(1.0_wp, ieee_quiet_nan)
            do b = 1, nbands
                if (do_band(b) == 0) cycle
                flux = fband(b, m)
                mags(b, m) = mag_from_flux(flux) - AB_ZEROPOINT - 2.5_wp * ABS_MAG_ZEROPOINT_LOG + &
                             ctx%state%batch%dist_term
            end do
            if (ctx%compute_vega_mags_val == 1) then
                do b = 2, nbands
                    if (mags(IDX_V_BAND, m) == mags(IDX_V_BAND, m)) then
                        mags(b, m) = mags(b, m) - ctx%state%magvega(b) + ctx%state%magvega(IDX_V_BAND)
                    else
                        mags(b, m) = ieee_value(1.0_wp, ieee_quiet_nan)
                    end if
                end do
            end if
        end do
    end subroutine compute_magnitudes_columns

    !> Grid-only data for redshift zred (observed-frame map, distance term, band weights and
    !> supports, trapezoid weights, IGM transmission), cached in ctx%state%batch and on the device
    !> until zred or igm_factor change; released with the context's device data.
    subroutine prepare_batch_grid(ctx, zred, igm_factor)
        type(fsps_context_t), intent(inout) :: ctx
        real(WP), intent(in)                :: zred, igm_factor

        integer  :: nspec, nbands, j, b, ilo, ihi, npad
        real(WP) :: z_factor, target, dm, sw
        logical  :: igm

        nspec = ctx%state%nspec
        nbands = ctx%state%nbands
        igm = ctx%add_igm_absorption_val == 1 .and. zred > SAFE_FLOOR
        associate(g => ctx%state%batch)
        if (g%ready) then
            if (g%zred == zred .and. (g%igm .eqv. igm) .and. g%igm_factor == igm_factor .and. &
                size(g%wl) == nspec .and. size(g%j_lo) == nbands) return
        end if
        call batch_grid_release(g)
        allocate(g%i_lo(nspec), g%t_obs(nspec), g%wt(nspec), g%wtl(nspec), g%wl(nspec), g%igm_t(nspec), &
                 g%bands_t(nbands, nspec), g%j_lo(nbands), g%j_hi(nbands))
        g%zred = zred
        g%igm = igm
        g%igm_factor = igm_factor

        ! Trapezoid weights in frequency (dust emission) and wavelength (filters), formed in FP64
        do j = 1, nspec - 1
            g%wt(j) = real(0.5_wp * abs(C_LIGHT / ctx%state%spec_lambda(j+1) - C_LIGHT / ctx%state%spec_lambda(j)), KP)
            g%wtl(j) = real(0.5_wp * abs(ctx%state%spec_lambda(j + 1) - ctx%state%spec_lambda(j)), KP)
        end do
        g%wt(nspec) = 0.0_KP
        g%wtl(nspec) = 0.0_KP

        ! Observed-frame interpolation (prepare_observed_spectrum): per pixel, the rest-frame
        ! interval and position, and the distance term
        if (abs(zred) > SAFE_FLOOR) then
            z_factor = 1.0_wp + zred
            do j = 1, nspec
                target = ctx%state%spec_lambda(j) / z_factor
                g%i_lo(j) = max(1, min(find_interval(ctx%state%spec_lambda, target), nspec - 1))
                g%t_obs(j) = (target - ctx%state%spec_lambda(g%i_lo(j))) / &
                             (ctx%state%spec_lambda(g%i_lo(j)+1) - ctx%state%spec_lambda(g%i_lo(j)))
            end do
            dm = interpolate_linear(ctx%state%cosmospl(:,1), ctx%state%cosmospl(:,3), zred)
            if (dm /= dm .or. dm <= SAFE_FLOOR) then
                g%dist_term = 0.0_wp
            else
                g%dist_term = 5.0_wp * log10(dm / 10.0_wp) + mag_from_flux(z_factor)
            end if
        else
            g%i_lo = [(j, j = 1, nspec)]
            g%t_obs = 0.0_wp
            g%dist_term = 0.0_wp
        end if
        do j = 1, nspec
            if (ctx%state%spec_lambda(j) > tiny(0.0_wp)) then
                g%wl(j) = real(1.0_wp / ctx%state%spec_lambda(j), KP)
            else
                g%wl(j) = 0.0_KP
            end if
        end do

        ! Bands: transposed at kernel precision, and the support of each
        do b = 1, nbands
            g%bands_t(b, :) = real(ctx%state%bands(1:nspec, b), KP)
            ilo = nspec; ihi = 1
            do j = 1, nspec
                if (ctx%state%bands(j, b) /= 0.0_wp) then
                    ilo = min(ilo, j); ihi = max(ihi, j)
                end if
            end do
            ! Trapezoid intervals touching the support: j = ilo-1 .. ihi
            g%j_lo(b) = max(1, ilo - 1)
            g%j_hi(b) = min(nspec - 1, ihi)
            if (ilo > ihi) then
                g%j_lo(b) = 1; g%j_hi(b) = 0
            end if
        end do
        g%igm_t = 1.0_KP
        !$acc enter data copyin(g%i_lo, g%j_lo, g%j_hi, g%t_obs, g%wt, g%wtl, g%wl, g%bands_t, g%igm_t)
        if (igm) call compute_igm_transmission(ctx%state%spec_lambda, zred, igm_factor, g%igm_t)

#ifdef FSPS_CUBLAS
        ! Filter matrix of the split-K band product, padded to whole pixel blocks
        npad = BAND_BLOCK * ((nspec + BAND_BLOCK - 1) / BAND_BLOCK)
        allocate(g%cmat(npad, nbands))
        do b = 1, nbands
            do j = 1, npad
                if (j > nspec) then
                    g%cmat(j, b) = 0.0_KP
                else
                    sw = real(g%wtl(j), WP)
                    if (j > 1) sw = sw + real(g%wtl(j - 1), WP)
                    g%cmat(j, b) = real(sw * ctx%state%bands(j, b), KP)
                end if
            end do
        end do
        !$acc enter data copyin(g%cmat)
#endif
        g%ready = .true.
        end associate
    end subroutine prepare_batch_grid

    !> Observed-frame spectrum / lambda (prepare_observed_spectrum), zero in the padding rows.
    subroutine observed_frame(n_mod, nspec, npad, i_lo, t_obs, wl, total, fobs)
        integer, intent(in)     :: n_mod, nspec, npad, i_lo(:)
        real(WP), intent(in)    :: t_obs(:)
        real(KP), intent(in)    :: wl(:), total(:,:)
        real(KP), intent(out)   :: fobs(:,:)

        integer  :: m, j, il
        real(WP) :: v

        !$acc parallel loop gang vector collapse(2) present(i_lo, t_obs, wl, total, fobs) private(il, v)
        do m = 1, n_mod
            do j = 1, npad
                if (j <= nspec) then
                    il = i_lo(j)
                    if (t_obs(j) == 0.0_wp) then   ! rest frame (i_lo = j: no right neighbour at j = nspec)
                        v = real(total(il, m), WP)
                    else
                        v = real(total(il, m), WP) + t_obs(j) * (real(total(il + 1, m), WP) - real(total(il, m), WP))
                    end if
                    if (v /= v) v = 0.0_wp
                    fobs(j, m) = real(max(v, 0.0_wp), KP) * wl(j)
                else
                    fobs(j, m) = 0.0_KP
                end if
            end do
        end do
    end subroutine observed_frame

    !> observed_frame, written transposed (models, pixels) for band_fluxes.
    subroutine observed_frame_t(n_mod, nspec, i_lo, t_obs, wl, total, fobs_t)
        integer, intent(in)     :: n_mod, nspec, i_lo(:)
        real(WP), intent(in)    :: t_obs(:)
        real(KP), intent(in)    :: wl(:), total(:,:)
        real(KP), intent(out)   :: fobs_t(:,:)

        integer  :: m, j, il
        real(WP) :: v

        !$acc parallel loop gang vector collapse(2) present(i_lo, t_obs, wl, total, fobs_t) private(il, v)
        do j = 1, nspec
            do m = 1, n_mod
                il = i_lo(j)
                if (t_obs(j) == 0.0_wp) then
                    v = real(total(il, m), WP)
                else
                    v = real(total(il, m), WP) + t_obs(j) * (real(total(il + 1, m), WP) - real(total(il, m), WP))
                end if
                if (v /= v) v = 0.0_wp
                fobs_t(m, j) = real(max(v, 0.0_wp), KP) * wl(j)
            end do
        end do
    end subroutine observed_frame_t

#ifdef FSPS_CUBLAS
    !> Band fluxes as a split-K matrix product: per block of BAND_BLOCK pixels,
    !> partial = C^T . fobs with C(j, b) = (dl(j-1) + dl(j)) / 2 * band(j, b) (the trapezoid rule
    !> regrouped by pixel, formed in FP64; prepare_batch_grid); the blocks are summed in FP64.
    subroutine band_fluxes_gemm(n_mod, npad, nbands, cmat, fobs, fband)
        integer, intent(in)     :: n_mod, npad, nbands
        real(KP), intent(in)    :: cmat(:,:), fobs(:,:)
        real(WP), intent(out)   :: fband(:,:)

        real(KP), allocatable :: partial(:,:,:)
        integer  :: b, m, q, nkb, istat
        real(WP) :: sw

        nkb = npad / BAND_BLOCK
        allocate(partial(nbands, n_mod, nkb))
        call blas_init()
        !$acc data create(partial) present(cmat, fobs, fband)
        !$acc host_data use_device(cmat, fobs, partial)
        do q = 1, nkb
#ifdef KERNEL_SINGLE_PRECISION
            istat = cublasSgemm_v2(blas, CUBLAS_OP_T, CUBLAS_OP_N, nbands, n_mod, BAND_BLOCK, 1.0_KP, &
                                   cmat((q - 1) * BAND_BLOCK + 1, 1), npad, fobs((q - 1) * BAND_BLOCK + 1, 1), npad, &
                                   0.0_KP, partial(1, 1, q), nbands)
#else
            istat = cublasDgemm_v2(blas, CUBLAS_OP_T, CUBLAS_OP_N, nbands, n_mod, BAND_BLOCK, 1.0_KP, &
                                   cmat((q - 1) * BAND_BLOCK + 1, 1), npad, fobs((q - 1) * BAND_BLOCK + 1, 1), npad, &
                                   0.0_KP, partial(1, 1, q), nbands)
#endif
            if (istat /= 0) exit
        end do
        !$acc end host_data
        if (istat /= 0) error stop 'compute_csp_batch: cuBLAS GEMM failed'
        !$acc parallel loop gang vector collapse(2) present(partial, fband) private(sw, q)
        do m = 1, n_mod
            do b = 1, nbands
                sw = 0.0_wp
                !$acc loop seq
                do q = 1, nkb
                    sw = sw + real(partial(b, m, q), WP)
                end do
                fband(b, m) = sw
            end do
        end do
        !$acc end data
    end subroutine band_fluxes_gemm

    subroutine blas_init()
        integer :: istat
        if (blas_ready) return
        istat = cublasCreate(blas)
        istat = cublasSetStream(blas, acc_get_cuda_stream(acc_async_sync))
        blas_ready = .true.
    end subroutine blas_init
#endif

    !> Band fluxes: the trapezoid filter integral of fobs over each band's support, compensated at
    !> FP32 kernel precision.
    subroutine band_fluxes(n_mod, nbands, bands_t, j_lo, j_hi, do_band, wt, fobs_t, fband)
        integer, intent(in)     :: n_mod, nbands, j_lo(:), j_hi(:), do_band(:)
        !> Observed spectrum / lambda, transposed (models, pixels): the threads of a gang read
        !> consecutive models
        real(KP), intent(in)    :: bands_t(:,:), wt(:), fobs_t(:,:)
        real(WP), intent(out)   :: fband(:,:)

        !> Models per gang: the threads of a gang share one band, hence one trip count (no
        !> divergence)
        integer, parameter :: MBLK = 128, SMALL_N = 16
        integer  :: m, j, b, mb, ml
        real(KP) :: s, c, term, tsum

#ifdef OMP_TARGET
        ! OpenMP target offload (-Daccel=omp_gpu): consecutive threads take consecutive models of
        ! one band (coalesced reads of fobs_t, one trip count), compensated at FP32
        !$omp target teams distribute parallel do collapse(2) private(s, c, term, tsum, j) &
        !$omp map(to: bands_t, j_lo, j_hi, do_band, wt, fobs_t) map(from: fband)
        do b = 1, nbands
            do m = 1, n_mod
                s = 0.0_KP
                c = 0.0_KP
                if (do_band(b) /= 0) then
                    do j = j_lo(b), j_hi(b)
                        if (KP == WP) then
                            s = s + wt(j) * (fobs_t(m, j) * bands_t(b, j) + fobs_t(m, j + 1) * bands_t(b, j + 1))
                        else
                            term = wt(j) * (fobs_t(m, j) * bands_t(b, j) + fobs_t(m, j + 1) * bands_t(b, j + 1)) - c
                            tsum = s + term
                            c = (tsum - s) - term
                            s = tsum
                        end if
                    end do
                end if
                fband(b, m) = real(s, WP)
            end do
        end do
        return
#endif
#ifdef _OPENACC
        if (n_mod < SMALL_N) then
            ! Few spectra: a gang per (band, spectrum), the pixels reduced across its threads
            ! (a tree reduction, pairwise like the dust integrals in fsps_dust). Measured on the
            ! P100: 1 spectrum 0.65 -> 0.03 ms; slower than the kernel below from ~100 spectra
            !$acc parallel loop gang collapse(2) present(bands_t, j_lo, j_hi, do_band, wt, fobs_t, fband) private(s)
            do b = 1, nbands
                do m = 1, n_mod
                    s = 0.0_KP
                    if (do_band(b) /= 0) then
                        !$acc loop vector reduction(+:s)
                        do j = j_lo(b), j_hi(b)
                            s = s + wt(j) * (fobs_t(m, j) * bands_t(b, j) + fobs_t(m, j + 1) * bands_t(b, j + 1))
                        end do
                    end if
                    fband(b, m) = real(s, WP)
                end do
            end do
            return
        end if
#endif
        !$acc parallel loop gang collapse(2) vector_length(MBLK) present(bands_t, j_lo, j_hi, do_band, wt, fobs_t, fband)
        do b = 1, nbands
            do mb = 1, (n_mod + MBLK - 1) / MBLK
                !$acc loop vector private(m, s, c, term, tsum, j)
                do ml = 1, MBLK
                    m = (mb - 1) * MBLK + ml
                    if (m <= n_mod) then
                        s = 0.0_KP
                        c = 0.0_KP
                        if (do_band(b) /= 0) then
                            !$acc loop seq
                            do j = j_lo(b), j_hi(b)
                                if (KP == WP) then
                                    s = s + wt(j) * (fobs_t(m, j) * bands_t(b, j) + fobs_t(m, j + 1) * bands_t(b, j + 1))
                                else
                                    term = wt(j) * (fobs_t(m, j) * bands_t(b, j) + fobs_t(m, j + 1) * bands_t(b, j + 1)) - c
                                    tsum = s + term
                                    c = (tsum - s) - term
                                    s = tsum
                                end if
                            end do
                        end if
                        fband(b, m) = real(s, WP)
                    end if
                end do
            end do
        end do
    end subroutine band_fluxes

end module fsps_batch_photometry
